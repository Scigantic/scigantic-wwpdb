"""Find PDB ENTRIES by structured criteria, then their ligands.

The rest of this package answers "what is component ATP?". This module answers the
question that comes first: "which structures should I even be looking at?" — e.g.
every X-ray protein-ligand complex at 2.5 A or better released after some date,
which is how you build a holdout set that a model cannot have trained on.

  search_entries(released_after="2025-06-30", resolution_max=2.5)  -> ["10AF", ...]
  count_entries(released_after="2025-06-30")                      -> 8765
  entry_ligands(["10AF"])  -> {"10AF": ["ADP", "MG"]}

Pipe the ligands into find() for names/SMILES, or into is_interesting_ligand() to
drop the ions and cryoprotectants that most entries carry alongside the real one.

Two things to know about RCSB's entry search, both of which have bitten us:

* A query with NO hits returns HTTP 204 and an empty body, not an empty result
  set. _http.post_json already maps that to {}, so a zero here is indistinguishable
  from a malformed query that the server rejected as unsatisfiable. ALWAYS sanity
  check a surprising zero by relaxing the most specific filter and confirming the
  count becomes non-zero; count_entries exists partly to make that cheap.
* `rows` is capped per request, so anything past the first page needs explicit
  pagination. search_entries pages for you up to `limit`; it does not silently
  truncate at RCSB's page size, which is the failure mode of a hand-rolled query.
"""
from __future__ import annotations

import json
import os
import re

from . import _http
from .rcsb import SEARCH_URL, find

# Entries live on the same GraphQL endpoint as components, but are kept on their own
# name so one can be repointed (a mirror, a test double) without the other.
DATA_GRAPHQL_ENTRIES = os.environ.get(
    "SCIGANTIC_PDB_GRAPHQL", "https://data.rcsb.org/graphql")

# RCSB rejects larger pages; also the point past which one response gets unwieldy.
_PAGE = 1000

# Ligand ids that appear in thousands of entries as buffers, ions, cryoprotectants
# and detergents. Not chemistry anyone is trying to predict, and they pollute any
# "protein-ligand complex" set you assemble by entry-level filters alone.
COMMON_ADDITIVES = frozenset("""
HOH DOD NA K MG CA MN FE FE2 ZN CU CU1 NI CO CD HG PT AU AG CL BR IOD F
SO4 PO4 NO3 CO3 ACT ACY EDO PEG PG4 PGE 1PE 2PE P6G 7PE XPE MPD GOL DMS
TRS EPE MES HEPES BTB BIS IMD TAM CIT FLC TLA MLA MLI SIN SCN AZI CYN
NH4 NH2 OH O OXY PER H2S BME DTT TCE GSH FMT EOH IPA ACN MOH CCN DMF
BOG LDA C8E OCT HEX HEZ PE4 PE5 BU1 BU2 BU3 MRD PDO TBU URE GLC NAG NDG
BMA MAN FUC GAL SIA XYS BGC SUC TRE EPD SRT PYR OAA AKG
""".split())


def _nodes(released_after, released_before, resolution_max, method,
           with_ligand, contains_ccd, sequence, identity_cutoff, extra):
    out = []
    if released_after:
        out.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date",
            "operator": "greater", "value": str(released_after)}})
    if released_before:
        out.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date",
            "operator": "less", "value": str(released_before)}})
    if resolution_max is not None:
        out.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": float(resolution_max)}})
    if method:
        out.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.experimental_method",
            "operator": "exact_match", "value": str(method)}})
    if with_ligand is not None:
        # "has a small molecule" and "has a protein" are separate counts; asking for
        # a protein-ligand complex means constraining both, or you also collect
        # nucleic-acid-only and protein-only entries.
        op, val = ("greater", 0) if with_ligand else ("equals", 0)
        out.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.nonpolymer_entity_count",
            "operator": op, "value": val}})
        if with_ligand:
            out.append({"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                "operator": "greater", "value": 0}})
    if contains_ccd:
        ids = [contains_ccd] if isinstance(contains_ccd, str) else list(contains_ccd)
        out.append({"type": "terminal", "service": "text_chem", "parameters": {
            "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
            "operator": "in", "value": [str(i).upper() for i in ids]}})
    if sequence:
        out.append({"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 1, "identity_cutoff": float(identity_cutoff),
            "sequence_type": "protein", "value": str(sequence)}})
    if extra:
        out.extend(extra)
    # `method` alone is not a filter in any useful sense — it still selects ~200,000
    # entries — so it does not count towards "you asked for something specific".
    selective = [n for n in out if n["parameters"].get("attribute")
                 != "rcsb_entry_info.experimental_method"]
    if not selective:
        raise ValueError(
            "refusing to ask RCSB for every entry it has: pass at least one filter "
            "besides method (a date, a resolution, with_ligand, a sequence, ...)")
    return out


def _payload(nodes, start, rows):
    query = nodes[0] if len(nodes) == 1 else {
        "type": "group", "logical_operator": "and", "nodes": nodes}
    return {
        "query": query,
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": int(start), "rows": int(rows)},
            "results_content_type": ["experimental"],
        },
    }


def count_entries(*, released_after=None, released_before=None, resolution_max=None,
                  method="X-ray", with_ligand=None, contains_ccd=None,
                  sequence=None, identity_cutoff=0.9, extra=None) -> int:
    """How many entries match, without pulling the ids. Cheap enough to use as the
    control for a surprising result: relax one filter, confirm the count moves."""
    nodes = _nodes(released_after, released_before, resolution_max, method,
                   with_ligand, contains_ccd, sequence, identity_cutoff, extra)
    data = _http.post_json(SEARCH_URL, _payload(nodes, 0, 1))
    return int(data.get("total_count") or 0)


def search_entries(*, released_after=None, released_before=None, resolution_max=None,
                   method="X-ray", with_ligand=None, contains_ccd=None,
                   sequence=None, identity_cutoff=0.9, limit=100, extra=None) -> list:
    """Entry ids matching the given filters, paginated for you.

    Dates are ISO strings ("2025-06-30") and both bounds are EXCLUSIVE, matching
    RCSB's greater/less operators. `method` defaults to X-ray because resolution
    only means the same thing within one method; pass method=None to drop it.
    `with_ligand=True` additionally requires a protein chain, so you get complexes
    rather than every entry that happens to carry an ion. `extra` takes raw RCSB
    query nodes for anything not wrapped here.

    Returns [] for no hits — see this module's note on HTTP 204, and use
    count_entries to tell "genuinely none" from "query the server would not accept".
    """
    nodes = _nodes(released_after, released_before, resolution_max, method,
                   with_ligand, contains_ccd, sequence, identity_cutoff, extra)
    want = int(limit)
    ids: list[str] = []
    start = 0
    while len(ids) < want:
        rows = min(_PAGE, want - len(ids))
        data = _http.post_json(SEARCH_URL, _payload(nodes, start, rows))
        batch = [r["identifier"] for r in data.get("result_set") or []]
        if not batch:
            break
        ids.extend(batch)
        start += len(batch)
        if len(batch) < rows:
            break
    return ids[:want]


_ENTRY_GQL = """{ entries(entry_ids: %s) {
  rcsb_id
  nonpolymer_entities { nonpolymer_comp { chem_comp { id } } }
} }"""


def entry_ligands(ids, additives=False) -> dict:
    """{entry_id: [ccd_id, ...]} for many entries in ONE GraphQL request.

    By default the ions, buffers and cryoprotectants in COMMON_ADDITIVES are
    stripped, since almost every entry carries some and they are never the ligand
    you meant. Pass additives=True to get the raw list. Entries RCSB does not
    return, and entries whose only ligands were additives, map to [].
    """
    if isinstance(ids, str):
        ids = [ids]
    ids = [str(i).strip().upper() for i in ids if str(i).strip()]
    if not ids:
        return {}
    out: dict[str, list] = {i: [] for i in ids}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        data = _http.post_json(DATA_GRAPHQL_ENTRIES,
                               {"query": _ENTRY_GQL % json.dumps(chunk)})
        for row in (data.get("data", {}) or {}).get("entries") or []:
            if not row:
                continue
            key = str(row.get("rcsb_id") or "").strip().upper()
            got = []
            for ne in row.get("nonpolymer_entities") or []:
                cc = ((ne or {}).get("nonpolymer_comp") or {}).get("chem_comp") or {}
                cid = str(cc.get("id") or "").strip().upper()
                if cid and (additives or cid not in COMMON_ADDITIVES):
                    got.append(cid)
            if key in out:
                out[key] = sorted(set(got))
    return out


def is_interesting_ligand(ccd_id, min_heavy_atoms=10) -> bool:
    """Is this CCD plausibly the ligand of interest rather than an additive?

    Cheap two-part test: not in COMMON_ADDITIVES, and big enough to matter. Heavy
    atoms are counted from the formula string, so this needs no CIF fetch; it does
    one find() call, which is cached.
    """
    cid = str(ccd_id).strip().upper()
    if not cid or cid in COMMON_ADDITIVES:
        return False
    hits = find([cid])
    if not hits:
        return False
    return _heavy_atoms(hits[0].formula) >= int(min_heavy_atoms)


def _heavy_atoms(formula) -> int:
    """Heavy (non-hydrogen) atom count from a CCD formula like 'C10 H16 N5 O13 P3'.

    Element symbols in CCD formulas are proper-cased ('Mg', 'Fe', 'O4 S'), which is
    what lets a two-letter element be told apart from two one-letter ones. Do not
    pass an upper-cased formula: 'MG' parses as M + G, i.e. two atoms. The CCD *id*
    is upper-case ('MG') and the formula is not, so don't swap them.
    """
    if not formula:
        return 0
    n = 0
    for el, cnt in re.findall(r"([A-Z][a-z]?)\s*(\d*)", str(formula)):
        if not el or el == "H" or el == "D":
            continue
        n += int(cnt) if cnt else 1
    return n

