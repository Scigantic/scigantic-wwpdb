"""Explore the CCD without downloading it: RCSB's search + batch data APIs.

  search("heme")     -> component ids by name        (RCSB Search API)
  find(ids)          -> metadata for MANY in one GET  (RCSB Data GraphQL)

These are the fast skim layer. Pull full structure (atoms/coords/bonds) only for
what you actually want, with component()/components() in .structure.
"""
from __future__ import annotations
import os
import json
from dataclasses import asdict

from . import _http, cache
from .model import Summary

DATA_GRAPHQL = os.environ.get("SCIGANTIC_CCD_GRAPHQL", "https://data.rcsb.org/graphql")
SEARCH_URL = os.environ.get(
    "SCIGANTIC_CCD_SEARCH", "https://search.rcsb.org/rcsbsearch/v2/query")

# RCSB's server-side cap on comp_ids per chem_comps query. Measured: 1000 is
# accepted, 1100 is not ("Input 'comp_ids' exceeds maximum allowed size of 1000").
MAX_IDS_PER_QUERY = 1000

# One GraphQL query, many components. RCSB returns null for unknown ids, which we
# drop, so find(["ATP", "NOPE"]) is safe.
_GQL = """{ chem_comps(comp_ids: %s) {
  rcsb_id
  chem_comp { id name formula type formula_weight }
  rcsb_chem_comp_descriptor { SMILES_stereo InChI InChIKey }
} }"""


def find(ids) -> list:
    """Metadata for one or many components, a request per 1000: a list of
    Summary(id, name, formula, type, weight, smiles, inchikey). No CIF fetch, no
    atoms/coords. This is how you skim a whole result set fast. Unknown ids are
    simply absent from the result.

    RCSB caps `comp_ids` at MAX_IDS_PER_QUERY (1000). Over that it answers HTTP 200
    with a GraphQL `errors` array and `chem_comps: null`, which used to make this
    function return [] — a silent empty result that reads exactly like "none of
    those exist". Ids are chunked to stay under the cap, and a GraphQL error that
    yields no rows now raises instead of being swallowed.

    Cached per id (see cache.py), on by default. name()/formula()/smiles()/
    inchi()/inchikey() each call this via _one(), so looking up several of
    those for the same id used to mean several identical GraphQL round trips;
    now only the first one for a given id (per cache.py's ttl_days, 14 by
    default) touches the network. Only ids RCSB actually returned a row for
    are cached; an id that comes back empty (unknown/typo) is re-queried
    every time rather than caching a negative result.
    """
    if isinstance(ids, str):
        ids = [ids]
    ids = [str(i).strip().upper() for i in ids if str(i).strip()]
    if not ids:
        return []

    hits: dict[str, Summary] = {}
    missing = []
    for i in ids:
        cached = cache.get(f"find:{i}")
        if cached is not None:
            hits[i] = Summary(**cached)
        else:
            missing.append(i)

    for chunk in (missing[i:i + MAX_IDS_PER_QUERY]
                  for i in range(0, len(missing), MAX_IDS_PER_QUERY)):
        data = _http.post_json(DATA_GRAPHQL, {"query": _GQL % json.dumps(chunk)})
        rows = (data.get("data", {}) or {}).get("chem_comps") or []
        if not rows and data.get("errors"):
            raise RuntimeError(
                "RCSB GraphQL rejected a chem_comps query for "
                f"{len(chunk)} ids: {data['errors'][0].get('message')}")
        for row in rows:
            if not row:
                continue
            cc = row.get("chem_comp") or {}
            desc = row.get("rcsb_chem_comp_descriptor") or {}
            s = Summary(
                id=str(row.get("rcsb_id") or cc.get("id") or ""),
                name=cc.get("name"), formula=cc.get("formula"),
                type=cc.get("type"), weight=cc.get("formula_weight"),
                smiles=desc.get("SMILES_stereo"),
                inchi=desc.get("InChI"), inchikey=desc.get("InChIKey"),
            )
            key = str(s.id).strip().upper()
            hits[key] = s
            cache.put(f"find:{key}", asdict(s))

    return [hits[i] for i in ids if i in hits]


def search(query: str, limit: int = 25) -> list:
    """Component ids whose name matches `query` (word match, case-insensitive),
    via RCSB's text search — no local catalog to build or ship. Pipe into
    find() to see them: `find(search("kinase inhibitor"))`."""
    payload = {
        "query": {
            "type": "terminal", "service": "text_chem",
            "parameters": {"attribute": "chem_comp.name",
                           "operator": "contains_words", "value": str(query)},
        },
        "return_type": "mol_definition",
        "request_options": {"paginate": {"start": 0, "rows": int(limit)}},
    }
    try:
        data = _http.post_json(SEARCH_URL, payload)
    except Exception:
        return []
    return [r["identifier"] for r in data.get("result_set", [])][:limit]


# ── one-liners (metadata, so they use the fast data API, not the full CIF) ──
def _one(ccd_id):
    hits = find([ccd_id])
    return hits[0] if hits else None


def name(ccd_id):
    s = _one(ccd_id)
    return s.name if s else None


def formula(ccd_id):
    s = _one(ccd_id)
    return s.formula if s else None


def smiles(ccd_id):
    s = _one(ccd_id)
    return s.smiles if s else None


def inchi(ccd_id):
    s = _one(ccd_id)
    return s.inchi if s else None


def inchikey(ccd_id):
    s = _one(ccd_id)
    return s.inchikey if s else None
