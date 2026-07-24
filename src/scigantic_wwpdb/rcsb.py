"""Explore the CCD without downloading it: RCSB's search + batch data APIs.

  search("heme")     -> component ids by name        (RCSB Search API)
  find(ids)          -> metadata for MANY in one GET  (RCSB Data GraphQL)

These are the fast skim layer. Pull full structure (atoms/coords/bonds) only for
what you actually want, with component()/components() in .structure.
"""
from __future__ import annotations
import os
import json

from . import _http
from .model import Summary

DATA_GRAPHQL = os.environ.get("SCIGANTIC_CCD_GRAPHQL", "https://data.rcsb.org/graphql")
SEARCH_URL = os.environ.get(
    "SCIGANTIC_CCD_SEARCH", "https://search.rcsb.org/rcsbsearch/v2/query")

# One GraphQL query, many components. RCSB returns null for unknown ids, which we
# drop, so find(["ATP", "NOPE"]) is safe.
_GQL = """{ chem_comps(comp_ids: %s) {
  rcsb_id
  chem_comp { id name formula type formula_weight }
  rcsb_chem_comp_descriptor { SMILES_stereo InChI InChIKey }
} }"""


def find(ids) -> list:
    """Metadata for one or many components in a SINGLE request: a list of
    Summary(id, name, formula, type, weight, smiles, inchikey). No CIF fetch, no
    atoms/coords. This is how you skim a whole result set fast. Unknown ids are
    simply absent from the result."""
    if isinstance(ids, str):
        ids = [ids]
    ids = [str(i).strip().upper() for i in ids if str(i).strip()]
    if not ids:
        return []
    data = _http.post_json(DATA_GRAPHQL, {"query": _GQL % json.dumps(ids)})
    out = []
    for row in (data.get("data", {}) or {}).get("chem_comps") or []:
        if not row:
            continue
        cc = row.get("chem_comp") or {}
        desc = row.get("rcsb_chem_comp_descriptor") or {}
        out.append(Summary(
            id=row.get("rcsb_id") or cc.get("id"),
            name=cc.get("name"), formula=cc.get("formula"),
            type=cc.get("type"), weight=cc.get("formula_weight"),
            smiles=desc.get("SMILES_stereo"),
            inchi=desc.get("InChI"), inchikey=desc.get("InChIKey"),
        ))
    return out


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
