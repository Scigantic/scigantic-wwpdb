"""Full structure: fetch a component's mmCIF and parse it (atoms, coords, bonds).

component(id) for one; components(ids) fetches many in parallel. Reach for these
once search()/find() have told you which components you actually want.
"""
from __future__ import annotations
import os

from . import _http
from ._cif import parse
from .model import Component

# Per-component MODERN mmCIF (atoms with ideal coords + SMILES/InChI) — RCSB
# serves one file per component, addressable by id, no directory listing.
LIGAND_BASE = os.environ.get(
    "SCIGANTIC_CCD_LIGAND_BASE", "https://files.rcsb.org/ligands/download")


def component_url(ccd_id: str) -> str:
    return f"{LIGAND_BASE}/{str(ccd_id).strip().upper()}.cif"


def fetch_cif(ccd_id: str, retries: int = 2) -> str:
    """Fetch one component's raw mmCIF text (a few KB)."""
    url = component_url(ccd_id)
    last = None
    for _ in range(retries + 1):
        try:
            r = _http.get(url)
            if r.status_code == 404:
                raise KeyError(f"CCD component {ccd_id!r} not found ({url})")
            r.raise_for_status()
            return r.text
        except KeyError:
            raise
        except Exception as exc:  # transient network / 5xx — retry
            last = exc
    raise last


def component(ccd_id: str) -> Component:
    """Fetch and parse one CCD component by id (full structure)."""
    comp = parse(fetch_cif(ccd_id))
    if not comp.id:
        comp.id = str(ccd_id).strip().upper()
    return comp


def components(ids, workers: int = 8, errors: str = "omit") -> list:
    """Fetch+parse full components for many ids in parallel (order preserved).
    Much faster than a loop for a search-result set.

    A single bad id must not sink the batch. errors="omit" (default) puts None in
    place of any id that fails (404 / network), so you can filter with
    `[c for c in components(ids) if c]`; errors="raise" restores fail-fast."""
    if errors not in ("omit", "raise"):
        raise ValueError("errors must be 'omit' or 'raise'")
    if errors == "raise":
        return _http.pmap(component, ids, workers=workers)

    def _safe(cid):
        try:
            return component(cid)
        except Exception:
            return None  # bad/unknown id or transient failure — don't kill the batch
    return _http.pmap(_safe, ids, workers=workers)
