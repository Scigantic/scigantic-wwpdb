"""scigantic_wwpdb — explore the wwPDB Chemical Component Dictionary (CCD) from Python.

The CCD is the reference chemistry for every small molecule, ligand, ion, and
modified residue in the PDB (~45,000 components) — the chemistry AlphaFold 3 uses
for ligand pose prediction, and what you reach for whenever a structure has a
HETATM to reason about.

Three tiers, so exploring stays fast and nothing downloads until you ask:

  search("heme")            -> component ids by name          (RCSB text search)
  find(["ATP", "HEM"])      -> metadata for MANY in one call  (RCSB data GraphQL)
  component("ATP")          -> full structure: atoms + ideal 3-D coords + bonds
  components(["ATP","ADP"]) -> full structures, fetched in parallel

Then, per component:
  smiles / inchi / inchikey / formula / name (id)   -> one-liners (fast metadata)
  to_rdkit / to_sdf (id)                            -> RDKit Mol / SDF (built from
                                                       atoms+bonds, so heme works)
  load_dictionary()                                 -> the whole CCD, one stream

The idiom for browsing is `find(search("kinase inhibitor"))`: search names it,
find skims the set in a single request, and you only pull full structures for the
handful you actually want. Every function is a few lines you can read and bend.

One tier above that, for when you do not yet know which structures matter:

  search_entries(released_after="2025-06-30", resolution_max=2.5, with_ligand=True)
  count_entries(...)                -> the same query's total, without the ids
  entry_ligands(ids)                -> {entry: [ccd, ...]}, additives stripped

That is how you assemble a holdout set a model cannot have trained on: filter
entries by release date, then keep the ones whose ligand is real chemistry rather
than a buffer. See entries.py for the two RCSB behaviours that make hand-rolling
this risky — a no-hit query returns HTTP 204, and results past the first page need
explicit pagination.

find()/component() responses are cached on disk per id, on by default (14-day
TTL): disable_cache() turns it off, enable_cache(ttl_days=...) changes the
directory or expiry. See cache.py for why this defaults on, unlike
scigantic-chembl's/scigantic-bindingdb's opt-in caching.
"""
from ._version import __version__
from .model import Summary, Component, Atom, Bond
from ._cif import parse, cif_block
from .rcsb import search, find, name, formula, smiles, inchi, inchikey
from .entries import (search_entries, count_entries, entry_ligands,
                      is_interesting_ligand, COMMON_ADDITIVES)
from .structure import component, components, component_url, fetch_cif
from .chem import to_rdkit, to_sdf
from .dictionary import load_dictionary, read_dictionary
from .cache import enable_cache, disable_cache, is_cache_enabled, cache_dir, clear as clear_cache

__all__ = [
    "__version__",
    # explore
    "search", "find",
    # which structures to look at in the first place
    "search_entries", "count_entries", "entry_ligands", "is_interesting_ligand",
    "COMMON_ADDITIVES",
    # full structure
    "component", "components", "parse", "component_url", "fetch_cif",
    # per-component one-liners
    "name", "formula", "smiles", "inchi", "inchikey",
    # chem
    "to_rdkit", "to_sdf",
    # bulk
    "load_dictionary", "read_dictionary",
    # caching (on by default, see cache.py)
    "enable_cache", "disable_cache", "is_cache_enabled", "cache_dir", "clear_cache",
    # models / low-level
    "Summary", "Component", "Atom", "Bond", "cif_block",
]
