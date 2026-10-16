# Changelog

## 0.3.0

- New `entries.py`: find PDB **entries** by structured criteria, which is the
  question that comes before "what is component ATP?".
  - `search_entries(released_after=, released_before=, resolution_max=, method=,
    with_ligand=, contains_ccd=, sequence=, identity_cutoff=, limit=, extra=)`
    returns entry ids, paginating past RCSB's per-request `rows` cap rather than
    silently stopping at the first page.
  - `count_entries(...)` returns the same query's total without pulling ids, so
    it is cheap to use as the control for a surprising zero.
  - `entry_ligands(ids)` maps many entries to their non-polymer CCD ids in one
    GraphQL request, stripping the ions/buffers/cryoprotectants in the new
    `COMMON_ADDITIVES` set (pass `additives=True` to keep them).
  - `is_interesting_ligand(ccd_id, min_heavy_atoms=10)` combines that set with a
    heavy-atom count read off the formula, so no CIF fetch is needed.
- `search_entries()`/`count_entries()` refuse a query whose only filter is
  `method`, since that still selects ~200,000 entries.
- Documented two RCSB behaviours that make hand-rolling these queries risky: a
  query with no hits answers HTTP 204 with an empty body (so a zero is
  indistinguishable from a query the server would not accept), and results past
  the first page need explicit pagination.

## 0.2.1

- Pinned dependency ceilings: `gemmi>=0.6,<1`, `requests>=2.25,<3` (previously
  unbounded on both).
- CI now also installs without the `rdkit` extra and runs the suite, so the
  "rdkit absent" code path (`chem.py`'s lazy import, already guarded in tests
  with `pytest.importorskip("rdkit")`) is actually exercised rather than
  always bundled in via the `test` extra.
- Added README badges (CI, PyPI version, license, Python versions).
- Added this changelog.

## 0.2.0

- `find()`/`component()` responses are now cached per id on disk, on by
  default (14-day TTL). `ccd.disable_cache()` / `ccd.enable_cache(cache_dir=,
  ttl_days=)` control it.
- Version now single-sourced from `src/scigantic_wwpdb/_version.py` via
  hatchling's dynamic version (previously hardcoded in `pyproject.toml`).

## 0.1.1

- `components()` tolerates per-id failures instead of failing the whole batch.
- Kept the single-stream `load_dictionary` path.
- Added CI publish-to-PyPI workflow (token via repo secret `PYPI_KEY`).
- Lifted the HTTP connection pool to 32 so `components(workers=N)` actually
  parallelizes.

## 0.1.0

First release.

- Restructured into modules: `rcsb` (search + batch metadata), `structure`
  (fetch + parse), `_cif` (gemmi), `chem` (RDKit/SDF), `dictionary` (the
  bundle), `model` (records).
- RCSB search and batch-metadata support for fast exploration (`search`,
  `find`).
