# Changelog

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
