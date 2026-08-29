"""Test isolation for the on-by-default response cache (see cache.py).

Without this, a test's mocked find()/fetch_cif() response would be written
to the real ~/.cache/scigantic-wwpdb, and a later test (in this run, or a
future dev's next `pytest` invocation) could silently read that stale
mocked value back instead of exercising its own mock. Every test gets its
own empty cache directory instead.
"""
import pytest

import scigantic_wwpdb as ccd


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path):
    ccd.enable_cache(cache_dir=str(tmp_path))
