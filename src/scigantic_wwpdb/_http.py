"""One shared HTTP session, plus a tiny parallel map for fetch-many."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor

import requests
from requests.adapters import HTTPAdapter

from ._version import __version__

_UA = {"User-Agent": f"scigantic-wwpdb/{__version__} "
                     "(+https://scigantic.com; mailto:support@scigantic.com)"}

session = requests.Session()
session.headers.update(_UA)
# Keep-alive is on by default; lift the connection pool past requests' default
# of 10 so components(ids, workers=N) actually runs N in parallel instead of
# queueing. Benchmarked from a us-east-1 pod: per-component CIF GETs plateau
# around 8 concurrent (files.rcsb.org is CloudFront-edge-cached, ~20 ms warm),
# so 32 leaves generous headroom without hammering the origin.
_adapter = HTTPAdapter(pool_connections=8, pool_maxsize=32, max_retries=0)
session.mount("https://", _adapter)
session.mount("http://", _adapter)


def get(url, **kw):
    kw.setdefault("timeout", 30)
    return session.get(url, **kw)


def post_json(url, payload, timeout=30):
    """POST JSON, return parsed JSON. Empty/204 bodies (e.g. a search with no
    hits) come back as {} rather than raising."""
    r = session.post(url, json=payload, timeout=timeout)
    if r.status_code == 204 or not r.content:
        return {}
    r.raise_for_status()
    return r.json()


def pmap(fn, items, workers=8):
    """Map fn over items concurrently (I/O-bound: many small HTTP GETs). Order
    preserved. Serial for 0/1 items so single lookups skip the pool."""
    items = list(items)
    if len(items) <= 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        return list(ex.map(fn, items))
