"""Local response cache, ON by default.

Unlike scigantic-chembl's and scigantic-bindingdb's caching (both default
OFF): those two read from a public S3 mirror with no meaningful rate limit,
so caching there is a pure convenience, opt-in because zero-setup-by-default
is the point. This package calls RCSB's live search/data/ligand-download
APIs for every lookup, the same category as scigantic-pubchem's PUG REST
calls, not a mirror. Concretely: name("HEM"), formula("HEM") and
smiles("HEM") each independently called find(["HEM"]) and re-fetched
identical metadata three times before this existed. Caching stays on unless
turned off explicitly, matching scigantic-pubchem's default.

Cached entries expire after ttl_days (14 by default). The CCD gets new
component ids fairly often as structures deposit, so this defaults shorter
than scigantic-pubchem's 30 days: long enough that a notebook session or a
same-day analysis re-using the same ids stays fast, short enough that a
newly-deposited id doesn't stay "not found" in a long-running process.
Pass ttl_days=None to disable expiry entirely.

Reads and writes of individual entries are safe to call concurrently,
including two threads racing to fill the same key: each write goes to its
own uniquely-named temp file before an atomic os.replace() into place, so a
reader never observes a partial write and two concurrent writers of the
same key never collide on the same temp path.
"""

from __future__ import annotations

import json
import os
import sys
import time
import uuid
from pathlib import Path
from typing import Any

_enabled = True
_cache_dir: Path | None = None
_ttl_seconds: float | None = 14 * 86400


def _default_cache_dir() -> Path:
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = str(Path.home() / "Library" / "Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "scigantic-wwpdb"


def _resolve_dir() -> Path:
    global _cache_dir
    if _cache_dir is None:
        env = os.environ.get("SCIGANTIC_WWPDB_CACHE")
        _cache_dir = Path(env) if env else _default_cache_dir()
        _cache_dir.mkdir(parents=True, exist_ok=True)
    return _cache_dir


def enable_cache(cache_dir: str | None = None, ttl_days: float | None = 14) -> Path:
    """Turn caching on (it already is, by default) and optionally point it at
    a specific directory and/or change how long an entry stays valid.

    ttl_days=None disables expiry: entries are reused forever until
    clear_cache() or disable_cache(). Returns the resolved directory.
    """
    global _enabled, _cache_dir, _ttl_seconds
    if cache_dir is not None:
        _cache_dir = Path(cache_dir)
        _cache_dir.mkdir(parents=True, exist_ok=True)
    else:
        _resolve_dir()
    _ttl_seconds = ttl_days * 86400 if ttl_days is not None else None
    _enabled = True
    return _cache_dir  # type: ignore[return-value]


def disable_cache() -> None:
    """Turn caching off. Every call hits RCSB fresh until re-enabled."""
    global _enabled
    _enabled = False


def is_cache_enabled() -> bool:
    return _enabled


def cache_dir() -> Path:
    return _resolve_dir()


def _path(key: str) -> Path:
    # Keys here are short ("cif:HEM", "find:ATP"): no hashing needed, and a
    # readable filename makes the cache directory easy to inspect by hand.
    # ":" is replaced too, not just "/" — it's a legal POSIX filename byte
    # but reserved on Windows (drive letters), and this package supports it
    # (see _default_cache_dir's win32 branch).
    safe = key.replace("/", "_").replace(":", "_")
    return _resolve_dir() / f"{safe}.json"


def get(key: str) -> Any | None:
    if not _enabled:
        return None
    file = _path(key)
    if not file.exists():
        return None
    try:
        entry = json.loads(file.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if _ttl_seconds is not None and time.time() - entry.get("cached_at", 0) > _ttl_seconds:
        file.unlink(missing_ok=True)
        return None
    return entry["value"]


def put(key: str, value: Any) -> None:
    if not _enabled:
        return
    file = _path(key)
    # Unique per call, not just per key: two threads racing to fill the same
    # key (components(ids, workers=N) hitting the same id twice, or two
    # notebooks sharing a cache dir) must not share a temp path, or the
    # second os.replace() raises FileNotFoundError once the first has
    # already consumed it.
    tmp = file.with_suffix(f".json.{uuid.uuid4().hex}.part")
    tmp.write_text(json.dumps({"cached_at": time.time(), "value": value}))
    os.replace(tmp, file)


def clear() -> int:
    """Delete every cached response, including any incomplete write left
    behind by a crash between a temp file's write and its rename. Returns
    how many files were removed."""
    d = _resolve_dir()
    n = 0
    for pattern in ("*.json", "*.part"):
        for f in d.glob(pattern):
            f.unlink()
            n += 1
    return n
