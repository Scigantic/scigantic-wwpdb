"""The whole CCD in one stream, for bulk sweeps (not 45k GETs)."""
from __future__ import annotations
import os
import gzip

from . import _http

# The full dictionary as one gzip'd mmCIF (~117 MB). EBI mirror serves the
# complete, untruncated file.
DICTIONARY_URL = os.environ.get(
    "SCIGANTIC_CCD_DICTIONARY_URL",
    "https://ftp.ebi.ac.uk/pub/databases/pdb/data/monomers/components.cif.gz")


def load_dictionary(dest: str | None = None, force: bool = False,
                    progress: bool = True) -> str:
    """Stream the full CCD (components.cif.gz, ~111 MB) to a local file once and
    return its path. Cached — re-runs are a no-op unless force=True.

    A single stream is deliberate: benchmarked from a us-east-1 pod, EBI serves
    this file at ~11 MB/s on one connection and *throttles* concurrent byte-range
    requests to ~3-4 MB/s, so splitting the download is ~3x slower, not faster."""
    dest = dest or os.path.join(
        os.environ.get("HOME", "/tmp"), "ccd", "components.cif.gz")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and not force and os.path.getsize(dest) > 0:
        return dest
    with _http.session.get(DICTIONARY_URL, stream=True, timeout=180) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        done = 0
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if progress and total:
                    print(f"\r  components.cif.gz {done >> 20}/{total >> 20} MiB "
                          f"({100 * done // total}%)", end="", flush=True)
        if progress:
            print()
    return dest


def read_dictionary(path: str | None = None):
    """Parse a downloaded components.cif.gz into a gemmi Document (all ~45k
    blocks). Memory-heavy — prefer component(id)/components(ids) for a handful.
    Downloads first if `path` is omitted."""
    import gemmi
    path = path or load_dictionary(progress=False)
    with gzip.open(path, "rt") as fh:
        return gemmi.cif.read_string(fh.read())
