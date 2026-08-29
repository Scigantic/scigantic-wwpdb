"""cache.py, and find()/fetch_cif() wired up to it. No live network: every
request-shaped call here is monkeypatched, matching this repo's existing
convention (see test_explore.py). conftest.py's autouse fixture already
points every test at its own empty cache directory."""
import json
from concurrent.futures import ThreadPoolExecutor

import pytest

import scigantic_wwpdb as ccd
from scigantic_wwpdb import cache, rcsb, structure


def test_cache_enabled_by_default():
    assert cache.is_cache_enabled() is True


_PAYLOAD = {"data": {"chem_comps": [
    {"rcsb_id": "ATP",
     "chem_comp": {"id": "ATP", "name": "ADENOSINE-5'-TRIPHOSPHATE",
                   "formula": "C10 H16 N5 O13 P3", "type": "non-polymer",
                   "formula_weight": 507.18},
     "rcsb_chem_comp_descriptor": {"SMILES_stereo": "c1nc...N",
                                   "InChI": "InChI=1S/...",
                                   "InChIKey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"}},
]}}


def test_find_second_call_does_not_hit_network(monkeypatch):
    calls = []
    monkeypatch.setattr(rcsb._http, "post_json",
                        lambda *a, **k: calls.append(1) or _PAYLOAD)
    first = ccd.find(["ATP"])
    assert len(calls) == 1

    second = ccd.find(["ATP"])
    assert len(calls) == 1  # no new request
    assert second[0].id == first[0].id == "ATP"
    assert second[0].formula == first[0].formula


def test_find_serves_one_liners_from_cache(monkeypatch):
    """name()/formula()/smiles()/inchi()/inchikey() each call find([id]) via
    _one(). Looking up several of those for the same id used to mean several
    identical GraphQL round trips; only the first should touch the network."""
    calls = []
    monkeypatch.setattr(rcsb._http, "post_json",
                        lambda *a, **k: calls.append(1) or _PAYLOAD)
    assert ccd.name("ATP") == "ADENOSINE-5'-TRIPHOSPHATE"
    assert ccd.formula("ATP") == "C10 H16 N5 O13 P3"
    assert ccd.smiles("ATP") == "c1nc...N"
    assert len(calls) == 1


def test_find_unknown_id_is_not_cached(monkeypatch):
    """A row RCSB never returned (bad id, typo) must not be remembered as
    permanently missing: re-querying it should hit the network again."""
    empty = {"data": {"chem_comps": [None]}}
    calls = []
    monkeypatch.setattr(rcsb._http, "post_json",
                        lambda *a, **k: calls.append(1) or empty)
    assert ccd.find(["NOPE"]) == []
    assert ccd.find(["NOPE"]) == []
    assert len(calls) == 2  # neither lookup was cached


def test_fetch_cif_second_call_does_not_hit_network(monkeypatch):
    class _Resp:
        status_code = 200
        text = "data_ATP\n#\n"

        def raise_for_status(self):
            pass

    calls = []
    monkeypatch.setattr(structure._http, "get",
                        lambda *a, **k: calls.append(1) or _Resp())
    first = structure.fetch_cif("ATP")
    assert len(calls) == 1
    second = structure.fetch_cif("ATP")
    assert len(calls) == 1
    assert first == second == "data_ATP\n#\n"


def test_fetch_cif_404_is_not_cached(monkeypatch):
    class _NotFound:
        status_code = 404
        text = ""

        def raise_for_status(self):
            pass

    monkeypatch.setattr(structure._http, "get", lambda *a, **k: _NotFound())
    with pytest.raises(KeyError):
        structure.fetch_cif("ZZZZZ", retries=0)
    assert cache.get("cif:ZZZZZ") is None


def test_entries_expire_after_ttl(tmp_path, monkeypatch):
    cache.enable_cache(cache_dir=str(tmp_path), ttl_days=30)
    calls = []
    monkeypatch.setattr(rcsb._http, "post_json",
                        lambda *a, **k: calls.append(1) or _PAYLOAD)
    ccd.find(["ATP"])
    assert len(calls) == 1

    entry_file = tmp_path / "find_ATP.json"
    entry = json.loads(entry_file.read_text())
    entry["cached_at"] -= 31 * 86400  # backdate past the TTL
    entry_file.write_text(json.dumps(entry))

    ccd.find(["ATP"])
    assert len(calls) == 2  # the expired entry forced a real request


def test_disable_cache_hits_network_every_time(monkeypatch):
    cache.disable_cache()
    try:
        calls = []
        monkeypatch.setattr(rcsb._http, "post_json",
                            lambda *a, **k: calls.append(1) or _PAYLOAD)
        ccd.find(["ATP"])
        ccd.find(["ATP"])
        assert len(calls) == 2
    finally:
        cache.enable_cache()


def test_concurrent_writes_to_the_same_key_never_raise(tmp_path):
    cache.enable_cache(cache_dir=str(tmp_path))
    errors = []

    def write(i):
        try:
            cache.put("cif:SAME", f"payload-{i}")
        except Exception as e:  # noqa: BLE001
            errors.append(e)

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(write, range(16)))

    assert errors == []
    assert list(tmp_path.glob("*.part")) == []  # every temp file was consumed
