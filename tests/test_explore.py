"""The explore layer: find() (batch metadata) and search() (RCSB text search).

Parsing of the API responses is tested offline with canned payloads; the real
network round-trips are opt-in via SCIGANTIC_LIVE_TESTS=1.
"""
import os

import pytest

import scigantic_wwpdb as ccd
from scigantic_wwpdb import rcsb, structure


# ── find(): parse a canned GraphQL payload (no network) ─────────────────────
def test_find_parses_graphql(monkeypatch):
    payload = {"data": {"chem_comps": [
        {"rcsb_id": "ATP",
         "chem_comp": {"id": "ATP", "name": "ADENOSINE-5'-TRIPHOSPHATE",
                       "formula": "C10 H16 N5 O13 P3", "type": "non-polymer",
                       "formula_weight": 507.18},
         "rcsb_chem_comp_descriptor": {"SMILES_stereo": "c1nc...N",
                                       "InChI": "InChI=1S/...",
                                       "InChIKey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"}},
        None,  # RCSB returns null for an unknown id — must be dropped
    ]}}
    monkeypatch.setattr(rcsb._http, "post_json", lambda *a, **k: payload)
    got = ccd.find(["ATP", "NOPE"])
    assert len(got) == 1
    s = got[0]
    assert isinstance(s, ccd.Summary)
    assert s.id == "ATP" and s.formula == "C10 H16 N5 O13 P3"
    assert s.weight == 507.18 and s.inchikey.endswith("KQYNXXCUSA-N")


def test_find_empty_input_no_call(monkeypatch):
    called = []
    monkeypatch.setattr(rcsb._http, "post_json",
                        lambda *a, **k: called.append(1) or {})
    assert ccd.find([]) == []
    assert not called  # empty id list must not hit the network


def test_find_accepts_single_id(monkeypatch):
    monkeypatch.setattr(rcsb._http, "post_json", lambda *a, **k: {"data": {
        "chem_comps": [{"rcsb_id": "GOL", "chem_comp": {"id": "GOL"},
                        "rcsb_chem_comp_descriptor": {}}]}})
    got = ccd.find("gol")
    assert len(got) == 1 and got[0].id == "GOL"


# ── search(): parse a canned Search API payload (no network) ────────────────
def test_search_parses_result_set(monkeypatch):
    monkeypatch.setattr(rcsb._http, "post_json", lambda *a, **k: {
        "total_count": 2,
        "result_set": [{"identifier": "GOL"}, {"identifier": "DGA"}]})
    assert ccd.search("glycerol") == ["GOL", "DGA"]


def test_search_no_hits_returns_empty(monkeypatch):
    # A search with no results comes back 204 -> {} from post_json.
    monkeypatch.setattr(rcsb._http, "post_json", lambda *a, **k: {})
    assert ccd.search("zzzznope") == []


# ── components(): parallel fetch preserves order (network stubbed) ──────────
def test_components_parallel_preserves_order(monkeypatch):
    monkeypatch.setattr(structure, "component", lambda cid: f"C:{cid}")
    assert ccd.components(["ATP", "HEM", "NAD"]) == ["C:ATP", "C:HEM", "C:NAD"]


def _flaky(cid):
    if cid == "BAD":
        raise KeyError("nope")
    return f"C:{cid}"


def test_components_omits_failures_by_default(monkeypatch):
    # One bad id must not sink the whole batch — it becomes None in place.
    monkeypatch.setattr(structure, "component", _flaky)
    assert ccd.components(["ATP", "BAD", "NAD"]) == ["C:ATP", None, "C:NAD"]


def test_components_errors_raise_is_fail_fast(monkeypatch):
    monkeypatch.setattr(structure, "component", _flaky)
    with pytest.raises(KeyError):
        ccd.components(["ATP", "BAD"], errors="raise")


# ── live smoke tests (network) — opt in with SCIGANTIC_LIVE_TESTS=1 ─────────
LIVE = pytest.mark.skipif(not os.environ.get("SCIGANTIC_LIVE_TESTS"),
                          reason="set SCIGANTIC_LIVE_TESTS=1 to run network tests")


@LIVE
def test_live_find_batch():
    got = {s.id: s for s in ccd.find(["ATP", "HEM", "NAD", "GOL"])}
    assert set(got) == {"ATP", "HEM", "NAD", "GOL"}
    assert got["HEM"].formula and "Fe" in got["HEM"].formula
    assert got["ATP"].smiles and got["ATP"].inchikey


@LIVE
def test_live_search_then_find():
    ids = ccd.search("glycerol", limit=5)
    assert "GOL" in ids
    assert any(s.id == "GOL" for s in ccd.find(ids))


@LIVE
def test_live_components_parallel():
    comps = ccd.components(["ATP", "ADP", "AMP"])
    assert [c.id for c in comps] == ["ATP", "ADP", "AMP"]
    assert all(c.atoms and c.smiles for c in comps)


@LIVE
def test_live_missing_id_raises():
    with pytest.raises(KeyError):
        ccd.component("ZZZZZ")
