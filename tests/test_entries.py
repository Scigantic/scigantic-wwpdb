"""Entry search: query construction, pagination, and the additive filter.

Every request-shaped call here is monkeypatched, matching this repo's existing
tests — nothing in this file touches the network.
"""
import pytest

from scigantic_wwpdb import entries


def _capture(monkeypatch, response):
    """Swap post_json for a recorder; returns the list it appends payloads to."""
    seen = []

    def fake(url, payload, *a, **k):
        seen.append((url, payload))
        return response(len(seen) - 1) if callable(response) else response

    monkeypatch.setattr(entries._http, "post_json", fake)
    return seen


def test_requires_a_filter_beyond_method(monkeypatch):
    seen = _capture(monkeypatch, {})
    # bare call, and method-only, both select ~200k entries and are refused
    with pytest.raises(ValueError):
        entries.search_entries()
    with pytest.raises(ValueError):
        entries.search_entries(method="X-ray")
    with pytest.raises(ValueError):
        entries.count_entries(method="electron microscopy")
    assert seen == [], "must not call RCSB when refusing the query"


def test_builds_expected_filter_nodes(monkeypatch):
    seen = _capture(monkeypatch, {"result_set": [], "total_count": 0})
    entries.search_entries(released_after="2025-06-30", resolution_max=2.5,
                           with_ligand=True, limit=5)
    _, payload = seen[0]
    nodes = payload["query"]["nodes"]
    attrs = [n["parameters"]["attribute"] for n in nodes]
    assert "rcsb_accession_info.initial_release_date" in attrs
    assert "rcsb_entry_info.resolution_combined" in attrs
    assert "rcsb_entry_info.experimental_method" in attrs
    # with_ligand must constrain BOTH a small molecule and a protein chain,
    # otherwise protein-free entries carrying an ion come back as "complexes"
    assert "rcsb_entry_info.nonpolymer_entity_count" in attrs
    assert "rcsb_entry_info.polymer_entity_count_protein" in attrs
    assert payload["return_type"] == "entry"
    assert payload["request_options"]["results_content_type"] == ["experimental"]


def test_single_filter_is_not_wrapped_in_a_group(monkeypatch):
    seen = _capture(monkeypatch, {"result_set": []})
    entries.search_entries(method=None, resolution_max=2.0, limit=1)
    _, payload = seen[0]
    assert payload["query"]["type"] == "terminal"


def test_method_can_be_dropped(monkeypatch):
    seen = _capture(monkeypatch, {"result_set": []})
    entries.search_entries(released_after="2020-01-01", method=None, limit=1)
    attrs = [n["parameters"]["attribute"] for n in seen[0][1]["query"]["nodes"]] \
        if seen[0][1]["query"].get("nodes") else \
        [seen[0][1]["query"]["parameters"]["attribute"]]
    assert "rcsb_entry_info.experimental_method" not in attrs


def test_paginates_past_one_page(monkeypatch):
    # three pages of 1000, then a short page that ends it
    def resp(i):
        n = 1000 if i < 2 else 400
        base = i * 1000
        return {"result_set": [{"identifier": f"E{base + k:05d}"} for k in range(n)]}

    seen = _capture(monkeypatch, resp)
    ids = entries.search_entries(released_after="2025-06-30", limit=5000)
    assert len(ids) == 2400, "short page must end pagination, not truncate at 1000"
    assert len(seen) == 3
    assert seen[0][1]["request_options"]["paginate"]["start"] == 0
    assert seen[1][1]["request_options"]["paginate"]["start"] == 1000
    assert len(set(ids)) == len(ids), "pages must not overlap"


def test_limit_is_respected_exactly(monkeypatch):
    _capture(monkeypatch, lambda i: {
        "result_set": [{"identifier": f"X{k}"} for k in range(1000)]})
    assert len(entries.search_entries(released_after="2025-01-01", limit=7)) == 7


def test_no_hits_returns_empty_not_error(monkeypatch):
    # RCSB answers an unsatisfiable query with HTTP 204 and no body, which
    # _http.post_json maps to {}. That must read as "none", not blow up.
    _capture(monkeypatch, {})
    assert entries.search_entries(released_after="2099-01-01", limit=10) == []


def test_count_entries_reads_total(monkeypatch):
    seen = _capture(monkeypatch, {"total_count": 8765, "result_set": []})
    assert entries.count_entries(released_after="2025-06-30") == 8765
    assert seen[0][1]["request_options"]["paginate"]["rows"] == 1, "count pulls no ids"


def test_count_entries_zero_on_empty_body(monkeypatch):
    _capture(monkeypatch, {})
    assert entries.count_entries(released_after="2099-01-01") == 0


def test_entry_ligands_strips_additives(monkeypatch):
    _capture(monkeypatch, {"data": {"entries": [
        {"rcsb_id": "1ABC", "nonpolymer_entities": [
            {"nonpolymer_comp": {"chem_comp": {"id": "ATP"}}},
            {"nonpolymer_comp": {"chem_comp": {"id": "MG"}}},
            {"nonpolymer_comp": {"chem_comp": {"id": "HOH"}}},
        ]},
        {"rcsb_id": "2DEF", "nonpolymer_entities": [
            {"nonpolymer_comp": {"chem_comp": {"id": "SO4"}}},
        ]},
    ]}})
    got = entries.entry_ligands(["1ABC", "2DEF"])
    assert got["1ABC"] == ["ATP"]
    assert got["2DEF"] == [], "an entry whose only ligands are additives is empty"


def test_entry_ligands_keeps_additives_when_asked(monkeypatch):
    _capture(monkeypatch, {"data": {"entries": [
        {"rcsb_id": "1ABC", "nonpolymer_entities": [
            {"nonpolymer_comp": {"chem_comp": {"id": "MG"}}},
        ]},
    ]}})
    assert entries.entry_ligands(["1abc"], additives=True)["1ABC"] == ["MG"]


def test_entry_ligands_unknown_entry_maps_to_empty(monkeypatch):
    _capture(monkeypatch, {"data": {"entries": [None]}})
    assert entries.entry_ligands(["9ZZZ"]) == {"9ZZZ": []}


def test_entry_ligands_no_ids_makes_no_call(monkeypatch):
    seen = _capture(monkeypatch, {})
    assert entries.entry_ligands([]) == {}
    assert seen == []


def test_heavy_atom_count_from_formula():
    # real CCD formulas, which are proper-cased: ATP, magnesium, sulfate, water
    assert entries._heavy_atoms("C10 H16 N5 O13 P3") == 31
    assert entries._heavy_atoms("Mg") == 1
    assert entries._heavy_atoms("Fe") == 1
    assert entries._heavy_atoms("O4 S") == 5
    assert entries._heavy_atoms("H2 O") == 1
    assert entries._heavy_atoms(None) == 0
    # deuterium is still hydrogen for this purpose
    assert entries._heavy_atoms("D2 O") == 1


def test_heavy_atoms_needs_proper_cased_formula():
    # documents the one way to misuse this: passing the upper-case CCD *id* as if it
    # were a formula. 'MG' is M + G, which is why the docstring warns against it.
    assert entries._heavy_atoms("MG") == 2
    assert entries._heavy_atoms("Mg") == 1


def test_is_interesting_ligand_rejects_additives(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("additives must be rejected without a lookup")

    monkeypatch.setattr(entries, "find", boom)
    assert entries.is_interesting_ligand("HOH") is False
    assert entries.is_interesting_ligand("SO4") is False


def test_is_interesting_ligand_uses_size(monkeypatch):
    class S:
        def __init__(self, f):
            self.formula = f

    monkeypatch.setattr(entries, "find", lambda ids: [S("C21 H28 N7 O14 P2")])
    assert entries.is_interesting_ligand("ATP") is True
    monkeypatch.setattr(entries, "find", lambda ids: [S("C2 H6 O")])
    assert entries.is_interesting_ligand("EOH") is False
    monkeypatch.setattr(entries, "find", lambda ids: [])
    assert entries.is_interesting_ligand("NOPE") is False
