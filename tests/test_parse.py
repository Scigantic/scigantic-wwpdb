import os
import pathlib

import pytest

import scigantic_wwpdb as ccd

FIX = pathlib.Path(__file__).parent / "fixtures"


def _load(fname):
    return (FIX / fname).read_text()


# ── offline parsing (real CIF fixtures, no network) ─────────────────────────
def test_parse_atp():
    c = ccd.parse(_load("ATP.cif"))
    assert c.id == "ATP"
    assert c.name == "ADENOSINE-5'-TRIPHOSPHATE"
    assert c.formula == "C10 H16 N5 O13 P3"
    assert len(c.atoms) == 47
    assert len(c.bonds) == 49
    assert c.smiles and "n" in c.smiles.lower()
    assert c.inchikey == "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"
    assert c.inchi and c.inchi.startswith("InChI=")
    first = c.atoms[0]
    assert first.element == "P"
    assert first.x is not None and first.y is not None and first.z is not None


def test_parse_hem_has_iron():
    c = ccd.parse(_load("HEM.cif"))
    assert c.id == "HEM"
    assert "Fe" in (c.formula or "")
    assert any(a.element == "Fe" for a in c.atoms)


def test_repr_is_informative():
    c = ccd.parse(_load("ATP.cif"))
    assert "ATP" in repr(c) and "atoms=47" in repr(c)


# ── rdkit conversion (built from atoms+bonds, so metals work) ───────────────
def test_to_rdkit_atp_carries_3d():
    pytest.importorskip("rdkit")
    m = ccd.to_rdkit(ccd.parse(_load("ATP.cif")))
    assert m.GetNumAtoms() == 47
    assert m.GetNumConformers() == 1  # ideal coords


def test_to_rdkit_hem_keeps_iron():
    pytest.importorskip("rdkit")
    m = ccd.to_rdkit(ccd.parse(_load("HEM.cif")))
    assert any(a.GetSymbol() == "Fe" for a in m.GetAtoms())


def test_to_sdf_is_molblock():
    pytest.importorskip("rdkit")
    sdf = ccd.to_sdf(ccd.parse(_load("ATP.cif")))
    assert "V2000" in sdf or "V3000" in sdf


# ── plumbing ────────────────────────────────────────────────────────────────
def test_component_url_uppercases():
    assert ccd.component_url("atp") == \
        "https://files.rcsb.org/ligands/download/ATP.cif"


def test_search_without_catalog_explains():
    with pytest.raises(RuntimeError, match="catalog"):
        ccd.search("adenosine")


# ── live smoke tests (network) — opt in with SCIGANTIC_LIVE_TESTS=1 ─────────
@pytest.mark.skipif(not os.environ.get("SCIGANTIC_LIVE_TESTS"),
                    reason="set SCIGANTIC_LIVE_TESTS=1 to run network tests")
def test_live_component_roundtrip():
    c = ccd.component("NAD")
    assert c.name and c.formula and c.smiles


@pytest.mark.skipif(not os.environ.get("SCIGANTIC_LIVE_TESTS"),
                    reason="set SCIGANTIC_LIVE_TESTS=1 to run network tests")
def test_live_missing_id_raises():
    with pytest.raises(KeyError):
        ccd.component("ZZZZZ")
