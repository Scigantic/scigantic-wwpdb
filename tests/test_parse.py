"""Offline parsing + RDKit conversion, using real CIF fixtures (no network)."""
import pathlib

import pytest

import scigantic_wwpdb as ccd

FIX = pathlib.Path(__file__).parent / "fixtures"


def _load(fname):
    return (FIX / fname).read_text()


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


def test_component_summary_view():
    s = ccd.parse(_load("ATP.cif")).summary()
    assert isinstance(s, ccd.Summary)
    assert s.id == "ATP" and s.formula == "C10 H16 N5 O13 P3"


def test_repr_is_informative():
    c = ccd.parse(_load("ATP.cif"))
    assert "ATP" in repr(c) and "atoms=47" in repr(c)


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


def test_component_url_uppercases():
    assert ccd.component_url("atp") == \
        "https://files.rcsb.org/ligands/download/ATP.cif"
