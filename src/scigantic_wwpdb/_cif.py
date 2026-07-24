"""Parse a component's mmCIF into a Component (gemmi is the reader)."""
from __future__ import annotations

from .model import Component, Atom, Bond


def cif_block(text: str):
    """Parse mmCIF text to a gemmi block (the raw handle, for anything below)."""
    import gemmi
    return gemmi.cif.read_string(text).sole_block()


def _val(block, tag):
    import gemmi
    raw = block.find_value(tag)
    if raw is None:
        return None
    s = gemmi.cif.as_string(raw)
    return None if s in ("", "?", ".") else s


def _rows(block, prefix, cols):
    import gemmi
    for row in block.find(prefix, cols):
        yield [
            (None if row[i] in ("?", ".") else gemmi.cif.as_string(row[i]))
            for i in range(len(cols))
        ]


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def parse(text: str) -> Component:
    """Parse CCD mmCIF text into a Component (name, formula, atoms with ideal
    3-D coords, bonds, SMILES/InChI). Use when you already have the CIF text;
    component(id) fetches then calls this."""
    block = cif_block(text)
    cid = _val(block, "_chem_comp.id") or ""

    atoms = []
    for r in _rows(block, "_chem_comp_atom.",
                   ["atom_id", "type_symbol", "charge", "pdbx_aromatic_flag",
                    "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal",
                    "pdbx_model_Cartn_z_ideal"]):
        atoms.append(Atom(
            id=r[0], element=(r[1] or "C").capitalize(),
            charge=int(_num(r[2]) or 0),
            aromatic=(r[3] == "Y"),
            x=_num(r[4]), y=_num(r[5]), z=_num(r[6]),
        ))

    bonds = []
    for r in _rows(block, "_chem_comp_bond.",
                   ["atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"]):
        bonds.append(Bond(a1=r[0], a2=r[1], order=(r[2] or "SING"),
                          aromatic=(r[3] == "Y")))

    smi = smi_oe = inch = inchk = None
    for r in _rows(block, "_pdbx_chem_comp_descriptor.",
                   ["type", "program", "descriptor"]):
        t, prog, d = (r[0] or ""), (r[1] or ""), r[2]
        if not d:
            continue
        if t == "SMILES_CANONICAL" and "OpenEye" in prog:
            smi_oe = d
        elif t == "SMILES_CANONICAL":
            smi = d
        elif t == "SMILES" and smi is None:
            smi = d
        elif t == "InChI":
            inch = d
        elif t == "InChIKey":
            inchk = d

    return Component(
        id=cid,
        name=_val(block, "_chem_comp.name"),
        formula=_val(block, "_chem_comp.formula"),
        type=_val(block, "_chem_comp.type"),
        formal_charge=int(_num(_val(block, "_chem_comp.pdbx_formal_charge")) or 0),
        atoms=atoms, bonds=bonds,
        smiles=smi or smi_oe, smiles_openeye=smi_oe,
        inchi=inch, inchikey=inchk, cif_text=text,
    )
