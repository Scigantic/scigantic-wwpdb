"""RDKit / SDF conversion (rdkit is imported lazily, so it's an optional dep)."""
from __future__ import annotations

from .model import Component
from .structure import component

_BOND_ORDER = {"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4, "AROM": 12}


def to_rdkit(ccd_id_or_component, sanitize: bool = True):
    """Build an RDKit Mol from the CCD atoms + bonds, with the ideal 3-D
    coordinates set as a conformer. Built from the dictionary's connectivity, not
    a SMILES round-trip, so organometallics like heme work where a SMILES parser
    returns None. Requires rdkit (`pip install "scigantic-wwpdb[rdkit]"`)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import RWMol
        from rdkit.Geometry import Point3D
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'to_rdkit needs rdkit — pip install "scigantic-wwpdb[rdkit]"') from exc

    comp = ccd_id_or_component if isinstance(ccd_id_or_component, Component) \
        else component(ccd_id_or_component)

    rw = RWMol()
    idx = {}
    for a in comp.atoms:
        atom = Chem.Atom(a.element)
        atom.SetFormalCharge(int(a.charge or 0))
        if a.aromatic:
            atom.SetIsAromatic(True)
        idx[a.id] = rw.AddAtom(atom)
    for b in comp.bonds:
        if b.a1 in idx and b.a2 in idx:
            bt = Chem.BondType.AROMATIC if b.aromatic else \
                Chem.BondType.values[_BOND_ORDER.get(b.order, 1)]
            rw.AddBond(idx[b.a1], idx[b.a2], bt)

    mol = rw.GetMol()
    if any(a.x is not None for a in comp.atoms):
        conf = Chem.Conformer(mol.GetNumAtoms())
        for a in comp.atoms:
            conf.SetAtomPosition(
                idx[a.id], Point3D(a.x or 0.0, a.y or 0.0, a.z or 0.0))
        mol.AddConformer(conf, assignId=True)
    mol.SetProp("_Name", comp.id)
    if sanitize:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass  # some CCD entries (radicals, odd valences) won't fully sanitize
    return mol


def to_sdf(ccd_id, path: str | None = None) -> str:
    """Render a component to an SDF (MDL molblock). Returns the SDF string; also
    writes it to `path` when given."""
    from rdkit import Chem
    sdf = Chem.MolToMolBlock(to_rdkit(ccd_id), kekulize=False)
    if path:
        with open(path, "w") as fh:
            fh.write(sdf)
    return sdf
