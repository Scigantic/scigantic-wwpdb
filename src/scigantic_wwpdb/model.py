"""Typed records the library returns."""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class Summary:
    """Lightweight metadata for a component (from RCSB's batch data API) — the
    fast skim layer: name/formula/type/weight/SMILES/InChIKey, no atoms or
    coordinates, and no per-component fetch."""
    id: str
    name: str | None = None
    formula: str | None = None
    type: str | None = None
    weight: float | None = None
    smiles: str | None = None
    inchi: str | None = None
    inchikey: str | None = None

    def __repr__(self):
        return f"Summary({self.id!r}, {self.name!r}, formula={self.formula!r})"


@dataclass
class Atom:
    id: str
    element: str
    charge: int = 0
    x: float | None = None  # ideal (computed) coordinates, Angstroms
    y: float | None = None
    z: float | None = None
    aromatic: bool = False


@dataclass
class Bond:
    a1: str
    a2: str
    order: str = "SING"  # SING | DOUB | TRIP | AROM
    aromatic: bool = False


@dataclass
class Component:
    """Full component parsed from its mmCIF: metadata + every atom (with ideal
    3-D coords), every bond, and the SMILES/InChI descriptors."""
    id: str
    name: str | None = None
    formula: str | None = None
    type: str | None = None
    formal_charge: int = 0
    atoms: list = field(default_factory=list)
    bonds: list = field(default_factory=list)
    smiles: str | None = None
    smiles_openeye: str | None = None
    inchi: str | None = None
    inchikey: str | None = None
    cif_text: str = ""

    def __repr__(self):
        return (f"Component({self.id!r}, name={self.name!r}, "
                f"formula={self.formula!r}, atoms={len(self.atoms)}, "
                f"bonds={len(self.bonds)}, smiles={bool(self.smiles)})")

    def summary(self) -> Summary:
        return Summary(id=self.id, name=self.name, formula=self.formula,
                       type=self.type, smiles=self.smiles, inchi=self.inchi,
                       inchikey=self.inchikey)
