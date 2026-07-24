"""scigantic_wwpdb — the wwPDB Chemical Component Dictionary (CCD), by id, no download.

The CCD defines every small molecule, ligand, ion, and modified residue that
appears in the PDB — the reference chemistry AlphaFold 3 leans on for ligand
pose prediction, and what you reach for whenever a structure has a HETATM you
need to reason about. It is ~45,000 components.

Rather than mirror or mount that (a single 45k-file directory whose public index
either truncates or is glacial to list), this library fetches exactly the
component you ask for, directly over HTTP, by id:

  component("ATP")       -> parsed CCD entry: name, formula, atoms+coords, bonds
  smiles("ATP")          -> canonical SMILES        inchi("ATP") / inchikey("ATP")
  to_rdkit("ATP")        -> an RDKit Mol with 3-D ideal coords (needs rdkit)
  to_sdf("ATP")          -> an SDF string / file
  parse(text)            -> parse mmCIF text you already have
  load_dictionary()      -> stream the whole CCD (components.cif.gz) once
  search("adenosine")    -> component ids whose name matches (uses a catalog)

One small GET per component (~a few KB), or one stream for the full dictionary.
Nothing here is a locked widget — every function is a few lines you can read,
copy, and bend.
"""
from __future__ import annotations
import os
import gzip
import functools
from dataclasses import dataclass, field

import requests

__version__ = "0.1.0"

__all__ = [
    "component", "component_url", "fetch_cif", "parse", "cif_block",
    "smiles", "inchi", "inchikey", "formula", "name",
    "to_rdkit", "to_sdf",
    "load_dictionary", "read_dictionary", "search", "catalog",
    "Component", "Atom", "Bond",
]

# ── locations (all overridable so a demo can point at a staging mirror) ──────
# Per-component MODERN mmCIF (atoms with ideal coords + SMILES/InChI descriptors)
# — RCSB serves one file per component, addressable by id, no listing needed.
LIGAND_BASE = os.environ.get(
    "SCIGANTIC_CCD_LIGAND_BASE", "https://files.rcsb.org/ligands/download"
)
# The whole dictionary as one gzip'd mmCIF (~117 MB). EBI mirror serves the
# complete, untruncated file.
DICTIONARY_URL = os.environ.get(
    "SCIGANTIC_CCD_DICTIONARY_URL",
    "https://ftp.ebi.ac.uk/pub/databases/pdb/data/monomers/components.cif.gz",
)
# Optional id -> {name, formula, type} index for search(). Absent is fine —
# search() just explains how to build/point at one.
CATALOG_URL = os.environ.get("SCIGANTIC_CCD_CATALOG")

_UA = {"User-Agent": f"scigantic-wwpdb/{__version__} (+https://scigantic.com; mailto:support@scigantic.com)"}
_session = requests.Session()
_session.headers.update(_UA)


def component_url(ccd_id: str) -> str:
    """RCSB URL for a single component's modern mmCIF."""
    return f"{LIGAND_BASE}/{str(ccd_id).strip().upper()}.cif"


def fetch_cif(ccd_id: str, retries: int = 2) -> str:
    """Fetch one component's raw mmCIF text (a few KB)."""
    url = component_url(ccd_id)
    last = None
    for _ in range(retries + 1):
        try:
            r = _session.get(url, timeout=30)
            if r.status_code == 404:
                raise KeyError(f"CCD component {ccd_id!r} not found ({url})")
            r.raise_for_status()
            return r.text
        except KeyError:
            raise
        except Exception as exc:  # transient network / 5xx — retry
            last = exc
    raise last


# ── parsing (gemmi is the CIF reader) ───────────────────────────────────────
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


def component(ccd_id: str) -> Component:
    """Fetch and parse one CCD component by id."""
    comp = parse(fetch_cif(ccd_id))
    if not comp.id:
        comp.id = str(ccd_id).upper()
    return comp


# ── one-liner accessors ─────────────────────────────────────────────────────
def smiles(ccd_id: str) -> str | None:
    """Canonical SMILES for a component."""
    return component(ccd_id).smiles


def inchi(ccd_id: str) -> str | None:
    return component(ccd_id).inchi


def inchikey(ccd_id: str) -> str | None:
    return component(ccd_id).inchikey


def formula(ccd_id: str) -> str | None:
    return component(ccd_id).formula


def name(ccd_id: str) -> str | None:
    return component(ccd_id).name


# ── RDKit / SDF (optional — rdkit is imported lazily) ───────────────────────
_BOND_ORDER = {"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4, "AROM": 12}


def to_rdkit(ccd_id_or_component, sanitize: bool = True):
    """Build an RDKit Mol from the CCD atoms + bonds, with the ideal 3-D
    coordinates set as a conformer. Faithful to the dictionary (not a SMILES
    round-trip, so organometallics like heme work). Requires rdkit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import RWMol
        from rdkit.Geometry import Point3D
    except ImportError as exc:  # pragma: no cover
        raise ImportError("to_rdkit needs rdkit — `pip install rdkit`") from exc

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
    mol = to_rdkit(ccd_id)
    sdf = Chem.MolToMolBlock(mol, kekulize=False)
    if path:
        with open(path, "w") as fh:
            fh.write(sdf)
    return sdf


# ── the whole dictionary (one stream, not 45k GETs) ─────────────────────────
def load_dictionary(dest: str | None = None, force: bool = False,
                    progress: bool = True) -> str:
    """Stream the full CCD (components.cif.gz, ~117 MB) to a local file once and
    return its path. Cached — re-runs are a no-op unless force=True."""
    dest = dest or os.path.join(
        os.environ.get("HOME", "/tmp"), "ccd", "components.cif.gz")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and not force and os.path.getsize(dest) > 0:
        return dest
    with _session.get(DICTIONARY_URL, stream=True, timeout=120) as r:
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
    blocks). Memory-heavy — prefer component(id) for a handful. Downloads first
    if `path` is omitted."""
    import gemmi
    path = path or load_dictionary(progress=False)
    with gzip.open(path, "rt") as fh:
        return gemmi.cif.read_string(fh.read())


# ── search (needs the optional catalog index) ───────────────────────────────
@functools.lru_cache(maxsize=1)
def catalog() -> dict:
    """Load the id -> {name, formula, type} catalog index (for search). Empty
    dict when no catalog is configured/reachable."""
    if not CATALOG_URL:
        return {}
    try:
        return _session.get(CATALOG_URL, timeout=30).json()
    except Exception:
        return {}


def search(query: str, limit: int = 25) -> list:
    """Component ids whose id/name matches `query` (substring, case-insensitive).
    Uses the catalog index; if none is configured, raises with how to build one."""
    cat = catalog()
    if not cat:
        raise RuntimeError(
            "search() needs the CCD catalog index (id -> name). None is "
            "configured. Set SCIGANTIC_CCD_CATALOG to a catalog.json, or look up "
            "a known id directly, e.g. component('ATP').")
    q = query.strip().lower()
    hits = [cid for cid, meta in cat.items()
            if q in cid.lower() or q in (meta.get("name", "") or "").lower()]
    return sorted(hits)[:limit]
