# scigantic-wwpdb

Explore the wwPDB [Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) (CCD) from Python, by id, no download.

The CCD is the reference chemistry for every small molecule, ligand, ion, and modified residue in the PDB, about 45,000 components. It is the small-molecule reference AlphaFold 3 uses for ligand pose prediction, and what you reach for whenever a structure has a `HETATM` to reason about.

The whole dictionary is a single 45,000-file directory. Mirroring or FUSE-mounting it is awkward: the public autoindex either truncates (files.wwpdb.org caps at 2,000 entries) or is slow to list (the complete mirror is 50k entries). So this library never lists the tree. It fetches exactly what you ask for, over HTTP, in three tiers so exploring stays fast:

```python
import scigantic_wwpdb as ccd

ids = ccd.search("heme")            # 1. ids by name           (RCSB text search)
ccd.find(ids)                       # 2. metadata for the set, ONE request:
#   [Summary('HEC', 'HEME C', formula='C34 H34 Fe N4 O4'), ...]
c = ccd.component("HEC")            # 3. full structure only for what you want
```

`find` batch-fetches name/formula/type/weight/SMILES/InChIKey for a whole result set in a single request (RCSB's data API), so you skim before you fetch. The browsing idiom is `find(search("kinase inhibitor"))`.

## Install

```
pip install scigantic-wwpdb          # core: gemmi + requests
pip install "scigantic-wwpdb[rdkit]" # + RDKit for to_rdkit / to_sdf / depiction
```

## Full structure, RDKit, SDF

```python
c = ccd.component("ATP")            # one small GET, parsed
c.name        # "ADENOSINE-5'-TRIPHOSPHATE"
c.formula     # "C10 H16 N5 O13 P3"
c.atoms[0]    # Atom(id='PG', element='P', charge=0, x=1.2, y=-0.226, z=-6.85, ...)

ccd.components(["ATP", "ADP", "AMP"])   # many, fetched in parallel
mol = ccd.to_rdkit("ATP")               # RDKit Mol with the ideal 3-D coordinates
ccd.to_sdf("ATP", "ATP.sdf")            # write an SDF
```

`to_rdkit` builds the molecule from the dictionary's atoms and bonds, not from a SMILES round-trip, so organometallics like heme (whose coordinate-bond SMILES a plain parser cannot read) still work.

One-liners `smiles / inchi / inchikey / formula / name(id)` use the fast metadata API, not a full CIF fetch.

## The whole dictionary

For bulk work, stream the full `components.cif.gz` once (cached) and parse it with gemmi:

```python
doc = ccd.read_dictionary()         # downloads ~117 MB once, then a gemmi Document
```

For anything less than a few thousand components, `component()` / `components()` are faster and use no local disk.

## Layout

Small, single-purpose modules: `rcsb` (search + batch metadata), `structure` (fetch + parse), `_cif` (gemmi), `chem` (RDKit/SDF), `dictionary` (the bundle), `model` (records). Per-component data is RCSB (`files.rcsb.org/ligands/download/<ID>.cif`); search and batch metadata are the RCSB search + data APIs; the bundle is EBI. All overridable via `SCIGANTIC_CCD_*` environment variables.

## Data and license

Library code: MIT (see `LICENSE`). The CCD data it fetches is wwPDB, released under [CC0 1.0](https://www.wwpdb.org/about/usage-policies); cite the wwPDB when you use it. Built by [Scigantic](https://scigantic.com).
