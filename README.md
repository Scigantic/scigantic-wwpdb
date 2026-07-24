# scigantic-wwpdb

The wwPDB [Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) (CCD), by id, no download.

The CCD is the reference chemistry for every small molecule, ligand, ion, and modified residue in the PDB, about 45,000 components. It is the small-molecule reference AlphaFold 3 uses for ligand pose prediction, and what you reach for whenever a structure has a `HETATM` to reason about.

The whole dictionary is a single 45,000-file directory. Mirroring or FUSE-mounting it is awkward: the public autoindex either truncates (files.wwpdb.org caps at 2,000 entries, hiding the bundle and most named ligands) or is slow to list (the complete mirror is 50k entries). So this library does the obvious thing instead: it fetches exactly the component you ask for, directly over HTTP, by id.

## Install

```
pip install scigantic-wwpdb          # core: gemmi + requests
pip install "scigantic-wwpdb[rdkit]" # + RDKit for to_rdkit / to_sdf / depiction
```

## Quickstart

```python
import scigantic_wwpdb as ccd

atp = ccd.component("ATP")          # one small GET, parsed
atp.name        # "ADENOSINE-5'-TRIPHOSPHATE"
atp.formula     # "C10 H16 N5 O13 P3"
atp.smiles      # canonical SMILES
atp.inchikey    # "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"
atp.atoms[0]    # Atom(id='PG', element='P', charge=0, x=1.2, y=-0.226, z=-6.85, ...)

ccd.smiles("HEM")                   # one-liners: smiles / inchi / inchikey / formula / name
mol = ccd.to_rdkit("ATP")           # RDKit Mol with the ideal 3-D coordinates
ccd.to_sdf("ATP", "ATP.sdf")        # write an SDF
```

`to_rdkit` builds the molecule from the dictionary's atoms and bonds, not from a SMILES round-trip, so organometallics like heme (whose coordinate-bond SMILES a plain parser cannot read) still work.

Per-component data is the modern mmCIF served by RCSB (`files.rcsb.org/ligands/download/<ID>.cif`). Overridable with the `SCIGANTIC_CCD_LIGAND_BASE` environment variable.

## The whole dictionary

For bulk work, stream the full `components.cif.gz` once (cached) and parse it with gemmi:

```python
path = ccd.load_dictionary()        # ~117 MB, one stream, cached under ~/ccd/
doc  = ccd.read_dictionary(path)    # a gemmi Document of all ~45k blocks
```

For anything less than a few thousand components, per-id `component()` calls are faster and use no local disk.

## Search

`search("adenosine")` needs an id-to-name catalog index (`SCIGANTIC_CCD_CATALOG`, a JSON map). Without one, look ids up directly with `component(id)`.

## Data and license

Library code: MIT (see `LICENSE`). The CCD data it fetches is wwPDB, released under [CC0 1.0](https://www.wwpdb.org/about/usage-policies); cite the wwPDB when you use it. Built by [Scigantic](https://scigantic.com).
