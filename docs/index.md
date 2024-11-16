# PiLAH

PiLAH stands for **Protein Ligand Abstraction and Hydrogenization tool**, it is a cheminformatic tool for protein and ligand preparation.
Given a PDB file it will separate the protein and the ligand, fix the bond order of the ligand, protonate protein and ligand separately, then write them to separate file.

It uses Biopython, RDKit, Meeko, and Openbabel for the general cheminformatics processing, then it will use Dimorphite-DL and pKAI to protonate small molecule and protein respectively.

## Features

### Ligand Processing
- Choose ligand using its chain ID and residue ID.
- Bond order correction using SMILES template.
- Bond order correction for ligand with missing atoms.
- Protonation using Dimorphite-DL.
- Choose the pH environment.

### Protein Processing
- Choose protein using its chain ID (multiple chain is possible).
- Able to include/exclude metal.
- Removing residue with missing atoms.
- Removing residue with incorrect bond length or angle.
- Renumbering residues for structure with Insertion Code.
- Choose residue coordinate with Alternate Location (altloc id).
- Protonation using pKAI model.
- Choose the pH environment.

### Input Output
- Input in PDB format.
- Ligand output in PDB, MOL2, SDF, and PDBQT format
- Protein output in PDB, MOL2, and PDBQT format.
- 2D image (PNG or SVG format) of generated ligand for easier inspection.

## License

`PiLAH` was created by Muhammad Radifar and Enade Perdana Istyastono. It is licensed under the terms of the Apache License 2.0 license.

## Credits

`PiLAH` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).


```{toctree}
:hidden:

self
```

```{toctree}
:hidden:
:caption: Usage

usage/installation
usage/gettingstarted
usage/usecases
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Information

changelog.md
contributing.md
conduct.md
```