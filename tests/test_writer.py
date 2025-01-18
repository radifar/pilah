from tempfile import NamedTemporaryFile

import pytest
from helpers import assert_atom_names_in_residues, assert_image_size, assert_image_svg
from openbabel import openbabel

from pilah.extract import extract
from pilah.protonate import process_ligand, process_protein
from pilah.writer import mol_drawer, mol_writer, rename_hydrogens, renumber_hydrogens


@pytest.fixture
def mol_1e66_protonated():
    """
    This one is special because of its many atom missing and
    incorrect bond length and angle. Therefore the pdb block
    must be generated from advanced filtering.
    """
    data = dict(
        input="tests/data/1e66.pdb",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="HUX",
        include_metal="yes",
        protein_out="protein_6hsh.pdbqt",
        ligand_out="HUX.pdbqt",
        ligand_image="GOK.png",
    )
    extraction_data = extract(data)
    pdb_block_protein_1e66 = extraction_data["protein"]
    mol_1e66, ionization_records = process_protein(data, pdb_block_protein_1e66)

    return mol_1e66, ionization_records


@pytest.fixture
def mol_5nzn_protonated(pdb_block_protein_5nzn):
    empty_dict = {"protein_out": "protein.mol2"}
    mol_5nzn, ionization_records = process_protein(empty_dict, pdb_block_protein_5nzn)

    return mol_5nzn, ionization_records


@pytest.fixture
def mol_6hsh_protonated(pdb_block_protein_6hsh):
    empty_dict = {"protein_out": "protein.mol2"}
    mol_6hsh, ionization_records = process_protein(empty_dict, pdb_block_protein_6hsh)

    return mol_6hsh, ionization_records


@pytest.fixture
def mol_6hsh_protonated_ph11(pdb_block_protein_6hsh):
    """
    To test writing negatively charged CYS in pdbqt format and
    make sure no negatively charged TYR in pdbqt format
    """
    conf_dict = {"ph": 11, "protein_out": "protein.pdbqt"}
    mol_6hsh, ionization_records = process_protein(conf_dict, pdb_block_protein_6hsh)

    return mol_6hsh, ionization_records


@pytest.fixture
def mol_6hsh_protonated_ph3(pdb_block_protein_6hsh):
    """
    To test writing neutral ASP and GLU in pdbqt format
    """
    conf_dict = {"ph": 3, "protein_out": "protein.pdbqt"}
    mol_6hsh, ionization_records = process_protein(conf_dict, pdb_block_protein_6hsh)

    return mol_6hsh, ionization_records


@pytest.fixture
def mol_G39_protonated(pdb_block_ligand_G39):
    empty_dict = {
        "protein_out": "protein.mol2",
        "ligand_smiles": "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O",
    }
    _, mol_G39, _ = process_ligand(empty_dict, pdb_block_ligand_G39)

    return mol_G39


@pytest.fixture
def mol_GOK_protonated(pdb_block_ligand_GOK):
    empty_dict = {
        "protein_out": "protein.mol2",
        "ligand_smiles": "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
    }
    _, mol_GOK, _ = process_ligand(empty_dict, pdb_block_ligand_GOK)

    return mol_GOK


@pytest.mark.parametrize(
    "mol_protonated",
    [
        ("mol_1e66_protonated"),
        ("mol_5nzn_protonated"),
        ("mol_6hsh_protonated"),
    ],
)
def test_renumber_hydrogens(mol_protonated, request):
    mol_protonated, _ = request.getfixturevalue(mol_protonated)
    mol_sorted = renumber_hydrogens(mol_protonated)

    heavy_atom_ids = dict()
    hydrogen_atom_ids = dict()

    for atom in mol_sorted.GetAtoms():
        residue_id = atom.GetPDBResidueInfo().GetResidueNumber()
        atom_id = atom.GetIdx()
        if atom.GetSymbol() == "H":
            if residue_id in hydrogen_atom_ids.keys():
                hydrogen_atom_ids[residue_id].append(atom_id)
            else:
                hydrogen_atom_ids[residue_id] = [atom_id]
        else:
            if residue_id in hydrogen_atom_ids.keys():
                heavy_atom_ids[residue_id].append(atom_id)
            else:
                heavy_atom_ids[residue_id] = [atom_id]

    for key in heavy_atom_ids.keys():
        last_heavy_atom_id = heavy_atom_ids[key][-1]
        for hydrogen_atom_id in hydrogen_atom_ids[key]:
            assert (hydrogen_atom_id - last_heavy_atom_id) < 20


@pytest.mark.parametrize(
    "mol_protonated",
    [
        ("mol_1e66_protonated"),
        ("mol_5nzn_protonated"),
        ("mol_6hsh_protonated"),
    ],
)
def test_rename_hydrogen(mol_protonated, request):
    mol_protonated, _ = request.getfixturevalue(mol_protonated)
    mol_renamed = rename_hydrogens(mol_protonated)

    residue_atom_names_dict = dict()
    for atom in mol_renamed.GetAtoms():
        atom_name = atom.GetMonomerInfo().GetName().strip()

        residue = atom.GetPDBResidueInfo()
        residue_name = residue.GetResidueName()
        residue_number = residue.GetResidueNumber()
        residue_id = (residue_name, residue_number)
        if residue_id in residue_atom_names_dict.keys():
            residue_atom_names_dict[residue_id].add(atom_name)
        else:
            residue_atom_names_dict[residue_id] = {atom_name}

    assert_atom_names_in_residues(residue_atom_names_dict)


@pytest.mark.parametrize("format", [(".pdb"), (".mol2"), (".pdbqt")])
@pytest.mark.parametrize(
    "mol_protonated",
    [
        ("mol_5nzn_protonated"),
        ("mol_6hsh_protonated"),
        ("mol_6hsh_protonated_ph11"),
        ("mol_6hsh_protonated_ph3"),
    ],
)
def test_mol_writer_protein(mol_protonated, format, request):
    with NamedTemporaryFile(suffix=format) as mol_file:
        filename = mol_file.name
        mol_protonated, ionization_records = request.getfixturevalue(mol_protonated)
        mol_renamed = rename_hydrogens(mol_protonated)

        mol_writer(mol_renamed, filename, ionization_records, receptor=True)

        mol_read_by_ob = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.SetInFormat(format[1:])
        converter.ReadFile(mol_read_by_ob, filename)

        assert mol_read_by_ob.NumResidues() > 200


@pytest.mark.parametrize("format", [(".pdb"), (".mol2"), (".sdf"), (".pdbqt")])
@pytest.mark.parametrize(
    "mol_protonated, num_heavy_atoms",
    [
        ("mol_G39_protonated", 20),
        ("mol_GOK_protonated", 29),
    ],
)
def test_mol_writer_ligand(mol_protonated, format, num_heavy_atoms, request):
    with NamedTemporaryFile(suffix=format) as mol_file:
        filename = mol_file.name
        mol_protonated = request.getfixturevalue(mol_protonated)

        mol_writer(mol_protonated, filename)

        mol_read_by_ob = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.SetInFormat(format[1:])
        converter.ReadFile(mol_read_by_ob, filename)

        assert mol_read_by_ob.NumHvyAtoms() == num_heavy_atoms


def test_mol_writer_file_format_not_recognized(mol_GOK_protonated, capfd):
    filename = "ligand.docx"
    mol_writer(mol_GOK_protonated, filename)

    out, err = capfd.readouterr()

    assert out == "Unrecognized file format: docx\n"


def test_mol_drawer_svg(mol_GOK_protonated):
    with NamedTemporaryFile(suffix=".svg") as mol_img:
        filename = mol_img.name
        mol_drawer(mol_GOK_protonated, filename, "small")

        assert_image_svg(filename)


@pytest.mark.parametrize(
    "size, real_size",
    [("small", "small"), ("medium", (420, 420)), ("large", (640, 640))],
)
def test_mol_drawer_png(mol_GOK_protonated, size, real_size):
    with NamedTemporaryFile(suffix=".png") as mol_img:
        filename = mol_img.name
        mol_drawer(mol_GOK_protonated, filename, size)

        assert_image_size(filename, size, real_size)


def test_mol_drawer_size_not_recognized(mol_GOK_protonated, capfd):
    with NamedTemporaryFile(suffix=".png") as mol_img:
        filename = mol_img.name
        mol_drawer(mol_GOK_protonated, filename, "tiny")

        out, err = capfd.readouterr()

    assert out == "Unrecognized image size: tiny, using small size instead\n"


def test_mol_drawer_file_format_not_recognized(mol_GOK_protonated, capfd):
    filename = "mol.gif"
    mol_drawer(mol_GOK_protonated, filename, "small")

    out, err = capfd.readouterr()

    assert out == "Unrecognized image format: gif\n"
