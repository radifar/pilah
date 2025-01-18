import io
import warnings
from copy import deepcopy
from tempfile import NamedTemporaryFile

import pytest
from Bio import BiopythonWarning
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from rdkit import Chem

from pilah.extract import (
    LigandSelect,
    LigandSelectByResNum,
    PreResidueSelect,
    ResidueSelect,
    extract,
    fix_insertion,
    get_healthy_residues,
)

parser = PDBParser(PERMISSIVE=True)


def test_fix_insertion():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        UPA_structure = parser.get_structure("1SQT", "tests/data/1sqt.pdb")[0]

    chain = ["A"]
    ligand_id = "UI3"

    total_insertion, renumber_residue_map = fix_insertion(UPA_structure, chain, ligand_id)

    assert total_insertion == {"A": 1}
    assert renumber_residue_map[("A", "GLY", 23, "A")] == 24
    assert UPA_structure["A"][(" ", 25, " ")].get_resname() == "GLY"
    assert UPA_structure["A"][(" ", 26, " ")].get_resname() == "SER"
    assert UPA_structure["A"][(" ", 27, " ")].get_resname() == "VAL"


def test_fix_insertion_multichain():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        UPA_structure = parser.get_structure("1YPE", "tests/data/1ype.pdb")[0]

    chain = ["L", "H", "I"]
    ligand_id = "UIP"

    total_insertion, renumber_residue_map = fix_insertion(UPA_structure, chain, ligand_id)
    assert total_insertion == {"L": 13, "H": 23, "I": 0}
    assert renumber_residue_map[("L", "ALA", 1, "B")] == 1
    assert renumber_residue_map[("L", "ASP", 1, "A")] == 2
    assert renumber_residue_map[("L", "CYS", 1, " ")] == 3
    assert renumber_residue_map[("L", "ILE", 14, "K")] == 27
    assert renumber_residue_map[("H", "PHE", 245, " ")] == 268
    assert renumber_residue_map[("I", "LEU", 10, " ")] == 10


def test_preresidue_select():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        athrombine_structure = parser.get_structure("1D3D", "tests/data/1d3d.pdb")[0]

    pdbio = PDBIO()
    pdbio.set_structure(athrombine_structure)

    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A"], "no")
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        athrombine_filtered = parser.get_structure("filtered", pre_protein_handle)[0]

    pre_protein_handle.truncate()
    pre_residue_filter = PreResidueSelect(["A"], "yes")
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        athrombine_filtered_metal = parser.get_structure("filtered", pre_protein_handle)[0]

    residue_list = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]

    metal_list = [
        "CA",
        "CO",
        "CU",
        "FE",
        "K",
        "MG",
        "MN",
        "NA",
        "NI",
        "ZN",
        "AG",
        "AL",
        "AU",
        "BA",
        "BE",
        "CD",
        "CR",
        "CS",
        "GA",
        "GD",
        "HG",
        "IR",
        "LI",
        "MO",
        "PB",
        "PD",
        "PR",
        "PT",
        "RB",
        "RE",
        "RH",
        "RU",
        "SR",
        "TB",
        "TL",
        "U",
        "V",
        "W",
        "Y",
        "YB",
    ]

    for residue in athrombine_filtered["A"]:
        assert residue.get_resname() in residue_list

    for residue in athrombine_filtered_metal["A"]:
        resname = residue.get_resname()
        assert (resname in residue_list) or (resname in metal_list)


@pytest.fixture
def ache_structure():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        structure = parser.get_structure("1E66", "tests/data/1e66.pdb")[0]

    return structure


@pytest.fixture
def hdac_structure():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        HDAC_structure = parser.get_structure("6HSH", "tests/data/6hsh.pdb")

    return HDAC_structure


@pytest.fixture
def akt1_structure():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        HDAC_structure = parser.get_structure("3CQW", "tests/data/3cqw.pdb")

    return HDAC_structure


@pytest.fixture
def hiv_protease_structure():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        HDAC_structure = parser.get_structure("1IZH", "tests/data/1izh.pdb")

    return HDAC_structure


@pytest.fixture
def hdac_io(hdac_structure):
    HDAC_io = PDBIO()
    HDAC_io.set_structure(hdac_structure)

    return HDAC_io


@pytest.fixture
def hdac_healthy_residues(hdac_io):
    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A", "B", "C", "D"])
    hdac_io.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()

    healthy_residue_dict = get_healthy_residues(pre_pdb_block)
    print(healthy_residue_dict["HIS"])

    return healthy_residue_dict


def test_get_healthy_residues(ache_structure):
    pdbio = PDBIO()
    pdbio.set_structure(ache_structure)

    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A"], "no")
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()

    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    residue_with_missing_atoms = {
        "ALA": [536],
        "ARG": [46, 47],
        "ASN": [42, 253, 257, 382],
        "ASP": [365],
        "GLN": [162, 488, 526],
        "GLU": [89, 260, 268, 299, 350, 434, 455, 489],
        "LEU": [7],
        "LYS": [52, 192, 270, 413, 454, 478, 511],
        "SER": [490],
    }

    residue_with_incorrect_bond_length_angle = {
        "ARG": [88, 468],
        "GLU": [247],
        "LYS": [325],
        "SER": [228],
    }

    missing_atoms_res_keys = residue_with_missing_atoms.keys()
    incorrect_bond_length_angle_keys = residue_with_incorrect_bond_length_angle.keys()
    for res_name, res_nums in healthy_residue_dict.items():
        if res_name in missing_atoms_res_keys:
            for res_num in res_nums:
                assert res_num not in residue_with_missing_atoms[res_name]
        if res_name in incorrect_bond_length_angle_keys:
            for res_num in res_nums:
                assert res_num not in residue_with_incorrect_bond_length_angle[res_name]


def test_residue_select(hdac_io, hdac_healthy_residues):
    # Based on PDB REMARK data of 6HSH
    residue_num = 447
    missing_residue_chain_A = 43
    missing_residue_chain_B = 37
    missing_residue_chain_C = 32
    missing_residue_chain_D = 42

    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        hdac_io.save(
            file_name,
            ResidueSelect(chain=["A"], healthy_residue_dict=hdac_healthy_residues),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            protein_chain_A = parser.get_structure("6hsh", file_name)
        residue_num_A = len(list(protein_chain_A.get_residues()))

        hdac_io.save(
            file_name,
            ResidueSelect(chain=["B"], healthy_residue_dict=hdac_healthy_residues),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            protein_chain_B = parser.get_structure("6hsh", file_name)
        residue_num_B = len(list(protein_chain_B.get_residues()))

        hdac_io.save(
            file_name,
            ResidueSelect(chain=["C"], healthy_residue_dict=hdac_healthy_residues),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            protein_chain_C = parser.get_structure("6hsh", file_name)
        residue_num_C = len(list(protein_chain_C.get_residues()))

        hdac_io.save(
            file_name,
            ResidueSelect(chain=["D"], healthy_residue_dict=hdac_healthy_residues),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            protein_chain_D = parser.get_structure("6hsh", file_name)
        residue_num_D = len(list(protein_chain_D.get_residues()))
        temp_pdb_file.close()

    assert residue_num_A == (residue_num - missing_residue_chain_A)
    assert residue_num_B == (residue_num - missing_residue_chain_B)
    assert residue_num_C == (residue_num - missing_residue_chain_C)
    assert residue_num_D == (residue_num - missing_residue_chain_D)


def test_residue_select_broken_residues(ache_structure):
    pdbio = PDBIO()
    pdbio.set_structure(ache_structure)

    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A"], "no")
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()

    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    residue_filter = ResidueSelect(["A"], "no", healthy_residue_dict)

    protein_handle = io.StringIO()
    pdbio.save(protein_handle, residue_filter)
    protein_handle.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        ache_ready = parser.get_structure("ache_ready", protein_handle)[0]

    missing_residue_chain_A = 10
    residue_with_missing_atom = 28
    residue_with_incorrect_bond_length_angle = 5
    residue_num_A = 543
    expected_residue_num_A = (
        residue_num_A - missing_residue_chain_A - residue_with_missing_atom - residue_with_incorrect_bond_length_angle
    )

    filtered_residue_num = len(list(ache_ready["A"].get_residues()))

    assert expected_residue_num_A == filtered_residue_num


def test_residue_select_with_metal(hdac_io, hdac_healthy_residues):
    residue_num = 447
    missing_residue_chain_A = 43

    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        residue_filter = ResidueSelect(chain=["A"], include_metal="yes", healthy_residue_dict=hdac_healthy_residues)
        hdac_io.save(file_name, residue_filter)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            protein_chain_A = parser.get_structure("6hsh", file_name)
        residue_list_A = list(protein_chain_A.get_residues())
        residue_num_A = len(residue_list_A)
        temp_pdb_file.close()

    assert residue_num_A == (residue_num - missing_residue_chain_A + 3)
    assert residue_list_A[-3].resname == "ZN"
    assert residue_list_A[-2].resname == "K"
    assert residue_list_A[-1].resname == "K"


def test_residue_select_broken_residues_with_metal(hdac_structure):
    pdbio = PDBIO()
    pdbio.set_structure(hdac_structure)

    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A"], "yes")
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()

    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    residue_filter = ResidueSelect(["A"], "yes", healthy_residue_dict)

    protein_handle = io.StringIO()
    pdbio.save(protein_handle, residue_filter)
    protein_handle.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        hdac_ready = parser.get_structure("hdac_ready", protein_handle)[0]

    missing_residue_chain_A = 43
    metal_residue = 3
    residue_num_A = 447
    expected_residue_num_A = residue_num_A - missing_residue_chain_A + metal_residue

    filtered_residue_num = len(list(hdac_ready["A"].get_residues()))

    assert expected_residue_num_A == filtered_residue_num


def test_residue_select_with_altloc(akt1_structure):
    altloc_list_a = "A:A182LYS A:A205SER".split()
    altloc_list_b = "B:A182LYS B:A205SER".split()

    akt1_structure_b = deepcopy(akt1_structure)
    pdbio_a = PDBIO()
    pdbio_a.set_structure(akt1_structure)
    pdbio_b = PDBIO()
    pdbio_b.set_structure(akt1_structure_b)

    pre_protein_handle = io.StringIO()
    pre_residue_filter = PreResidueSelect(["A"], "no")
    pdbio_a.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()

    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    residue_filter_a = ResidueSelect(["A"], "no", healthy_residue_dict, altloc_list_a)
    residue_filter_b = ResidueSelect(["A"], "no", healthy_residue_dict, altloc_list_b)

    protein_handle_a = io.StringIO()
    protein_handle_b = io.StringIO()
    pdbio_a.save(protein_handle_a, residue_filter_a)
    pdbio_b.save(protein_handle_b, residue_filter_b)
    protein_handle_a.seek(0)
    protein_handle_b.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        akt1_ready_a = parser.get_structure("akt1_ready_a", protein_handle_a)[0]["A"]
        akt1_ready_b = parser.get_structure("akt1_ready_b", protein_handle_b)[0]["A"]

    resid_ASN148_a = akt1_ready_a[148]
    resid_ASN148_b = akt1_ready_b[148]
    resid_LYS182_a = akt1_ready_a[182]
    resid_LYS182_b = akt1_ready_b[182]
    resid_205SER_a = akt1_ready_a[205]
    resid_205SER_b = akt1_ready_b[205]

    asn_148_expected = [7.154, 1.079, -6.201]
    asn_148_actual_a = resid_ASN148_a["OD1"].get_coord()
    asn_148_actual_b = resid_ASN148_b["OD1"].get_coord()

    lys_182_expected_a = [16.908, -4.690, 4.425]
    lys_182_expected_b = [19.416, -4.104, 8.118]
    lys_182_actual_a = resid_LYS182_a["NZ"].get_coord()
    lys_182_actual_b = resid_LYS182_b["NZ"].get_coord()

    ser_205_expected_a = [-5.409, -7.680, 20.298]
    ser_205_expected_b = [-4.735, -6.415, 18.583]
    ser_205_actual_a = resid_205SER_a["OG"].get_coord()
    ser_205_actual_b = resid_205SER_b["OG"].get_coord()

    assert all([axis in asn_148_actual_a for axis in asn_148_expected])
    assert all([axis in asn_148_actual_b for axis in asn_148_expected])
    assert all([axis in lys_182_actual_a for axis in lys_182_expected_a])
    assert all([axis in lys_182_actual_b for axis in lys_182_expected_b])
    assert all([axis in ser_205_actual_a for axis in ser_205_expected_a])
    assert all([axis in ser_205_actual_b for axis in ser_205_expected_b])


def test_ligand_select(hdac_io):
    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        hdac_io.save(file_name, LigandSelect(res_name="GOK", chain="A"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            quisinostat_chain_A = parser.get_structure("6hsh", file_name)
            quisinostat_atom_num = len(list(quisinostat_chain_A.get_atoms()))
        temp_pdb_file.close()

    assert quisinostat_atom_num == 29


def test_ligand_select_with_altloc(hiv_protease_structure):
    altloc_list_b = ["B:B1001Q50"]

    hiv_protease_structure_b = deepcopy(hiv_protease_structure)
    pdbio_a = PDBIO()
    pdbio_a.set_structure(hiv_protease_structure)
    pdbio_b = PDBIO()
    pdbio_b.set_structure(hiv_protease_structure_b)

    ligand_filter_a = LigandSelect("Q50", "B")
    ligand_filter_b = LigandSelect("Q50", "B", altloc_list_b)

    q50_handle_a = io.StringIO()
    q50_handle_b = io.StringIO()
    pdbio_a.save(q50_handle_a, ligand_filter_a)
    pdbio_b.save(q50_handle_b, ligand_filter_b)
    q50_handle_a.seek(0)
    q50_handle_b.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        q50_ready_a = parser.get_structure("q50_ready_a", q50_handle_a)[0]["B"]
        q50_ready_b = parser.get_structure("q50_ready_b", q50_handle_b)[0]["B"]

    resid_Q50_a = q50_ready_a["H_Q50", 1001, " "]
    resid_Q50_b = q50_ready_b["H_Q50", 1001, " "]

    q50_expected_a = [11.236, 24.803, 5.276]
    q50_expected_b = [15.089, 21.732, 5.390]
    q50_actual_a = resid_Q50_a["N8"].get_coord()
    q50_actual_b = resid_Q50_b["N8"].get_coord()

    assert all([axis in q50_actual_a for axis in q50_expected_a])
    assert all([axis in q50_actual_b for axis in q50_expected_b])


def test_ligand_select_by_resnum():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        double_4HI_structure = parser.get_structure("4HI", "tests/data/4hi_multi.pdb")

    ligand_pdbio = PDBIO()
    ligand_pdbio.set_structure(double_4HI_structure)
    ligand_handle = io.StringIO()

    ligand_pdbio.save(ligand_handle, LigandSelectByResNum("2"))
    ligand_handle.seek(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        selected_4HI_structure = parser.get_structure("4HI", ligand_handle)[0]["B"]

    resid_4HI_2 = selected_4HI_structure["H_4HI", 2, " "]

    assert len(list(selected_4HI_structure.get_residues())) == 1
    assert len(list(resid_4HI_2.get_atoms())) == 36


def test_extract(mock_config_default):
    extraction_data = extract(mock_config_default.data)

    protein_block = extraction_data["protein"]
    ligand_block = extraction_data["ligand"]
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)

    assert protein_mol.GetNumAtoms() > 3200
    assert ligand_mol.GetNumAtoms() == 29


def test_extract_cif(mock_config_cif):
    extraction_data = extract(mock_config_cif.data)

    protein_block = extraction_data["protein"]
    ligand_block = extraction_data["ligand"]
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)

    assert protein_mol.GetNumAtoms() > 3200
    assert ligand_mol.GetNumAtoms() == 29


def test_extract_multichain_protein(mock_config_multichain_protein):
    extraction_data = extract(mock_config_multichain_protein.data)

    protein_block = extraction_data["protein"]
    protein_mol = Chem.MolFromPDBBlock(protein_block)

    protein_chains = Chem.SplitMolByPDBChainId(protein_mol)

    assert all([chain in protein_chains.keys() for chain in ["A", "B"]])
    assert protein_mol.GetNumAtoms() > 6400


def test_extract_multiligand(mock_config_multiligand, capfd):
    extraction_data = extract(mock_config_multiligand.data)

    ligand_block = extraction_data["ligand"]
    multiple_ligand = extraction_data["multiple_ligand"]

    ligand_mol = Chem.MolFromPDBBlock(ligand_block)
    ligand_resnum = ligand_mol.GetAtomWithIdx(1).GetPDBResidueInfo().GetResidueNumber()

    out, err = capfd.readouterr()

    assert ligand_resnum == 1
    assert ligand_mol.GetNumAtoms() == 36
    assert multiple_ligand[0] == "B"
    assert multiple_ligand[1] == "4HI"
    assert multiple_ligand[2] == "1, 2"
    assert "Multiple ligand detected" in out
    assert "chain B" in out
    assert "1, 2" in out


def test_extract_multiligand_select_resnum(mock_config_multiligand_select_resnum):
    extraction_data = extract(mock_config_multiligand_select_resnum.data)

    ligand_block = extraction_data["ligand"]
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)
    ligand_resnum = ligand_mol.GetAtomWithIdx(1).GetPDBResidueInfo().GetResidueNumber()

    assert ligand_resnum == 2
    assert ligand_mol.GetNumAtoms() == 36
    assert "multiple_ligand" not in extraction_data.keys()


def test_extract_with_insertion_code(mock_config_with_insertion):
    extraction_data = extract(mock_config_with_insertion.data)

    total_insertion = extraction_data["total_insertion"]
    renumber_residue_map = extraction_data["renumber_residue_map"]

    assert total_insertion["A"] == 1
    assert renumber_residue_map[("A", "GLY", 23, "A")] == 24
    assert renumber_residue_map[("A", "SER", 25, " ")] == 26


def test_extract_wrong_ligand_chain_or_id(mock_config_wrong_ligand_chain):
    err = "\nError: Ligand empty make sure that ligand_chain and ligand_id are correct"
    with pytest.raises(SystemExit, match=err):
        extract(mock_config_wrong_ligand_chain.data)


def test_extract_unknown_format(mock_config_unknown_format):
    err = "Input file format mol2 is not recognized"
    with pytest.raises(SystemExit, match=err):
        extract(mock_config_unknown_format.data)
