from unittest.mock import Mock

import pytest
from rdkit import Chem

pytest.register_assert_rewrite("helpers")


@pytest.fixture
def mock_config_default():
    config_default = Mock()
    config_default.load.return_value = "configuration"
    config_default.data = dict(
        input="tests/data/6hsh.pdb",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="GOK",
        protein_out="protein_6hsh.pdb",
        ligand_out="GOK.pdb",
        ligand_image="GOK.png",
    )

    return config_default


@pytest.fixture
def mock_config():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/6hsh.pdb",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="GOK",
        include_metal="yes",
        protein_out="protein_6hsh.pdb",
        ligand_out="GOK.pdb",
        ligand_image="GOK.png",
        image_size="large",
        ligand_smiles="Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_multichain_protein():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/6hsh.pdb",
        protein_chain="A B",
        ligand_chain="A",
        ligand_id="GOK",
        include_metal="yes",
        protein_out="protein_6hsh.pdb",
        ligand_out="GOK.pdb",
        ligand_image="GOK.png",
        image_size="large",
        ligand_smiles="Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_multiligand():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/3ccw.pdb",
        protein_chain="A",
        ligand_chain="B",
        ligand_id="4HI",
        include_metal="yes",
        protein_out="protein_3ccw.pdb",
        ligand_out="4HI.pdb",
        ligand_image="4HI.png",
        image_size="large",
        ligand_smiles="CC(C)c1c(nc(n1CC[C@H](C[C@H](CC(=O)O)O)O)c2ccc(cc2)F)C(=O)NCc3ccccc3",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_multiligand_select_resnum():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/3ccw.pdb",
        protein_chain="A",
        ligand_chain="B",
        ligand_id="4HI",
        ligand_res_num="2",
        include_metal="yes",
        protein_out="protein_3ccw.pdb",
        ligand_out="4HI.pdb",
        ligand_image="4HI.png",
        image_size="large",
        ligand_smiles="CC(C)c1c(nc(n1CC[C@H](C[C@H](CC(=O)O)O)O)c2ccc(cc2)F)C(=O)NCc3ccccc3",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_with_insertion():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/1sqt.pdb",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="UI3",
        include_metal="yes",
        protein_out="protein_1sqt.pdb",
        ligand_out="UI3.pdb",
        ligand_image="UI3.png",
        image_size="large",
        ligand_smiles="[H]/N=C(/c1ccc2ccc(c(c2c1)c3cnn(c3)S(=O)(=O)C)OC)\\N",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_wrong_ligand_chain():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input="tests/data/1sqt.pdb",
        protein_chain="A",
        ligand_chain="B",
        ligand_id="UI3",
        include_metal="yes",
        protein_out="protein_1sqt.pdb",
        ligand_out="UI3.pdb",
        ligand_image="UI3.png",
        image_size="large",
        ligand_smiles="[H]/N=C(/c1ccc2ccc(c(c2c1)c3cnn(c3)S(=O)(=O)C)OC)\\N",
        ph=7.4,
    )

    return config


@pytest.fixture
def mock_config_cif():
    config_cif = Mock()
    config_cif.load.return_value = "configuration"
    config_cif.data = dict(
        input="tests/data/6hsh.cif",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="GOK",
        protein_out="protein_6hsh.pdb",
        ligand_out="GOK.pdb",
        ligand_image="GOK.png",
    )

    return config_cif


@pytest.fixture
def mock_config_unknown_format():
    config_cif = Mock()
    config_cif.load.return_value = "configuration"
    config_cif.data = dict(
        input="tests/data/6hsh.mol2",
        protein_chain="A",
        ligand_chain="A",
        ligand_id="GOK",
        protein_out="protein_6hsh.pdb",
        ligand_out="GOK.pdb",
        ligand_image="GOK.png",
    )

    return config_cif


@pytest.fixture
def pdb_block_1j3f():
    with open("tests/data/1j3f.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_1xmk():
    with open("tests/data/1xmk.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_2xji():
    with open("tests/data/2xji.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_3ucy():
    with open("tests/data/3ucy.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_protein_1e66():
    with open("tests/data/protein_1e66.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_protein_5nzn():
    with open("tests/data/protein_5nzn.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_protein_6hsh():
    with open("tests/data/protein_6hsh.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_ligand_HUX():
    with open("tests/data/ligand_HUX_1e66.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_ligand_G39():
    with open("tests/data/ligand_G39_5nzn.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_ligand_GOK():
    with open("tests/data/ligand_GOK_6hsh.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_ligand_B49():
    with open("tests/data/ligand_B49_3g0e.pdb") as f:
        return f.read()


@pytest.fixture
def pdb_block_ligand_UI3():
    with open("tests/data/ligand_UI3_1sqt.pdb") as f:
        return f.read()


@pytest.fixture
def cr_apo_mbs(pdb_block_1j3f):
    mol = Chem.MolFromPDBBlock(pdb_block_1j3f)
    mol_withH = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_noH = Chem.RemoveHs(mol_withH, updateExplicitCount=True)
    return mol_noH


@pytest.fixture
def zb_ADAR1(pdb_block_1xmk):
    mol = Chem.MolFromPDBBlock(pdb_block_1xmk)
    mol_withH = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_noH = Chem.RemoveHs(mol_withH, updateExplicitCount=True)
    return mol_noH


@pytest.fixture
def methanobactin(pdb_block_2xji):
    mol = Chem.MolFromPDBBlock(pdb_block_2xji)
    mol_withH = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_noH = Chem.RemoveHs(mol_withH, updateExplicitCount=True)
    return mol_noH


@pytest.fixture
def calmodulin(pdb_block_3ucy):
    mol = Chem.MolFromPDBBlock(pdb_block_3ucy)
    mol_withH = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_noH = Chem.RemoveHs(mol_withH, updateExplicitCount=True)
    return mol_noH


@pytest.fixture
def ionizable_peptide():
    return Chem.MolFromSequence("KRHDEYC")
