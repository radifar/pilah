from unittest.mock import Mock

import pytest


@pytest.fixture
def mock_config_default():
    config_default = Mock()
    config_default.load.return_value = "configuration"
    config_default.data = dict(
        input = "tests/data/6hsh.pdb",
        ligand_id = "GOK",
        protein_out = "protein_6hsh.pdb",
        ligand_out = "GOK.pdb",
        ligand_image = "GOK.png",
    )

    return config_default

@pytest.fixture
def mock_config():
    config = Mock()
    config.load.return_value = "configuration"
    config.data = dict(
        input = "tests/data/6hsh.pdb",
        ligand_id = "GOK",
        include_metal = "yes",
        protein_out = "protein_6hsh.pdb",
        ligand_out = "GOK.pdb",
        ligand_image = "GOK.png",
        image_size = "large",
        ligand_smiles = "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
        ph = 7.4
    )

    return config

@pytest.fixture
def mock_config_cif():
    config_cif = Mock()
    config_cif.load.return_value = "configuration"
    config_cif.data = dict(
        input = "tests/data/6hsh.cif",
        ligand_id = "GOK",
        protein_out = "protein_6hsh.pdb",
        ligand_out = "GOK.pdb",
        ligand_image = "GOK.png",
    )

    return config_cif

@pytest.fixture
def mock_config_unknown_format():
    config_cif = Mock()
    config_cif.load.return_value = "configuration"
    config_cif.data = dict(
        input = "tests/data/6hsh.mol2",
        ligand_id = "GOK",
        protein_out = "protein_6hsh.pdb",
        ligand_out = "GOK.pdb",
        ligand_image = "GOK.png",
    )

    return config_cif