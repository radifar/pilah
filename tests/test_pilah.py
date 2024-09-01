from unittest.mock import Mock, patch

from typer.testing import CliRunner

from pilah.pilah import app


runner = CliRunner()

mock_config_default = Mock()
mock_config_default.load.return_value = "configuration"
mock_config_default.data = dict(
    input = "6hsh.pdb",
    ligand_id = "GOK",
    protein_out = "protein_6hsh.pdb",
    ligand_out = "GOK.pdb",
    ligand_image = "GOK.png",
)

mock_config = Mock()
mock_config.load.return_value = "configuration"
mock_config.data = dict(
    input = "6hsh.pdb",
    ligand_id = "GOK",
    include_metal = "yes",
    protein_out = "protein_6hsh.pdb",
    ligand_out = "GOK.pdb",
    ligand_image = "GOK.png",
    image_size = "large",
    ligand_smiles = "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
    ph = 7.4
)

complex_pdb_block = dict(
    ligand = "ligand_pdb_block",
    protein = "protein_pdb_block"
)

def test_app():
    with (
        patch('pilah.pilah.Config', return_value=mock_config),
        patch('pilah.pilah.extract', return_value=complex_pdb_block),
        patch('pilah.pilah.process_ligand', return_value=["ligand_mol", "ligand_mol_Hs"]),
        patch('pilah.pilah.process_protein', return_value=["protein_mol","ionization_records"]),
        patch('pilah.pilah.renumber_hydrogens'),
        patch('pilah.pilah.rename_hydrogens'),
        patch('pilah.pilah.mol_writer'),
        patch('pilah.pilah.log_writer'),
        patch('pilah.pilah.mol_drawer')
    ):
        result = runner.invoke(app, ["run", "tests/data/config_pdb_gok.txt"])
        assert result.exit_code == 0

def test_app_default_config():
    with (
        patch('pilah.pilah.Config', return_value=mock_config_default),
        patch('pilah.pilah.extract', return_value=complex_pdb_block),
        patch('pilah.pilah.process_ligand', return_value=["ligand_mol", "ligand_mol_Hs"]),
        patch('pilah.pilah.process_protein', return_value=["protein_mol","ionization_records"]),
        patch('pilah.pilah.renumber_hydrogens'),
        patch('pilah.pilah.rename_hydrogens'),
        patch('pilah.pilah.mol_writer'),
        patch('pilah.pilah.log_writer'),
        patch('pilah.pilah.mol_drawer')
    ):
        result = runner.invoke(app, ["run", "tests/data/config_pdb_gok.txt"])
        assert result.exit_code == 0
