from unittest.mock import patch

from typer.testing import CliRunner

from pilah.pilah import app


runner = CliRunner()

complex_pdb_block = dict(ligand="ligand_pdb_block", protein="protein_pdb_block")


def test_app(mock_config):
    with (
        patch("pilah.pilah.Config", return_value=mock_config),
        patch("pilah.pilah.extract", return_value=complex_pdb_block),
        patch(
            "pilah.pilah.process_ligand",
            return_value=["ligand_mol", "ligand_mol_Hs", "ligand_processing_log"],
        ),
        patch(
            "pilah.pilah.process_protein",
            return_value=["protein_mol", "ionization_records"],
        ),
        patch("pilah.pilah.renumber_hydrogens"),
        patch("pilah.pilah.rename_hydrogens"),
        patch("pilah.pilah.mol_writer"),
        patch("pilah.pilah.log_writer"),
        patch("pilah.pilah.mol_drawer"),
    ):
        result = runner.invoke(app, ["run", "tests/data/config_pdb_gok.txt"])
        assert result.exit_code == 0


def test_app_default_config(mock_config_default):
    with (
        patch("pilah.pilah.Config", return_value=mock_config_default),
        patch("pilah.pilah.extract", return_value=complex_pdb_block),
        patch(
            "pilah.pilah.process_ligand",
            return_value=["ligand_mol", "ligand_mol_Hs", "ligand_processing_log"],
        ),
        patch(
            "pilah.pilah.process_protein",
            return_value=["protein_mol", "ionization_records"],
        ),
        patch("pilah.pilah.renumber_hydrogens"),
        patch("pilah.pilah.rename_hydrogens"),
        patch("pilah.pilah.mol_writer"),
        patch("pilah.pilah.log_writer"),
        patch("pilah.pilah.mol_drawer"),
    ):
        result = runner.invoke(app, ["run", "tests/data/config_pdb_gok.txt"])
        assert result.exit_code == 0
