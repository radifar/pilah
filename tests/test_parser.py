import pytest

from pilah.parser import Config


def test_parser_all(capfd):
    config_all = Config()
    config_all.load("tests/data/config_pdb_gok_all.txt")
    out, err = capfd.readouterr()

    assert config_all.data["input"] == "tests/data/6hsh.pdb"
    assert out == ""


def test_parser_warning(capfd):
    config_warning = Config()
    config_warning.load("tests/data/config_pdb_gok_warning.txt")
    out, err = capfd.readouterr()

    assert out == "Warning: 'bogus_option' option is not recognized\n"

def test_parser_missing_mandatory():
    with pytest.raises(SystemExit) as excinfo:
        config_missing = Config()
        config_missing.load("tests/data/config_pdb_gok_missing.txt")

    assert excinfo.value.code == "Missing mandatory option: ligand_chain"