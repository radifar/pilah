from tempfile import NamedTemporaryFile
import warnings

from Bio import BiopythonWarning
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
import pytest
from rdkit import Chem

from pilah.extract import ResidueSelect, LigandSelect, extract


parser = PDBParser(PERMISSIVE=True)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    HDAC_structure = parser.get_structure("6HSH", "tests/data/6hsh.pdb")
HDAC_io = PDBIO()
HDAC_io.set_structure(HDAC_structure)

# Based on PDB REMARK data of 6HSH
residue_num = 447
missing_residue_chain_A = 43
missing_residue_chain_B = 37
missing_residue_chain_C = 32
missing_residue_chain_D = 42


def test_residue_select():
    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        HDAC_io.save(file_name, ResidueSelect())
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            protein_chain_A = parser.get_structure("6hsh", file_name)
        residue_num_A = len(list(protein_chain_A.get_residues()))

        HDAC_io.save(file_name, ResidueSelect(chain="B"))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            protein_chain_B = parser.get_structure("6hsh", file_name)
        residue_num_B = len(list(protein_chain_B.get_residues()))

        HDAC_io.save(file_name, ResidueSelect(chain="C"))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            protein_chain_C = parser.get_structure("6hsh", file_name)
        residue_num_C = len(list(protein_chain_C.get_residues()))

        HDAC_io.save(file_name, ResidueSelect(chain="D"))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            protein_chain_D = parser.get_structure("6hsh", file_name)
        residue_num_D = len(list(protein_chain_D.get_residues()))
        temp_pdb_file.close()
    
    assert residue_num_A == (residue_num - missing_residue_chain_A)
    assert residue_num_B == (residue_num - missing_residue_chain_B)
    assert residue_num_C == (residue_num - missing_residue_chain_C)
    assert residue_num_D == (residue_num - missing_residue_chain_D)

def test_residue_select_with_metal():
    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        HDAC_io.save(file_name, ResidueSelect(include_metal="yes"))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            protein_chain_A = parser.get_structure("6hsh", file_name)
        residue_list_A = list(protein_chain_A.get_residues())
        residue_num_A = len(residue_list_A)
        temp_pdb_file.close()
    
    assert residue_num_A == (residue_num - missing_residue_chain_A + 3)
    assert residue_list_A[-3].resname == "ZN"
    assert residue_list_A[-2].resname == "K"
    assert residue_list_A[-1].resname == "K"

def test_ligand_select():
    with NamedTemporaryFile() as temp_pdb_file:
        file_name = temp_pdb_file.name
        HDAC_io.save(file_name, LigandSelect("GOK"))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            quisinostat_chain_A = parser.get_structure("6hsh", file_name)
            quisinostat_atom_num = len(list(quisinostat_chain_A.get_atoms()))
        temp_pdb_file.close()
    
    assert quisinostat_atom_num == 29

def test_extract(mock_config_default):
    complex_block = extract(mock_config_default.data)
    
    protein_block = complex_block["protein"]
    ligand_block = complex_block["ligand"]
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)

    assert protein_mol.GetNumAtoms() == 3235
    assert ligand_mol.GetNumAtoms() == 29

def test_extract_cif(mock_config_cif):
    complex_block = extract(mock_config_cif.data)
    
    protein_block = complex_block["protein"]
    ligand_block = complex_block["ligand"]
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)

    assert protein_mol.GetNumAtoms() == 3235
    assert ligand_mol.GetNumAtoms() == 29

def test_extract_unknown_format(mock_config_unknown_format):
    with pytest.raises(SystemExit) as excinfo:
        extract(mock_config_unknown_format.data)

    assert excinfo.value.code == 'Input file format mol2 is not recognized'
