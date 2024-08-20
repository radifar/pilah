import io
import sys
import warnings

from Bio import BiopythonWarning
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import FastMMCIFParser


# Metal list retrieved from
# Lin, Geng-Yu, et al. "MESPEUS: a database of metal coordination groups
# in proteins." Nucleic Acids Research 52.D1 (2024): D483-D493.
# https://academic.oup.com/nar/article/52/D1/D483/7370098

metal_list = ["CA", "CO", "CU", "FE", "K", "MG", "MN", "NA", "NI",
              "ZN", "AG", "AL", "AU", "BA", "BE", "CD", "CR", "CS",
              "GA", "GD", "HG", "IR", "LI", "MO", "PB", "PD", "PR",
              "PT", "RB", "RE", "RH", "RU", "SR", "TB", "TL", "U",
              "V", "W", "Y", "YB"]

class ResidueSelect(Select):
    def __init__(self, chain="A", include_metal="no"):
        self.chain_select = chain
        self.include_metal = include_metal

    def accept_chain(self, chain):
        return chain.id == self.chain_select

    def accept_residue(self, residue):
        if self.include_metal == "no":
            return residue.id[0] == " "
        elif self.include_metal == "yes":
            return residue.id[0] == " " or residue.get_resname() in metal_list

class LigandSelect(Select):
    def __init__(self, res_name, chain="A"):
        self.res_name = res_name
        self.chain_select = chain

    def accept_chain(self, chain):
        return chain.id == self.chain_select

    def accept_residue(self, residue):
        return residue.get_resname() == self.res_name

def extract(data):
    ligand_id = data["ligand_id"]
    chain = data.get("chain", "A")
    include_metal = data.get("include_metal", "no")

    input_file = data["input"]
    extension = input_file.split(".")[-1]
    if extension == "pdb":
        parser = PDBParser()
    elif extension == "cif":
        parser = FastMMCIFParser()
    else:
        sys.exit(f"Input file format {extension} is not recognized")

    # Suppress warning from Biopython
    # https://www.biostars.org/p/251583/
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        structure = parser.get_structure('structure', input_file)[0]

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    
    protein_handle = io.StringIO()
    ligand_handle = io.StringIO()

    pdbio.save(protein_handle, ResidueSelect(chain=chain, include_metal=include_metal))
    pdbio.save(ligand_handle, LigandSelect(ligand_id, chain=chain))

    protein_handle.seek(0)
    ligand_handle.seek(0)

    complex_block = {
        "protein": protein_handle.read(),
        "ligand": ligand_handle.read()
        }

    return complex_block
