import io
import sys
import warnings

from Bio import BiopythonWarning
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import FastMMCIFParser
from rdkit import Chem
from rich.console import Console

from pilah.atom_names import atom_names_dict


# Metal list retrieved from
# Lin, Geng-Yu, et al. "MESPEUS: a database of metal coordination groups
# in proteins." Nucleic Acids Research 52.D1 (2024): D483-D493.
# https://academic.oup.com/nar/article/52/D1/D483/7370098

metal_list = ["CA", "CO", "CU", "FE", "K", "MG", "MN", "NA", "NI",
              "ZN", "AG", "AL", "AU", "BA", "BE", "CD", "CR", "CS",
              "GA", "GD", "HG", "IR", "LI", "MO", "PB", "PD", "PR",
              "PT", "RB", "RE", "RH", "RU", "SR", "TB", "TL", "U",
              "V", "W", "Y", "YB"]

# side chain smarts patterns retrieved from:
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
residue_smarts_dict = {
    "ARG": Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]"),
    "ASN": Chem.MolFromSmarts("[CH2X4][CX3](=[OX1])[NX3H2]"), # also hits GLN 
    "ASP": Chem.MolFromSmarts("[CH2X4][CX3](=[OX1])[OH0-,OH]"), # also hits GLU
    "CYS": Chem.MolFromSmarts("[CH2X4][$([SX2H,SX1H0-]),$([#16X2H0][#16X2H0])]"),
    "GLN": Chem.MolFromSmarts("[CH2X4][CH2X4][CX3](=[OX1])[NX3H2]"),
    "GLU": Chem.MolFromSmarts("[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]"),
    "HIS": Chem.MolFromSmarts("[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1"),
    "ILE": Chem.MolFromSmarts("[CHX4]([CH3X4])[CH2X4][CH3X4]"),
    "LEU": Chem.MolFromSmarts("[CH2X4][CHX4]([CH3X4])[CH3X4]"),
    "LYS": Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]"),
    "MET": Chem.MolFromSmarts("[CH2X4][CH2X4][SX2][CH3X4]"),
    "PHE": Chem.MolFromSmarts("[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1"),
    "PRO": Chem.MolFromSmarts("N1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[O,N]"),
    "SER": Chem.MolFromSmarts("[CH2X4][OX2H]"),
    "THR": Chem.MolFromSmarts("[CHX4]([CH3X4])[OX2H]"),
    "TRP": Chem.MolFromSmarts("[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12"),
    "TYR": Chem.MolFromSmarts("[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1"),
}

class PreResidueSelect(Select):
    def __init__(self, chain="A", include_metal="no"):
        self.chain_select = chain
        self.include_metal = include_metal

    def accept_chain(self, chain):
        return chain.id == self.chain_select
    
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return True
        elif self.include_metal == "yes":
            return residue.get_resname() in metal_list
        else:
            return False

class ResidueSelect(Select):
    def __init__(self, chain="A", include_metal="no", out_format="", healthy_residue_dict=dict()):
        self.chain_select = chain
        self.include_metal = include_metal
        self.out_format = out_format
        self.healthy_residue_dict = healthy_residue_dict
        self.res_w_missing_atoms = []
        self.res_w_incorrect_bond_length_angle = []

    def accept_chain(self, chain):
        return chain.id == self.chain_select

    def accept_residue(self, residue):
        # remove residues with missing atom if the output format is pdbqt
        if (self.out_format == "pdbqt") and residue.id[0] == " ":
            residue_atom_names = set()
            for atom in residue:
                atom_name = atom.get_name()
                residue_atom_names.add(atom_name)
            residue_name = residue.get_resname()
            residue_number = residue.get_id()[1]
            ref_atom_names = atom_names_dict[residue_name]
            if residue_atom_names == ref_atom_names or residue_atom_names == ref_atom_names.union({"OXT"}):
                if residue_name in ["GLY", "VAL", "ALA"]:
                    return True
                elif residue_number in self.healthy_residue_dict[residue_name]:
                    return True
                else:
                    residue_id = (residue_name, residue_number)
                    self.res_w_incorrect_bond_length_angle.append(residue_id)
                    return False
            else:
                residue_id = (residue_name, residue_number)
                self.res_w_missing_atoms.append(residue_id)
                return False
        elif residue.id[0] == " ":
            return True
        elif self.include_metal == "yes":
            return residue.get_resname() in metal_list
        else:
            return False

class LigandSelect(Select):
    def __init__(self, res_name, chain="A"):
        self.res_name = res_name
        self.chain_select = chain

    def accept_chain(self, chain):
        return chain.id == self.chain_select

    def accept_residue(self, residue):
        return residue.get_resname() == self.res_name

def display_removed_residues(residue_filter):
    console = Console()

    if residue_filter.res_w_missing_atoms:
        console.print("[bold deep_pink2]\n     Some residues with missing atoms were removed:[/bold deep_pink2]")
        console.print("[red3]     residue_name residue_number[/red3]")
        for residue_name, residue_number in residue_filter.res_w_missing_atoms:
            console.print(f"[deep_pink1]         {residue_name}          {residue_number}[/deep_pink1]")
    if residue_filter.res_w_incorrect_bond_length_angle:
        console.print("[bold deep_pink2]\n     Some residues with incorrect bond length and angle were removed:[/bold deep_pink2]")
        console.print("[red3]     residue_name residue_number[/red3]")
        for residue_name, residue_number in residue_filter.res_w_incorrect_bond_length_angle:
            console.print(f"[deep_pink1]         {residue_name}          {residue_number}[/deep_pink1]")

def get_healthy_residues(pre_pdb_block):
    """
    This function only works on filtered PDB (using Biopython).
    Using this function directly on original PDB could give inaccurate result
    as some residues may covalently linked with other residue (phosphate, heme, NAG)
    """
    healthy_residue_dict = {
        "ARG": [],
        "ASN": [],
        "ASP": [],
        "CYS": [],
        "GLN": [],
        "GLU": [],
        "HIS": [],
        "ILE": [],
        "LEU": [],
        "LYS": [],
        "MET": [],
        "PHE": [],
        "PRO": [],
        "SER": [],
        "THR": [],
        "TRP": [],
        "TYR": [],
    }

    mol = Chem.MolFromPDBBlock(pre_pdb_block)
    mol = Chem.AddHs(mol)
    for residue_key, smarts in residue_smarts_dict.items():
        match_indices = mol.GetSubstructMatches(smarts)
        for indices in match_indices:
            first_atom = mol.GetAtomWithIdx(indices[0])
            residue_info = first_atom.GetPDBResidueInfo()
            residue_name = residue_info.GetResidueName()
            residue_number = residue_info.GetResidueNumber()
            if residue_name == residue_key:
                healthy_residue_dict[residue_name].append(residue_number)

    return healthy_residue_dict

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
    with warnings.catch_warnings(): # pragma: no cover
        warnings.simplefilter('ignore', BiopythonWarning)
        structure = parser.get_structure('structure', input_file)[0]

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    
    pre_protein_handle = io.StringIO()
    protein_handle = io.StringIO()
    ligand_handle = io.StringIO()

    pre_residue_filter = PreResidueSelect(chain, include_metal)
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()
    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    out_format = data["protein_out"].split(".")[-1]
    residue_filter = ResidueSelect(chain, include_metal, out_format, healthy_residue_dict)

    pdbio.save(protein_handle, residue_filter)
    pdbio.save(ligand_handle, LigandSelect(ligand_id, chain=chain))

    if residue_filter.res_w_missing_atoms or residue_filter.res_w_incorrect_bond_length_angle:
        display_removed_residues(residue_filter)

    protein_handle.seek(0)
    ligand_handle.seek(0)

    extraction_data = {
        "protein": protein_handle.read(),
        "ligand": ligand_handle.read(),
        "res_w_missing_atoms": residue_filter.res_w_missing_atoms,
        "res_w_incorrect_bond_length_angle": residue_filter.res_w_incorrect_bond_length_angle
        }
    pre_protein_handle.close()
    protein_handle.close()
    ligand_handle.close()

    return extraction_data
