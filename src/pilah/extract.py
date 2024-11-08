from copy import deepcopy
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


console = Console()

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
    """
    After fixing insertion code (if exist), it will extract the protein only.
    This is required to strip any post translation modification prior to scanning
    the healthy residues.
    """
    def __init__(self, chain, include_metal="no"):
        self.chain_select = chain
        self.include_metal = include_metal

    def accept_chain(self, chain):
        return chain.id in self.chain_select
    
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return True
        elif self.include_metal == "yes":
            return residue.get_resname() in metal_list
        else:
            return False

class ResidueSelect(Select):
    def __init__(self,
                 chain,
                 include_metal="no",
                 out_format="",
                 healthy_residue_dict=dict(),
                 altloc_list=list()):
        self.chain_select = chain
        self.include_metal = include_metal
        self.out_format = out_format
        self.healthy_residue_dict = healthy_residue_dict
        self.altloc = dict()
        for altloc_id in altloc_list:
            value, res_id = altloc_id.split(":")
            self.altloc[res_id] = value
        self.res_w_missing_atoms = []
        self.res_w_incorrect_bond_length_angle = []

    def accept_chain(self, chain):
        return chain.id in self.chain_select

    def accept_residue(self, residue):
        # remove residues with missing atom if the output format is pdbqt
        if (self.out_format == "pdbqt") and residue.id[0] == " ":
            residue_atom_names = set()
            for atom in residue:
                atom_name = atom.get_name()
                residue_atom_names.add(atom_name)
            residue_name = residue.get_resname()
            chain_id = residue.get_full_id()[2]
            residue_number = residue.get_id()[1]
            residue_id = (chain_id, residue_number)
            residue_full_id = (chain_id, residue_name, residue_number)
            if residue_name not in metal_list:
                ref_atom_names = atom_names_dict[residue_name]
            else:
                ref_atom_names = set()
                ref_atom_names.add(residue_name)
            if (residue_atom_names == ref_atom_names) or (residue_atom_names == ref_atom_names.union({"OXT"})):
                if residue_name in ["GLY", "VAL", "ALA"]:
                    return True
                elif residue_name in metal_list:
                    if self.include_metal == "yes":
                        return True
                    else:
                        return False
                elif residue_id in self.healthy_residue_dict[residue_name]:
                    return True
                else:
                    self.res_w_incorrect_bond_length_angle.append(residue_full_id)
                    return False
            else:
                self.res_w_missing_atoms.append(residue_full_id)
                return False
        elif residue.id[0] == " ":
            return True
        elif self.include_metal == "yes":
            return residue.get_resname() in metal_list
        else:
            return False
    
    def accept_atom(self, atom):
        atom_altloc = atom.get_altloc()
        parent_residue_id = atom.get_parent().get_full_id()
        chain_id = parent_residue_id[2]
        res_number = parent_residue_id[3][1]
        res_name = atom.get_parent().get_resname()
        altloc_id = f"{chain_id}{res_number}{res_name}"
        if not atom.is_disordered():
            return True
        elif altloc_id in self.altloc.keys():
            altloc_match = atom_altloc == self.altloc[altloc_id]
            atom.set_altloc(" ")
            return altloc_match
        elif atom_altloc == "A": # disordered but not specified in config, use default altloc
            atom.set_altloc(" ")
            return True
        else:
            return False

class LigandSelect(Select):
    def __init__(self, res_name, chain, altloc_list=list()):
        self.res_name = res_name
        self.chain_select = chain
        self.ligand_res_num = []
        self.altloc = dict()
        for altloc_id in altloc_list:
            value, res_id = altloc_id.split(":")
            self.altloc[res_id] = value

    def accept_chain(self, chain):
        return chain.id == self.chain_select

    def accept_residue(self, residue):
        is_ligand = residue.get_resname() == self.res_name
        if is_ligand:
            self.ligand_res_num.append(residue.get_id()[1])
        return is_ligand
    
    def accept_atom(self, atom):
        atom_altloc = atom.get_altloc()
        parent_residue_id = atom.get_parent().get_full_id()
        chain_id = parent_residue_id[2]
        res_number = parent_residue_id[3][1]
        res_name = atom.get_parent().get_resname()
        altloc_id = f"{chain_id}{res_number}{res_name}"
        if not atom.is_disordered():
            return True
        elif altloc_id in self.altloc.keys():
            altloc_match = atom_altloc == self.altloc[altloc_id]
            atom.set_altloc(" ")
            return altloc_match
        elif atom_altloc == "A": # disordered but not specified in config, use default altloc
            atom.set_altloc(" ")
            return True
        else:
            return False

class LigandSelectByResNum(Select):
    def __init__(self, res_num):
        self.res_num = int(res_num)
    
    def accept_residue(self, residue):
        residue_number = residue.get_id()[1]
        return residue_number == self.res_num

def display_removed_residues(residue_filter):
    if residue_filter.res_w_missing_atoms:
        console.print("[bold deep_pink2]\n     Some residues with missing atoms were removed:[/bold deep_pink2]")
        console.print("[red3]     chain_id residue_name residue_number[/red3]")
        for chain_id, residue_name, residue_number in residue_filter.res_w_missing_atoms:
            console.print(f"[deep_pink1]         {chain_id}       {residue_name}          {residue_number}[/deep_pink1]")
    if residue_filter.res_w_incorrect_bond_length_angle:
        console.print("[bold deep_pink2]\n     Some residues with incorrect bond length and angle were removed:[/bold deep_pink2]")
        console.print("[red3]     chain_id residue_name residue_number[/red3]")
        for chain_id, residue_name, residue_number in residue_filter.res_w_incorrect_bond_length_angle:
            console.print(f"[deep_pink1]         {chain_id}       {residue_name}          {residue_number}[/deep_pink1]")

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
            chain_id = residue_info.GetChainId()
            residue_name = residue_info.GetResidueName()
            residue_number = residue_info.GetResidueNumber()
            residue_id = (chain_id, residue_number)
            if residue_name == residue_key:
                healthy_residue_dict[residue_name].append(residue_id)

    return healthy_residue_dict

def fix_insertion(structure, protein_chain, ligand_id):
    """
    This function will mutate the structure object by detaching every
    residue inside chain(s), then renumbering and delete the insertion code,
    finally it will reattach the residue to the structure object.
    This is because directly changing the id of a residue might cause error
    because the new id still belong to other residue.
    """
    renumber_residue_map = dict()
    total_insertion = dict()
    for chain in protein_chain:
        total_insertion[chain] = 0
        residue_list = []
        for residue in structure[chain].get_list():
            structure[chain].detach_child(residue.id)
            hetero, residue_number, ins_code = residue.id
            residue_name = residue.get_resname()
            old_id = (chain, residue_name, residue_number, ins_code)

            is_metal = residue_name in metal_list
            is_ligand = residue_name == ligand_id
            if hetero == " " or is_metal:
                if ins_code != " ":
                    total_insertion[chain] += 1
                residue_number += total_insertion[chain]
                residue.id = (" ", residue_number, " ")
                residue_list.append(residue)
                renumber_residue_map[old_id] = residue_number
            elif is_ligand:
                residue_number += total_insertion[chain]
                residue.id = (hetero, residue_number, " ")
                residue_list.append(residue)
                renumber_residue_map[old_id] = residue_number
        for residue in residue_list:
            structure[chain].add(residue)
    
    return total_insertion, renumber_residue_map

def extract(data):
    ligand_id = data["ligand_id"]
    protein_chain = data["protein_chain"].split()
    ligand_chain = data["ligand_chain"]
    include_metal = data.get("include_metal", "no")
    altloc = data.get("altloc", "")
    altloc_list = altloc.split()

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

    total_insertion, renumber_residue_map = fix_insertion(structure, protein_chain, ligand_id)
    if any(list(total_insertion.values())):
        console.print("[bold deep_pink1]\n     Insertion code detected, residue number on residues with insertion code[/bold deep_pink1]")
        console.print("[bold deep_pink1]     and the residues following those residues will be renumbered.[/bold deep_pink1]")
        console.print("[bold deep_pink1]     The mapping between original residue number and the new one will be logged into the Log file[/bold deep_pink1]")

    # create copy after residue renumbering
    structure_copy = deepcopy(structure)

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    
    pre_protein_handle = io.StringIO()
    protein_handle = io.StringIO()
    ligand_handle = io.StringIO()

    pre_residue_filter = PreResidueSelect(protein_chain, include_metal)
    pdbio.save(pre_protein_handle, pre_residue_filter)
    pre_protein_handle.seek(0)
    pre_pdb_block = pre_protein_handle.read()
    healthy_residue_dict = get_healthy_residues(pre_pdb_block)

    out_format = data["protein_out"].split(".")[-1]
    residue_filter = ResidueSelect(protein_chain,
                                   include_metal,
                                   out_format,
                                   healthy_residue_dict,
                                   altloc_list)

    pdbio.save(protein_handle, residue_filter)
    pdbio.save("test.pdb", residue_filter)
    if residue_filter.res_w_missing_atoms or residue_filter.res_w_incorrect_bond_length_angle:
        display_removed_residues(residue_filter)
    protein_handle.seek(0)

    ligand_filter = LigandSelect(ligand_id, ligand_chain, altloc_list)
    ligand_pdbio = PDBIO()
    ligand_pdbio.set_structure(structure_copy)
    ligand_pdbio.save(ligand_handle, ligand_filter)
    ligand_handle.seek(0)

    total_ligand = len(ligand_filter.ligand_res_num)
    res_num_not_provided = False
    if total_ligand > 1:
        lowest_ligand_res_num = ligand_filter.ligand_res_num[0]
        res_num_not_provided = "ligand_res_num" not in data.keys()
        selected_res_num = data.get("ligand_res_num", lowest_ligand_res_num)

        ligand_structure = parser.get_structure('ligand', ligand_handle)[0]
        ligand_pdbio = PDBIO()
        ligand_pdbio.set_structure(ligand_structure)
        ligand_handle.seek(0)
        ligand_handle.truncate()
        ligand_filter_by_resnum = LigandSelectByResNum(selected_res_num)
        ligand_pdbio.save(ligand_handle, ligand_filter_by_resnum)
        ligand_handle.seek(0)

        if res_num_not_provided:
            ligand_res_num = ', '.join(str(i) for i in ligand_filter.ligand_res_num)
            console.print(f"[bold deep_pink1]\n     Multiple ligand detected in chain {ligand_chain}. The default is to choose the ligand with the lowest residue number.[/bold deep_pink1]")
            console.print("[bold deep_pink1]     To choose ligand with specific residue number use 'ligand_seq_num' option.[/bold deep_pink1]")
            console.print(f"[bold deep_pink1]     The residue number of ligand detected in chain {ligand_chain} are: {ligand_res_num}.[/bold deep_pink1]")
    

    extraction_data = {
        "protein": protein_handle.read(),
        "ligand": ligand_handle.read(),
        "res_w_missing_atoms": residue_filter.res_w_missing_atoms,
        "res_w_incorrect_bond_length_angle": residue_filter.res_w_incorrect_bond_length_angle,
        "total_insertion": total_insertion,
        "renumber_residue_map": renumber_residue_map
        }
    
    if res_num_not_provided:
        extraction_data["multiple_ligand"] = [ligand_chain, ligand_id, ligand_res_num]
    
    pre_protein_handle.close()
    protein_handle.close()
    ligand_handle.close()

    return extraction_data
