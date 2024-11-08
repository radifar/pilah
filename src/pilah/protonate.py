from rdkit.Chem import AllChem as Chem
from rich.console import Console

from pilah.dimorphite_dl import dimorphite_dl
from pilah.pKAI import pKAI
from pilah.element import positive_one, positive_two, positive_three


def process_ligand(data: dict, ligand_block: str):
    """
    Provide RDKit mol object in corrected form using provided SMILES template (optional).
    When template not provided, the ligand will be provided as is.
    It will return ligand in protonated and unprotonated form (for image generation)
    """
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)
    ligand_smiles = data["ligand_smiles"]
    pH = float(data.get("ph", 7.4))
    ligand_processing_log = dict()

    dimorphite_args = {
        "smiles": ligand_smiles,
        "min_ph": pH,
        "max_ph": pH,
        "pka_precision": 0
    }

    dimorphite_results = dimorphite_dl.Protonate(dimorphite_args)
    protonated_smiles = next(dimorphite_results).strip()
    template_mol = Chem.MolFromSmiles(protonated_smiles)
    try:
        corrected_ligand = Chem.AssignBondOrdersFromTemplate(template_mol, ligand_mol)
    except ValueError:
        console = Console()
        warning = """\n     Ligand bond assignment failed, trying to force remove
     hydrogens from the template before bond assignment.
     Notice that it might cause incorrect stereochemistry."""
        console.print(f"[bold deep_pink2]{warning}[/bold deep_pink2]")
        template_mol = Chem.RemoveAllHs(template_mol)
        corrected_ligand = Chem.AssignBondOrdersFromTemplate(template_mol, ligand_mol)
        ligand_processing_log["force_remove_hyd"] = warning
    corrected_ligand_with_Hs = Chem.AddHs(corrected_ligand,
                                          addCoords=True,
                                          addResidueInfo=True)
    
    return (corrected_ligand, corrected_ligand_with_Hs, ligand_processing_log)

def get_atoms_from_pattern(mol, pattern):
    smarts = Chem.MolFromSmarts(pattern)
    
    return mol.GetSubstructMatches(smarts)

def process_protein(data: dict, protein_block: str):
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    protein_with_Hs = Chem.AddHs(protein_mol,
                                 addCoords=True,
                                 addResidueInfo=True)
    
    protein_output_format = data["protein_out"].split(".")[-1]
    model = data.get("pkai_model", "pKAI")
    pkai_result = pKAI.pKAI(protein_block, model_name=model)

    # The hydrogenization then dehydrogenization is necessary
    # because it requires correct explicitCount for later
    # hydrogenization with explicitOnly option
    protein_without_Hs = Chem.RemoveHs(protein_with_Hs, updateExplicitCount=True)
    ionized_mol = AA_modifier(protein_without_Hs, pkai_result)

    pH = float(data.get("ph", 7.4))
    ptreshold = float(data.get("ptreshold", 1.0))

    if protein_output_format == "pdbqt":
        ionized_mol.no_tyr_ionization = True
    ionized_mol.ionize_aa(pH, ptreshold)
    ionized_mol_with_Hs = ionized_mol.get_protonated_mol()
    ionization_records = ionized_mol.ionization_records

    return ionized_mol_with_Hs, ionization_records

ionizable_AA_Smarts = {
    "carboxylate": "C(=O)[OH]",
    "imidazole": "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1",
    "thiol": "[SH]",
    "disulfide": "[#16X2H0][#16X2H0]",
    "phenol": "c[OH]",
    "guanidinium": "NC(=N)",
    "amonium": "[N;H2&+0][C;!$(C=*)]",
    # non-organic SMARTS pattern from the following discussion
    # https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/000001d83185$4e3aaff0$eab00fd0$@san.rr.com/
    "non-organic": "[!#5;!#6;!#7;!#8;!#16;!#15;!F;!Cl;!Br;!I;!#1]"
}

class AA_modifier():
    def __init__(self, mol, pkai_result) -> None:
        self.mol = mol
        self.positive = positive_one + positive_two + positive_three
        # pdbqt can not recognize deprotonated TYR, when output format is pdbqt
        # self.no_tyr_ionization will set to True
        self.no_tyr_ionization = False
        
        self.ionization_records = {}
        for aa in pkai_result:
            chain_id, residue_number, residue_name, pKa = aa
            residue_id = (chain_id, residue_number, residue_name)
            self.ionization_records[residue_id] = [pKa, set(), "Neutral"]
        
        self.residue_ids = self.ionization_records.keys()
        
        for moiety, pattern in ionizable_AA_Smarts.items():
            filtered_atoms = get_atoms_from_pattern(mol, pattern)
            if moiety == "disulfide":
                for atom_ids in filtered_atoms:
                    for atom_id in atom_ids:
                        atom = mol.GetAtomWithIdx(atom_id)
                        residue = atom.GetPDBResidueInfo()
                        chain_id = residue.GetChainId()
                        residue_number = residue.GetResidueNumber()
                        residue_name = residue.GetResidueName()

                        residue_id = (chain_id, residue_number, residue_name)
                        self.ionization_records[residue_id][1].update([atom_id])
                        self.ionization_records[residue_id][2] = "SS_bridge"
            else:
                for atom_ids in filtered_atoms:
                    atom = mol.GetAtomWithIdx(atom_ids[0])
                    residue = atom.GetPDBResidueInfo()
                    chain_id = residue.GetChainId()
                    residue_number = residue.GetResidueNumber()
                    residue_name = residue.GetResidueName()

                    residue_id = (chain_id, residue_number, residue_name)

                    if residue_id in self.residue_ids:
                        self.ionization_records[residue_id][1].update(atom_ids)
                    elif residue_name == "ARG":
                        self.ionization_records[residue_id] = [14.0, set(atom_ids), "Positive"]
                    elif residue_name.strip() in self.positive:
                        self.ionization_records[residue_id] = [14.0, set(atom_ids), "Positive"]
            
            # terminal carboxylate always deprotonated to match Meeko template
            if moiety == "carboxylate":
                atom_loop_break = False
                for atom_ids in filtered_atoms:
                    for atom_id in atom_ids:
                        atom = self.mol.GetAtomWithIdx(atom_id)
                        atom_name = atom.GetMonomerInfo().GetName()
                        if atom_name.strip() == "OXT":
                            atom.SetNumExplicitHs(0)
                            atom.SetFormalCharge(-1)
                            
                            atom_loop_break = True
                            break
                    if atom_loop_break:
                        break
    
    def ionize_aa(self, pH, ptreshold) -> None:
        for residue_id, ionization_data in self.ionization_records.items():
            residue_name = residue_id[2]
            pKa, atom_idx, ionization_state = ionization_data
            atom_idx_list = sorted(atom_idx)
            
            if (residue_name == "LYS") and (pKa - ptreshold >= pH):
                ionization_data[2] = "Positive"

                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "NZ":
                        atom.SetNumExplicitHs(3)
                        atom.SetFormalCharge(1)
            
            elif (residue_name == "ARG") and (pKa - ptreshold >=pH):
                ionization_data[2] = "Positive"
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "NH2":
                        atom.SetNumExplicitHs(2)
                        atom.SetFormalCharge(1)
            
            elif (residue_name == "HIS") and (pKa - ptreshold >=pH):
                ionization_data[2] = "Positive"
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "ND1":
                        atom.SetNumExplicitHs(1)
                        atom.SetFormalCharge(1)
            
            elif residue_name == "CYS":
                if ionization_data[2] == "SS_bridge":
                    continue
                elif (pKa + ptreshold <= pH):
                    ionization_data[2] = "Negative"
                    for index in atom_idx_list:
                        atom = self.mol.GetAtomWithIdx(index)
                        atom_name = atom.GetMonomerInfo().GetName()
                        if atom_name.strip() == "SG":
                            atom.SetNumExplicitHs(0)
                            atom.SetFormalCharge(-1)
            
            elif (residue_name == "TYR") and (pKa + ptreshold <= pH):
                if self.no_tyr_ionization:
                    continue
                ionization_data[2] = "Negative"
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OH":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif (residue_name == "ASP") and (pKa + ptreshold <= pH):
                ionization_data[2] = "Negative"
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OD2":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif (residue_name == "GLU") and (pKa + ptreshold <= pH):
                ionization_data[2] = "Negative"
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OE2":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif residue_name.strip() in positive_one:
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(1)

            elif residue_name.strip() in positive_two:
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(2)

            elif residue_name.strip() in positive_three:
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(3)

    
    def get_protonated_mol(self):
        protonated_mol = Chem.AddHs(self.mol,
                                    explicitOnly=True,
                                    addCoords=True,
                                    addResidueInfo=True)
        
        return protonated_mol

