from rdkit.Chem import AllChem as Chem

from pilah.dimorphite_dl import dimorphite_dl
from pilah.pKAI import pKAI
from pilah.element import positive_one, positive_two, positive_three


def process_ligand(data: dict, ligand_block: str):
    ligand_mol = Chem.MolFromPDBBlock(ligand_block)
    ligand_smiles = Chem.MolToSmiles(ligand_mol)
    true_smiles = data.get("ligand_smiles", ligand_smiles)
    pH = float(data.get("ph", 7.4))

    dimorphite_args = {
        "smiles": true_smiles,
        "min_ph": pH,
        "max_ph": pH,
        "pka_precision": 0
    }

    dimorphite_results = dimorphite_dl.Protonate(dimorphite_args)
    protonated_smiles = next(dimorphite_results).strip()
    template_mol = Chem.MolFromSmiles(protonated_smiles)
    corrected_ligand = Chem.AssignBondOrdersFromTemplate(template_mol, ligand_mol)
    corrected_ligand_with_Hs = Chem.AddHs(corrected_ligand,
                                          addCoords=True,
                                          addResidueInfo=True)
    
    return (corrected_ligand, corrected_ligand_with_Hs)

def get_atoms_from_pattern(mol, pattern):
    smarts = Chem.MolFromSmarts(pattern)
    return mol.GetSubstructMatches(smarts)

def process_protein(data: dict, protein_block: str):
    protein_mol = Chem.MolFromPDBBlock(protein_block)
    protein_with_Hs = Chem.AddHs(protein_mol,
                                 addCoords=True,
                                 addResidueInfo=True)
    
    model = data.get("pkai_model", "pKAI")
    pkai_result = pKAI.pKAI(protein_block, model_name=model)

    protein_without_Hs = Chem.RemoveHs(protein_with_Hs, updateExplicitCount=True)
    ionized_mol = AA_modifier(protein_without_Hs, pkai_result)

    pH = float(data.get("ph", 7.4))
    ptreshold = float(data.get("ptreshold", 1.0))
    
    ionized_mol.ionize_aa(pH, ptreshold)
    ionized_mol_with_Hs = ionized_mol.get_protonated_mol()
    ionization_records = ionized_mol.ionizable_AA

    return ionized_mol_with_Hs, ionization_records

class AA_modifier():
    def __init__(self, mol, pkai_result) -> None:
        self.mol = mol
        self.positive = positive_one + positive_two + positive_three
        ionizable_AA_Smarts = {
            "carboxylate": "C(=O)[OH]",
            "imidazole": "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1",
            "thiol": "[SH]",
            "phenol": "c[OH]",
            "guanidinium": "NC(=N)",
            "amonium": "[N;H2&+0][C;!$(C=*)]",
            # non-organic SMARTS pattern from the following discussion
            # https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/000001d83185$4e3aaff0$eab00fd0$@san.rr.com/
            "non-organic": "[!#5;!#6;!#7;!#8;!#16;!#15;!F;!Cl;!Br;!I;!#1]"
        }
        self.ionizable_AA = {}
        for aa in pkai_result:
            chain_id, residue_number, residue_name, pKa = aa
            residue_id = (chain_id, residue_number, residue_name)
            self.ionizable_AA[residue_id] = [pKa, set(), "Neutral"]
        
        self.residue_ids = self.ionizable_AA.keys()
        
        for pattern in ionizable_AA_Smarts.values():
            filtered_atoms = get_atoms_from_pattern(mol, pattern)
            for atom_ids in filtered_atoms:
                atom = mol.GetAtomWithIdx(atom_ids[0])
                residue = atom.GetPDBResidueInfo()
                chain_id = residue.GetChainId()
                residue_number = residue.GetResidueNumber()
                residue_name = residue.GetResidueName()

                residue_id = (chain_id, residue_number, residue_name)

                if residue_id in self.residue_ids:
                    self.ionizable_AA[residue_id][1].update(atom_ids)
                elif residue_name == "ARG":
                    self.ionizable_AA[residue_id] = [14.0, set(atom_ids), "Positive"]
                elif residue_name.strip() in self.positive:
                    self.ionizable_AA[residue_id] = [14.0, set(atom_ids), "Positive"]
    
    def ionize_aa(self, pH, ptreshold) -> None:
        for residue_id, ionization_data in self.ionizable_AA.items():
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
                # print(f"ARG with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "NH2":
                        atom.SetNumExplicitHs(2)
                        atom.SetFormalCharge(1)
            
            elif (residue_name == "HIS") and (pKa - ptreshold >=pH):
                ionization_data[2] = "Positive"
                # print(f"HIS with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "ND1":
                        atom.SetNumExplicitHs(1)
                        atom.SetFormalCharge(1)
            
            elif (residue_name == "CYS") and (pKa - ptreshold <= pH):
                ionization_data[2] = "Negative"
                # print(f"CYS with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "SG":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif (residue_name == "TYR") and (pKa - ptreshold <= pH):
                ionization_data[2] = "Negative"
                # print(f"TYR with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OH":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif (residue_name == "ASP") and (pKa - ptreshold <= pH):
                ionization_data[2] = "Negative"
                # print(f"ASP with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OD2":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif (residue_name == "GLU") and (pKa - ptreshold <= pH):
                ionization_data[2] = "Negative"
                # print(f"GLU with residue number {res_number}")
                for index in atom_idx_list:
                    atom = self.mol.GetAtomWithIdx(index)
                    atom_name = atom.GetMonomerInfo().GetName()
                    if atom_name.strip() == "OE2":
                        atom.SetNumExplicitHs(0)
                        atom.SetFormalCharge(-1)
            
            elif residue_name.strip() in positive_one:
                # print(f"positive one {residue_name} with residue number {res_number}")
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(1)

            elif residue_name.strip() in positive_two:
                # print(f"positive two {residue_name} with residue number {res_number}")
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(2)

            elif residue_name.strip() in positive_three:
                # print(f"positive three {residue_name} with residue number {res_number}")
                atom = self.mol.GetAtomWithIdx(atom_idx_list[0])
                atom.SetNumExplicitHs(0)
                atom.SetFormalCharge(3)

    
    def get_protonated_mol(self):
        protonated_mol = Chem.AddHs(self.mol,
                                    explicitOnly=True,
                                    addCoords=True,
                                    addResidueInfo=True)
        
        return protonated_mol
