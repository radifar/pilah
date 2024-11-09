import json
import pathlib
import sys
from tempfile import NamedTemporaryFile

from meeko import (MoleculePreparation,
                   PDBQTReceptor,
                   PDBQTWriterLegacy)
from openbabel import openbabel as ob
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw


def check(success, error_msg): # pragma: no cover
    if not success:
        print("Error: " + error_msg, file=sys.stderr)
        sys.exit(2)

def renumber_hydrogens(mol):
    bonds = [bond for bond in mol.GetBonds()]
    interresidual_bond = []
    for index, bond in enumerate(bonds):
        atom_a_residue_number = bond.GetBeginAtom().GetPDBResidueInfo().GetResidueNumber()
        atom_b_residue_number = bond.GetEndAtom().GetPDBResidueInfo().GetResidueNumber()
        if atom_a_residue_number != atom_b_residue_number:
            interresidual_bond.append(index)
    
    fragments = Chem.FragmentOnBonds(mol, interresidual_bond, addDummies=False)
    mol_fragments = Chem.GetMolFrags(fragments, asMols=True)
    pdb_block = ""

    for fragment in mol_fragments:
        pdb_block += Chem.MolToPDBBlock(fragment, flavor=10).removesuffix("END\n")

    pdb_block += "END\n"

    return Chem.MolFromPDBBlock(pdb_block, removeHs=False)

def rename_hydrogens(mol):
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol != "H":
            hydrogens = []
            name = atom.GetMonomerInfo().GetName()
            residue_name = atom.GetPDBResidueInfo().GetResidueName()
            suffix = name.strip()[1:]
            suffix_len = len(suffix)
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == "H": hydrogens.append(neighbor)

            if len(hydrogens) == 1:
                # Treat hydrogen attached to C carboxyl whose neighbouring residue is missing
                # as hydrogen attached to C alpha (HA)
                # PDB ID 1e66 A 298 GLY
                if name.strip() == "C":
                    if residue_name == "GLY":
                        h_name_1 = "HA2"
                    else:
                        h_name_1 = "HA"
                # PDB ID 1e66 A 271 PRO
                elif name.strip() == "N" and residue_name == "PRO":
                    h_name_1 = "HA"
                else:
                    h_name_1 = "H" + suffix
                h_name_1 = f" {h_name_1:3}"
                hydrogens[0].GetMonomerInfo().SetName(h_name_1)
            elif len(hydrogens) == 2:
                # Treat hydrogen attached to N of peptide bond whose neighbouring residue is
                # missing as hydrogen attached to N of peptide bond (H).
                # This is important as following the common rule will name these hydrogen as
                # H1 and H2, which is did not matched with any residue template.
                if atom_symbol == "N":
                    if name.strip() == "N":
                        h_name_1 = "H"
                        h_name_2 = "H"
                    elif name.strip() == "NZ":
                        h_name_1 = "H" + suffix + "2"
                        h_name_2 = "H" + suffix + "3"
                    else:
                        h_name_1 = "H" + suffix + "1"
                        h_name_2 = "H" + suffix + "2"
                else:
                    h_name_1 = "H" + suffix + "2"
                    h_name_2 = "H" + suffix + "3"

                if suffix_len < 2:
                    h_name_1 = f" {h_name_1:3}"
                    h_name_2 = f" {h_name_2:3}"
                hydrogens[0].GetMonomerInfo().SetName(h_name_1)
                hydrogens[1].GetMonomerInfo().SetName(h_name_2)
            elif len(hydrogens) == 3:
                h_name_1 = "H" + suffix + "1"
                h_name_2 = "H" + suffix + "2"
                h_name_3 = "H" + suffix + "3"
                if suffix_len < 2:
                    h_name_1 = f" {h_name_1:3}"
                    h_name_2 = f" {h_name_2:3}"
                    h_name_3 = f" {h_name_3:3}"
                hydrogens[0].GetMonomerInfo().SetName(h_name_1)
                hydrogens[1].GetMonomerInfo().SetName(h_name_2)
                hydrogens[2].GetMonomerInfo().SetName(h_name_3)

    return mol

def mol_writer(mol, filename, ionization_records=None, receptor=False):
    mol_format = filename.split(".")[-1]

    if mol_format == "pdb":
        Chem.MolToPDBFile(mol, filename)
    elif mol_format == "mol2":
        pdb_block = Chem.MolToPDBBlock(mol)
        converter = ob.OBConversion()
        converter.SetInAndOutFormats("pdb", "mol2")
        mol = ob.OBMol()
        converter.ReadString(mol, pdb_block)
        converter.WriteFile(mol, filename)
    elif mol_format == "sdf":
        writer = Chem.SDWriter(filename)
        writer.write(mol)
        writer.close()
    elif mol_format == "pdbqt":
        pkg_dir = pathlib.Path(__file__).parents[0]
        with open(pkg_dir / "meeko" / "data" / "residue_params.json") as f:
            residue_params = json.load(f)
        if receptor:
            with NamedTemporaryFile() as temp_pdb_file:
                pdb_block = Chem.MolToPDBBlock(mol)
                pdb_block_renamed_residue = ""
                for line in pdb_block.splitlines():
                    if ("ATOM" not in line) and ("HETATM" not in line):
                        line += "\n"
                        pdb_block_renamed_residue += line
                        continue
                    chain_id = line[21].strip()
                    residue_number = int(line[22:26].strip())
                    residue_name = line[17:20].strip()
                    residue_id = (chain_id, residue_number, residue_name)
                    if residue_name == "LYS":
                        ionization_data = ionization_records[residue_id]
                        lys_charge = ionization_data[2]
                        if lys_charge == "Neutral":
                            line = line.replace("LYS", "LYN")
                    elif residue_name == "HIS":
                        ionization_data = ionization_records[residue_id]
                        his_charge = ionization_data[2]
                        if his_charge == "Neutral":
                            line = line.replace("HIS", "HIE")
                        elif his_charge == "Positive":
                            line = line.replace("HIS", "HIP")
                    elif residue_name == "CYS":
                        ionization_data = ionization_records[residue_id]
                        cys_charge = ionization_data[2]
                        if cys_charge == "Negative":
                            line = line.replace("CYS", "CYM")
                        elif cys_charge == "SS_bridge":
                            line = line.replace("CYS", "CYX")
                    elif residue_name == "ASP":
                        ionization_data = ionization_records[residue_id]
                        asp_charge = ionization_data[2]
                        if asp_charge == "Neutral":
                            line = line.replace("ASP", "ASH")
                    elif residue_name == "GLU":
                        ionization_data = ionization_records[residue_id]
                        glu_charge = ionization_data[2]
                        if glu_charge == "Neutral":
                            line = line.replace("GLU", "GLH")
                    line += "\n"
                    pdb_block_renamed_residue += line

                binary_pdb_block = pdb_block_renamed_residue.encode("ascii")
                temp_pdb_file.write(binary_pdb_block)
                # Chem.MolToPDBFile(mol, temp_pdb_file.name)
                receptor = PDBQTReceptor(temp_pdb_file.name, skip_typing=True)
                temp_pdb_file.close()

            is_ok, err = receptor.assign_types_charges(residue_params=residue_params)
            check(is_ok, err)

            all_flexres = set()
            pdbqt_string, is_ok, err = receptor.write_pdbqt_string(flexres=all_flexres)
            check(is_ok, err)
            if is_ok:
                with open(filename,'w') as w:
                    w.write(pdbqt_string["rigid"])
            else: # pragma: no cover 
                print(error_msg)
        else:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            for setup in mol_setups:
                pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
                if is_ok:
                    with open(filename,'w') as w:
                        w.write(pdbqt_string)
                else: # pragma: no cover 
                    print(error_msg)
    else:
        print(f"Unrecognized file format: {mol_format}")

def mol_drawer(mol, filename, size):
    image_format = filename.split(".")[-1]
    Chem.Compute2DCoords(mol)

    if image_format == "png":
        if size == "small":
            Draw.MolToFile(mol, filename=filename)
        elif size == "medium":
            Draw.MolToFile(mol, filename=filename, size=(420, 420))
        elif size == "large":
            Draw.MolToFile(mol, filename=filename, size=(640, 640))
        else:
            print(f"Unrecognized image size: {size}, using small size instead")
            Draw.MolToFile(mol, filename=filename)
    elif image_format == "svg":
        Draw.MolToFile(mol, filename=filename, imageType="svg")
    else:
        print(f"Unrecognized image format: {image_format}")
