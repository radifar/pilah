import re

from PIL import Image
from data import atom_name_per_residue


ionizable_atom_dict = {
    # LYS
    "NZ": {
        "Neutral": {"Hyd": 2, "Charge": 0},
        "Positive": {"Hyd": 3, "Charge": 1}
    },
    # ARG
    "NH2": {
        "Positive": {"Hyd": 2, "Charge": 1}
    },
    # HIS
    "ND1": {
        "Neutral": {"Hyd": 0, "Charge": 0},
        "Positive": {"Hyd": 1, "Charge": 1}
    },
    # CYS
    "SG": {
        "Neutral": {"Hyd": 1, "Charge": 0},
        "Negative": {"Hyd": 0, "Charge": -1},
        "SS_bridge": {"Hyd": 0, "Charge": 0}
    },
    # TYR
    "OH": {
        "Neutral": {"Hyd": 1, "Charge": 0},
        "Negative": {"Hyd": 0, "Charge": -1}
    },
    # ASP
    "OD2": {
        "Neutral": {"Hyd": 1, "Charge": 0},
        "Negative": {"Hyd": 0, "Charge": -1}
    },
    # GLU
    "OE2": {
        "Neutral": {"Hyd": 1, "Charge": 0},
        "Negative": {"Hyd": 0, "Charge": -1}
    }
}

def assert_protonated_mol(protonated_mol, whole_i_records, ids):
    atom_names = ionizable_atom_dict.keys()
    for id in ids:
        _, atom_id_list, ionization_state = whole_i_records[id]

        for atom_id in atom_id_list:
            atom = protonated_mol.GetAtomWithIdx(atom_id)
            atom_name = atom.GetMonomerInfo().GetName().strip()
            if atom_name in atom_names:
                hyd_num = atom.GetTotalNumHs(includeNeighbors=True)
                charge = atom.GetFormalCharge()

                expected_hyd_num = ionizable_atom_dict[atom_name][ionization_state]["Hyd"]

                assert hyd_num == expected_hyd_num
                if ionization_state == "Neutral":
                    assert charge == 0
                elif ionization_state == "Positive":
                    assert charge == 1
                elif ionization_state == "Negative":
                    assert charge == -1

def assert_atom_names_in_residues(residue_atom_names_dict):
    ref_atom_names_dict = atom_name_per_residue.atom_names_dict
    resname_prefix = ["", "N", "C"]

    for residue_id, atom_names in residue_atom_names_dict.items():
        residue_flag = []
        residue_name = residue_id[0]
        for prefix in resname_prefix:
            full_res_name = prefix + residue_name
            ref_atom_names = ref_atom_names_dict[full_res_name]
            if atom_names <= ref_atom_names:
                num_atoms = len(atom_names)
                num_ref_atoms = len(ref_atom_names)
                num_atoms_diff = num_ref_atoms - num_atoms
                # allow slight differences in number of atoms
                # in case the neighbor is missing and replaced with hydrogen by RDKit
                if num_atoms_diff <= 1:
                    residue_flag.append(True)
                else:
                    residue_flag.append(False)
            else:
                residue_flag.append(False)

        assert True in residue_flag

def assert_image_svg(filename):
    # https://stackoverflow.com/questions/63419010/check-if-an-image-file-is-a-valid-svg-file-in-python
    SVG_R = r'(?:<\?xml\b[^>]*>[^<]*)?(?:<!--.*?-->[^<]*)*(?:<svg|<!DOCTYPE svg)\b'
    SVG_RE = re.compile(SVG_R, re.DOTALL)

    with open(filename) as f:
        file_contents = f.read()

    assert SVG_RE.match(file_contents) is not None

def assert_image_size(filename, size, real_size):
    if size == "small":
        with Image.open(filename) as img:
            assert img.verify() is None
    else:
        with Image.open(filename) as img:
            assert img.size == real_size

