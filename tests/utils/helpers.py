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
        "Negative": {"Hyd": 0, "Charge": -1}
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
    # TODO: change ionizable_AA to ionization records so that it can be reused by
    # test_process_protein
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
