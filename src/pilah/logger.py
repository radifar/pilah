from datetime import datetime


def log_writer(
    config_data,
    extraction_data,
    ionization_records,
    ligand_processing_log,
    pilah_version,
):  # pragma: no cover
    res_w_missing_atoms = extraction_data["res_w_missing_atoms"]
    res_w_incorrect_bond_length_angle = extraction_data["res_w_incorrect_bond_length_angle"]
    total_insertion = extraction_data["total_insertion"]
    renumber_residue_map = extraction_data["renumber_residue_map"]
    multiple_ligand = extraction_data.get("multiple_ligand", "")
    force_remove_hyd = ligand_processing_log.get("force_remove_hyd", "")
    ligand_missing_atoms = ligand_processing_log.get("ligand_missing_atoms", "")

    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_txt = f"PiLAH version: {pilah_version}\n\n"
    log_txt += "----- PiLAH Configuration -----\n\n"
    logfile_name = f"log_{now}.txt"
    with open(logfile_name, "w") as log_handle:
        for key, value in config_data.items():
            log_txt += f"{key:16}: {value}\n"

        if multiple_ligand:
            ligand_chain, ligand_id, ligand_res_num = multiple_ligand
            log_txt += "\n\n----- Multiple Ligand Detected -----"
            log_txt += (
                f"\n\nThe residue number of ligand with id {ligand_id} in chain {ligand_chain}:\n{ligand_res_num}"
            )
            log_txt += "\n\nBy default it will choose the ligand with the smallest residue number."

        if force_remove_hyd:
            log_txt += "\n\n----- Force Remove Hydrogens of Ligand -----\n"
            log_txt += force_remove_hyd

        if ligand_missing_atoms:
            log_txt += "\n\n----- Ligand Atoms Are Missing -----\n"
            log_txt += ligand_missing_atoms

        log_txt += "\n\n----- Ionization Records -----\n"
        log_txt += "\nChain   Residue Number   Residue Name     pKa  Charge\n"

        ionization_records = sort_records(ionization_records)
        for residue_id, record in ionization_records.items():
            chain, resnum, resname = residue_id
            pKa, _, charge = record
            log_txt += f"  {chain}  {resnum:11} {'':12} {resname:9} {pKa:5}  {charge}\n"

        if res_w_missing_atoms or res_w_incorrect_bond_length_angle:
            log_txt += "\n\n----- Removed Residues -----"
        if res_w_missing_atoms:
            log_txt += "\n\nSome residues with missing atoms were removed:"
            log_txt += "\nchain_id residue_name residue_number"
            for chain_id, residue_name, residue_number in res_w_missing_atoms:
                log_txt += f"\n     {chain_id}       {residue_name}          {residue_number}"
        if res_w_incorrect_bond_length_angle:
            log_txt += "\n\nSome residues with incorrect bond length and angle were removed:"
            log_txt += "\nchain_id residue_name residue_number"
            for (
                chain_id,
                residue_name,
                residue_number,
            ) in res_w_incorrect_bond_length_angle:
                log_txt += f"\n     {chain_id}       {residue_name}          {residue_number}"

        if any(list(total_insertion.values())):
            log_txt += "\n\n\n----- Renumbered Residues -----"
            for chain, total_num in total_insertion.items():
                log_txt += f"\n\nTotal residue with insertion code in chain {chain} = {total_num}"
            log_txt += "\n\nThe mapping is formatted as follows:\n chain_id residue_name residue_number insertion_code => residue_name new_residue_number\n"
            for old_id, new_resnum in renumber_residue_map.items():
                log_txt += f"\n    {old_id[0]} {old_id[1]} {old_id[2]:4} {old_id[3]} => {old_id[1]} {new_resnum:4}"

        log_handle.write(log_txt)

    print(f"\n   Log file written to {logfile_name}\n")


def sort_records(records):
    temp_records = dict()
    sorted_records = dict()

    for residue_id, record in records.items():
        chain, resnum, _ = residue_id
        temp_records[(chain, resnum)] = [residue_id, record]

    temp_records_list = sorted(temp_records.items(), key=lambda x: (x[0][0], x[0][1]))

    for _, record_list in temp_records_list:
        residue_id, record = record_list
        sorted_records[residue_id] = record

    return sorted_records
