from datetime import datetime

def log_writer(config_data, pilah_version, ionization_records): # pragma: no cover
    now = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_txt = f"PiLAH version: {pilah_version}\n\n"
    log_txt += "----- PiLAH Configuration -----\n\n"
    logfile_name = f"log_{now}.txt"
    with open(logfile_name, "w") as log_handle:
        for key, value in config_data.items():
            log_txt += f"{key:16}: {value}\n"
        
        log_txt += "\n\n----- Ionization records -----\n"
        log_txt += "\nChain   Residue Number   Residue Name     pKa  Charge\n"

        ionization_records = sort_records(ionization_records)
        for residue_id, record in ionization_records.items():
            chain, resnum, resname = residue_id
            pKa, _, charge = record
            log_txt += f"  {chain}  {resnum:11} {'':12} {resname:9} {pKa:5}  {charge}\n"

        log_handle.write(log_txt)
    
    print(f"\n   Log file written to {logfile_name}\n")

def sort_records(records):
    temp_records = dict()
    sorted_records = dict()

    for residue_id, record in records.items():
        _, resnum, _ = residue_id
        temp_records[resnum] = [residue_id, record]
    
    temp_records_list = sorted(temp_records.items())

    for resnum, record_list in temp_records_list:
        residue_id, record = record_list
        sorted_records[residue_id] = record

    return sorted_records

    
