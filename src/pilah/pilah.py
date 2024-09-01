import typer
from rich.console import Console

from pilah import __version__ as pilah_version
from pilah.extract import extract
from pilah.logger import log_writer
from pilah.parser import Config
from pilah.protonate import process_ligand, process_protein
from pilah.writer import mol_writer, mol_drawer, renumber_hydrogens, rename_hydrogens


console = Console()
app = typer.Typer()

@app.command(short_help="run PDB extraction, protonation and correction")
def run(config_file: str):
    console.print(f"[bold deep_pink3]   Running {config_file}[/bold deep_pink3]")
    print(f"   on PiLAH version: {pilah_version}")
    
    config = Config()
    config.load(config_file)
    
    complex_pdb_block = extract(config.data)
    
    ligand_pdb_block = complex_pdb_block["ligand"]
    protein_pdb_block = complex_pdb_block["protein"]
    ligand_mol, ligand_with_Hs = process_ligand(config.data, ligand_pdb_block)

    
    protein_protonation_results = process_protein(config.data, protein_pdb_block)
    protein_with_Hs, ionization_records = protein_protonation_results
    protein_with_Hs_renumbered = renumber_hydrogens(protein_with_Hs)
    protein_with_Hs_renamed = rename_hydrogens(protein_with_Hs_renumbered)

    ligand_id = config.data["ligand_id"]
    ligand_out = config.data.get("ligand_out", f"{ligand_id}_out.pdb")
    protein_out = config.data.get("protein_out")

    mol_writer(ligand_with_Hs, ligand_out)
    mol_writer(protein_with_Hs_renamed, protein_out, ionization_records, receptor=True)

    log_writer(config.data, pilah_version, ionization_records)

    if "ligand_image" in config.data.keys():
        ligand_image = config.data["ligand_image"]

        if "image_size" in config.data.keys():
            image_size = config.data["image_size"]
        else:
            image_size = "small"

        mol_drawer(ligand_mol, ligand_image, image_size)

# one command one callback, just a temporary helper command
# delete this when there are multiple commands
# https://typer.tiangolo.com/tutorial/commands/one-or-multiple/
@app.callback()
def callback():
    pass

if __name__ == "__main__": # pragma: no cover
    app()
