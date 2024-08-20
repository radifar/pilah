from configparser import ConfigParser

from rich import print


options = [
    # mandatory
    "input",
    "ligand_id",
    "protein_out",
    "ligand_out",
    # optional
    "ligand_image",
    "image_size",
    "include_metal",
    "chain",
    "ligand_smiles",
    "pkai_model",
    "ph",
    "ptreshold"
]


class Config(ConfigParser):
    def __init__(self):
        super().__init__()
    
    def load(self, config_file):
        with open(config_file) as f:
            text = f.read()
        
        self.read_string("[pxpc]\n" + text)
        self.data = dict(self["pxpc"])

        for key in self.data.keys():
            if key not in options:
                print(f"Warning: '{key}' option is not recognized")
        