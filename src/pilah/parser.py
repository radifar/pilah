import sys
from configparser import ConfigParser

from rich import print

mandatory_opts = [
    "input",
    "ligand_id",
    "protein_out",
    "ligand_out",
    "protein_chain",
    "ligand_chain",
    "ligand_smiles",
]

optional_opts = [
    "altloc",
    "ligand_res_num",
    "ligand_image",
    "image_size",
    "include_metal",
    "pkai_model",
    "ph",
    "ptreshold",
]


class Config(ConfigParser):
    def __init__(self):
        super().__init__()

    def load(self, config_file):
        with open(config_file) as f:
            text = f.read()

        self.read_string("[pxpc]\n" + text)
        self.data = dict(self["pxpc"])

        data_keys = self.data.keys()
        for key in data_keys:
            options = mandatory_opts + optional_opts
            if key not in options:
                print(f"Warning: '{key}' option is not recognized")

        for option in mandatory_opts:
            if option not in data_keys:
                sys.exit(f"Missing mandatory option: {option}")
