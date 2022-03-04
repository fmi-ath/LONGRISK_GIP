"""Common utility functions. This file does not import any other modules from
this folder.
"""
from pathlib import Path
import yaml

def create_config_dictionary_from_config_file(config_filename):
    with open(config_filename) as file:
        config = yaml.full_load(file)
    return config

def ensure_folders_exist(config) -> None:
    """Create output folders if they don't exist"""
    for value in (('folders','temporary_files'), ('folders','processed_files'), ('folders','itzi_output_files'), ('grass_info','grass_db')):
        p = Path(config[value[0]][value[1]])
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)