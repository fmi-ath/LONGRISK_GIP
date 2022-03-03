"""Common utility functions. This file does not import any other modules from
this folder.
"""
from pathlib import Path
import yaml

def create_config_dictionary_from_config_file(config_filename):
    with open(config_filename) as file:
        config = yaml.full_load(file)
    return config

def get_path_for(p: str, config) -> Path :
    """Helper function to get path settings from configuration (*_files parameters)

    Example: if p=temporary, return value for temporary_files

    Args:
        p (str): tag for the configuration parameter (that ends with _files) or
                 `mygisdb` or `grass_db` for GRASS database path

    Returns:
        Path: value from configuration file
    """
    if p in {'mygisdb', 'grass_db'}:
        return Path(config['grass_info']['grass_db'])
    path = config['folders'][f'{p}_files']
    return Path(path)

def ensure_folders_exist(config) -> None:
    """Create output folders if they don't exist"""
    for value in ('temporary_files', 'processed_files', 'itzi_output_files', 'grass_db_files'):
        p = Path(config['folders'][value])
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)