"""Common utility functions. This file does not import any other modules from
this folder.
"""
import sys
from configparser import ConfigParser
from pathlib import Path

# SETUP: read config file to provide values for other scripts
try:
    config_file_name = sys.argv[1]
except IndexError as e:
    raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e

CONFIG = ConfigParser()
CONFIG.read(config_file_name)

def get_path_for(p: str) -> Path :
    """[summary]

    Args:
        p (str): [description]

    Returns:
        Path: [description]
    """
    if p in {'mygisdb', 'grass_db'}:
        return Path(CONFIG.get('grass_info', 'mygisdb'))
    path = CONFIG.get('folders', f'{p}_files')
    return Path(path)

def ensure_folders_exist() -> None:
    """Create output folders if they don't exist"""
    for value in ('temporary', 'grass_output', 'itzi_output', 'grass_db'):
        p = get_path_for(value)
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':
    ensure_folders_exist()
