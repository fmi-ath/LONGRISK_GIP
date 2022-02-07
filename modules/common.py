"""Common utility functions. This file does not import any other modules from
this folder.
"""
import sys
from configparser import ConfigParser
from pathlib import Path




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
        return Path(config.get('grass_info', 'grass_db'))
    path = config.get('folders', f'{p}_files')
    return Path(path)

def ensure_folders_exist(config) -> None:
    """Create output folders if they don't exist"""
    for value in ('temporary', 'processed', 'itzi_output', 'grass_db'):
        p = get_path_for(value, config)
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':
    # SETUP: read config file to provide values for other scripts
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e

    CONFIG = ConfigParser()
    CONFIG.read(config_file_name)
    ensure_folders_exist(CONFIG)
