"""Extract the Itzi output from gisdb and store them to GeoTIFF files"""
import sys
import modules.GRASS_utils as grutl
from modules import common
import yaml

def create_config_dictionary_from_config_file(config_filename):
    with open(config_filename) as file:
        config = yaml.full_load(file)
    return config
#* ---
#* We read the input information from the ini file and then the magic happens.
#*
#* This script exports the rasters generated from itzi simulation to .tif files
#* ---
def main(config_filename):
    config = create_config_dictionary_from_config_file(config_filename)

    mygisdb = config['grass_info']['grass_db']
    mylocation = config['grass_info']['location']
    mymapset = config['grass_info']['mapset']
    CRS = config['grass_info']['CRS']
    prefix = config['grass_output']['prefix']
    if not prefix:
        prefix = f'{mymapset}_itzi'

    user = grutl.initiate_GRASS_sesion(str(mygisdb), mylocation, mymapset, CRS_code=CRS)

    output_path = config['folders']['itzi_output_files']

    grutl.GRASS_export_rasters(str(output_path), mymapset, search_criteria = f'{prefix}*')

    grutl.end_GRASS_sesion(user)

    print('\nDONE\n')

if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)
