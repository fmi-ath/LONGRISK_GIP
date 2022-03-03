"""Extract the Itzi output from gisdb and store them to GeoTIFF files"""
import sys
import modules.GRASS_utils as grutl
import yaml
import modules.common as cm
#* This script exports the rasters generated from itzi simulation to .tif files

def main(config_filename):
    config = cm.create_config_dictionary_from_config_file(config_filename)
    mymapset = config['grass_info']['mapset']
    prefix = config['grass_output']['prefix']
    if not prefix:
        prefix = f'{mymapset}_itzi'

    user = grutl.initiate_GRASS_session(config['grass_info'])
    grutl.GRASS_export_rasters(str(config['folders']['itzi_output_files']), mymapset, search_criteria = f'{prefix}*')
    grutl.end_GRASS_session(user)
    print('\nDONE\n')

if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)
