import sys
from configparser import ConfigParser

import modules.GRASS_utils as grutl

#* ---
#* We read the input information from the ini file and then the magic happens.
#*
#* This script exports the rasters generated from itzi simulation to .tif files
#* ---

ini_config_file = sys.argv[1]

# instantiate
config = ConfigParser()

# parse existing file
config.read(ini_config_file)

mygisdb = config.get('grass_info', 'mygisdb')
mylocation = config.get('grass_info', 'mylocation')
mymapset = config.get('grass_info', 'mymapset')
CRS = config.getint('grass_info', 'CRS')
prefix = config.get('grass_output', 'prefix')
if not prefix:
    prefix = f'{mymapset}_itzi'

user = grutl.initiate_GRASS_sesion(mygisdb, mylocation, mymapset, CRS_code = CRS)

output_path = 'GRASS_itzi/itzi_files'

grutl.GRASS_export_rasters(output_path, mymapset, search_criteria = f'{prefix}*')

grutl.end_GRASS_sesion(user)

print('\nDONE\n')
