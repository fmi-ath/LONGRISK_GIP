"""Extract the Itzi output from gisdb and store them to GeoTIFF files"""

import modules.GRASS_utils as grutl
from modules import common

#* ---
#* We read the input information from the ini file and then the magic happens.
#*
#* This script exports the rasters generated from itzi simulation to .tif files
#* ---

config = common.CONFIG

mygisdb = common.get_path_for('grass_db')
mylocation = config.get('grass_info', 'location')
mymapset = config.get('grass_info', 'mapset')
CRS = config.get('grass_info', 'CRS')
prefix = config.get('grass_output', 'prefix')
if not prefix:
    prefix = f'{mymapset}_itzi'

user = grutl.initiate_GRASS_sesion(str(mygisdb), mylocation, mymapset, CRS_code=CRS)

output_path = common.get_path_for('itzi_output')

grutl.GRASS_export_rasters(str(output_path), mymapset, search_criteria = f'{prefix}*')

grutl.end_GRASS_sesion(user)

print('\nDONE\n')
