import os
import sys
import subprocess
from configparser import ConfigParser
from pathlib import Path

# pylint: disable=import-error
import modules.GRASS_utils as grutl
# import some convenient GRASS GIS Python API parts
from grass.script import core as gcore
import grass.script as gscript
import grass.script.setup as gsetup
# import grass python libraries
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import temporal as t

#* ---
#* We read the input information from the ini file and then the magic happens
#* ---

ini_config_file = sys.argv[1]

# instantiate
config = ConfigParser()

# parse existing file
config.read(ini_config_file)
grass_info = config['grass_info']
grass_time = config['grass_time']
grass_input = config['grass_input']
grass_output = config['grass_output']
grass_statistics = config['grass_statistics']
grass_options = config['grass_options']
storage_conf = config['storage']

#* Create a rc file for grass if it doesn't exist yet
# If GISRC is set, assume valid path. Else default to grass's default: $HOME/.grass7/rc
gisrc = os.environ.get('GISRC', None)
if gisrc is None or not gisrc:
    raise RuntimeError('Environment variable GISRC is not set! Cannot run grass.')

gisrc = Path(gisrc).expanduser().resolve()

if not gisrc.exists():
    print(f'GISRC file {gisrc} does not exist! Creating one based on {ini_config_file}:')
    rcstr = (f"GISDBASE: {grass_info.get('mygisdb')}\n"
             f"LOCATION_NAME: {grass_info.get('mylocation')}\n"
             f"MAPSET: {grass_info.get('mymapset')}\n"
             "GUI: text\n")
    print(f'\n{rcstr}')
    with open(gisrc, 'w', encoding='utf-8') as f:
        f.write(rcstr)

#* ---
#* 1. We define the paths in which we would like to work out the simulation
#*     - mygisdb: str. Path of working location
#*     - mylocation: str. Name of the location (Ex. Helsinki)
#*     - mymapset: str. User that will be working on that location
#*
#* 2. We initiate the grass session with 'initiate_GRASS_sesion'. If sessions does
#*    not exists, it creates it. If it exists, it opens it and load the files
#* ---

mygisdb = grass_info.get('mygisdb', 'GRASS_itzi/grassdata')
mylocation = grass_info.get('mylocation')
mymapset = grass_info.get('mymapset')
CRS = grass_info.get('CRS')

mapset_path = Path(os.path.join(mygisdb, mylocation, mymapset))

if os.path.exists(mapset_path):
    subprocess.call(["rm", "-r", mapset_path])

user = grutl.initiate_GRASS_sesion(mygisdb, mylocation, mymapset, CRS_code = CRS)

#* ---
#* 3. Add all the relevant rasters for the simulation.
#*    For conviniency, the name of the rasters in GRASS will be the same as the
#*    .tif file names but without the extension (i.e. input = Helsinki_cropped.tif,
#*    output = Helsinki_cropped)
#* ---

rasters_path = Path(mygisdb)

g.list(flags = 'p', type = 'raster')

grutl.import_multiple_raster_files(rasters_path, search_criteria = '*.tif')

# Set the region to match the raster of interest
g.region(raster = f'DEM_cropped@{mymapset}')

# We would like to mask data that falls outside the boundaries for the simulation
r.mask(raster = f'DEM_cropped@{mymapset}')

# If you would like to check the imported files
g.list(flags = 'p', type = 'raster')

#* ---
#* 5. Create time and space dataset with rain data required for the simulation.
#* ---

stds = 'rain_minutely'

t.create(output = stds, semantictype = 'mean', title = 'Rain Rate', description = 'Rain rate data for itzi')

#* ---
#* 6. Import the minutely recorded radar rain events.
#* ---

rain_path = rasters_path / 'rain'

grutl.import_multiple_raster_files(rain_path, search_criteria = '*.tif')

g.list(flags = 'p', type = 'raster')

#* ---
#* 7. Register the minutely rain radar data in the time and space dataset created previously.
#*    As there are many files and the only way grass allows registering many files at once is with a
#*    .txt file indicating the name of the files, we will first create the file and then register.
#* ---

rain_txt_file = rain_path / 'rain_registering_data.txt'

start_time = grass_time.get('start_time')
increment_number = int(grass_time.get('increment_number'))
increment_unit = grass_time.get('increment_unit')

grutl.create_rain_raster_text_file(rain_path, rain_txt_file, search_criteria = '*.tif', start_time = start_time,
                                    increment_number = increment_number, increment_unit = increment_unit)

t.register(type = 'raster', input = stds, file = str(rain_txt_file))

#* ---
#* 8. We now create the izti configuration file to run the simulation. Here we set all the parameters we need
#*    and call the respective rasters and dataset we would like to use for the simulation. All the parameters
#*    must be strings
#* ---

output_itzi_file = Path(storage_conf.get('store_root')) / storage_conf.get('itzi_config_file')
grass_bin_path = grass_info.get('grass_bin_path')
grassdata_path = mygisdb
location = mylocation
mapset = mymapset
end_time = grass_time.get('end_time')
record_step = grass_time.get('record_step')
duration = grass_time.get('duration')
dem = grass_input.get('dem') or 'DEM_cropped'
friction = grass_input.get('friction') or 'friction'
rain = stds
start_h = grass_input.get('start_h', '')
start_y = grass_input.get('start_y', '')
inflow = grass_input.get('inflow', '')
bctype = grass_input.get('bctype', '')
bcval = grass_input.get('bcval', '')
infiltration = grass_input.get('infiltration')
if not infiltration:

    if config.getboolean('rain', 'infiltration_rate'):

        inf_stds = 'infiltration_minutely'

        t.create(output=inf_stds, semantictype='mean', title='Infiltration Rate',
                 description='Infiltration rate data for itzi')

        infiltration_path = rasters_path / 'infiltration'

        grutl.import_multiple_raster_files(infiltration_path, search_criteria = '*.tif')

        infiltration_txt_file = infiltration_path / 'infiltration_registering_data.txt'
        start_time = grass_time.get('start_time')
        increment_number = int(grass_time.get('increment_number'))
        increment_unit = grass_time.get('increment_unit')

        grutl.create_rain_raster_text_file(infiltration_path, infiltration_txt_file,
                                           search_criteria='*.tif', start_time=start_time,
                                           increment_number=increment_number,
                                           increment_unit=increment_unit)

        t.register(type = 'raster', input = inf_stds, file = str(infiltration_txt_file))

        infiltration = inf_stds

    else:

        infiltration = 'infiltration_0'

losses = grass_input.get('losses') or 'losses'

prefix = grass_output.get('prefix') or f'{mymapset}_itzi'

values = grass_output.get('values')

stats_file = grass_statistics.get('stats_file') or f'GRASS_itzi/{mymapset}_itzi.csv'

hmin = grass_options.get('hmin')
slmax = grass_options.get('slmax')
cfl = grass_options.get('cfl')
theta = grass_options.get('theta')
vrouting = grass_options.get('vrouting')
dtmax = grass_options.get('dtmax')
dtinf = grass_options.get('dtinf')

grutl.create_itzi_config_file(output_itzi_file, record_step=record_step, dem=dem, friction=friction,
                              grass_bin=grass_bin_path, grassdata=grassdata_path,
                              location=location, mapset=mapset, start_time=start_time,
                              end_time=end_time, duration=duration, start_h=start_h,start_y=start_y,
                              rain=rain, inflow=inflow, bctype=bctype, bcval=bcval,
                              infiltration=infiltration, losses=losses, prefix=prefix,
                              values=values, stats_file=stats_file, hmin=hmin, slmax=slmax, cfl=cfl,
                              theta=theta, vrouting=vrouting, dtmax=dtmax, dtinf=dtinf)

#* ---
#* 9. End the current grass session
#* ---

grutl.end_GRASS_sesion(user)

#* ---
#* X. Run the simulation in the terminal:
#*
#*     $ itzi run output_itzi_file
#*
#* XI. Export the rasters of interest to GeoTiff file so that they can be analysed separately using:
#*     'export_tif.py'
#* ---

print('\nGRASS has been set correctly. Ready to run itzi\n')
