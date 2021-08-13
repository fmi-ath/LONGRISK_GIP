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
import os, sys, glob
from configparser import ConfigParser
from subprocess import call

"""
We read the input information from the ini file and then the magic happens
"""

ini_config_file = sys.argv[1]

# instantiate
config = ConfigParser()

# parse existing file
config.read(ini_config_file)

"""
1. We define the paths in which we would like to work out the simulation
    - mygisdb: str. Path of working location
    - mylocation: str. Name of the location (Ex. Helsinki)
    - mymapset: str. User that will be working on that location

2. We initiate the grass session with 'initiate_GRASS_sesion'. If sessions does
   not exists, it creates it. If it exists, it opens it and load the files
"""

mygisdb = config.get('grass_info', 'mygisdb')
mylocation = config.get('grass_info', 'mylocation')
mymapset = config.get('grass_info', 'mymapset')
CRS = config.get('grass_info', 'CRS')

mapset_path = os.path.join(mygisdb, mylocation, mymapset)

if os.path.exists(mapset_path):

    call(["rm", "-r", mapset_path])

user = grutl.initiate_GRASS_sesion(mygisdb, mylocation, mymapset, CRS_code = CRS)

"""
3. Add all the relevant rasters for the simulation.
   For conviniency, the name of the rasters in GRASS will be the same as the
   .tif file names but without the extension (i.e. input = Helsinki_cropped.tif,
   output = Helsinki_cropped)
"""

rasters_path = 'GRASS_itzi/grassdata'

g.list(flags = 'p', type = 'raster')

grutl.import_multiple_raster_files(rasters_path, search_criteria = '*.tif')

# Set the region to match the raster of interest
g.region(raster = 'DEM_cropped' + '@' + mymapset)

# We would like to mask data that falls outside the boundaries for the simulation
r.mask(raster = 'DEM_cropped' + '@' + mymapset)

# If you would like to check the imported files
g.list(flags = 'p', type = 'raster')

"""
5. Create time and space dataset with rain data required for the simulation.
"""

stds = 'rain_minutely'

t.create(output = stds, semantictype = 'mean', title = 'Rain Rate', description = 'Rain rate data for itzi')

"""
6. Import the minutely recorded radar rain events.
"""

rain_path = 'GRASS_itzi/grassdata/rain'

grutl.import_multiple_raster_files(rain_path, search_criteria = '*.tif')

g.list(flags = 'p', type = 'raster')

"""
7. Register the minutely rain radar data in the time and space dataset created previously.
   As there are many files and the only way grass allows registering many files at once is with a
   .txt file indicating the name of the files, we will first create the file and then register.
"""

rain_txt_file = 'GRASS_itzi/grassdata/rain/rain_registering_data.txt'
start_time = config.get('grass_time', 'start_time')
increment_number = config.getint('grass_time', 'increment_number')
increment_unit = config.get('grass_time', 'increment_unit')

grutl.create_rain_raster_text_file(rain_path, rain_txt_file, search_criteria = '*.tif', start_time = start_time,
                                    increment_number = increment_number, increment_unit = increment_unit)

t.register(type = 'raster', input = stds, file = rain_txt_file)

"""
8. We now create the izti configuration file to run the simulation. Here we set all the parameters we need
   and call the respective rasters and dataset we would like to use for the simulation. All the parameters
   must be strings
"""

output_itzi_file = 'GRASS_itzi/itzi_config_file.ini'
grass_bin_path = config.get('grass_info', 'grass_bin_path')
grassdata_path = mygisdb
location = mylocation
mapset = mymapset
end_time = config.get('grass_time', 'end_time')
record_step = config.get('grass_time', 'record_step')
duration = config.get('grass_time', 'duration')
dem = config.get('grass_input', 'dem')
if dem == '':
    dem = 'DEM_cropped'
friction = config.get('grass_input', 'friction')
if friction == '':
    friction = 'friction'
rain = stds
start_h = config.get('grass_input', 'start_h')
start_y = config.get('grass_input', 'start_y')
inflow = config.get('grass_input', 'inflow')
bctype = config.get('grass_input', 'bctype')
bcval = config.get('grass_input', 'bcval')
infiltration = config.get('grass_input', 'infiltration')
if infiltration == '':

    if config.getboolean('rain', 'infiltration_rate'):

        inf_stds = 'infiltration_minutely'

        t.create(output = inf_stds, semantictype = 'mean', title = 'Infiltration Rate', description = 'Infiltration rate data for itzi')

        infiltration_path = 'GRASS_itzi/grassdata/infiltration'

        grutl.import_multiple_raster_files(infiltration_path, search_criteria = '*.tif')

        infiltration_txt_file = 'GRASS_itzi/grassdata/infiltration/infiltration_registering_data.txt'
        start_time = config.get('grass_time', 'start_time')
        increment_number = config.getint('grass_time', 'increment_number')
        increment_unit = config.get('grass_time', 'increment_unit')

        grutl.create_rain_raster_text_file(infiltration_path, infiltration_txt_file, search_criteria = '*.tif', start_time = start_time,
                                            increment_number = increment_number, increment_unit = increment_unit)

        t.register(type = 'raster', input = inf_stds, file = infiltration_txt_file)

        infiltration = inf_stds

    else:

        infiltration = 'infiltration_0'

losses = config.get('grass_input', 'losses')
if losses == '':
    losses = 'losses'
prefix = config.get('grass_output', 'prefix')
if prefix == '':
    prefix = mymapset + '_itzi'
values = config.get('grass_output', 'values')
stats_file = config.get('grass_statistics', 'stats_file')
if stats_file == '':
    stats_file = 'GRASS_itzi/' + mymapset + '_itzi.csv'
hmin = config.get('grass_options', 'hmin')
slmax = config.get('grass_options', 'slmax')
cfl = config.get('grass_options', 'cfl')
theta = config.get('grass_options', 'theta')
vrouting = config.get('grass_options', 'vrouting')
dtmax = config.get('grass_options', 'dtmax')
dtinf = config.get('grass_options', 'dtinf')

grutl.create_itzi_config_file(output_itzi_file, grass_bin = grass_bin_path, grassdata = grassdata_path, location = location, mapset = mapset,
                            start_time = start_time, end_time = end_time, duration = duration, record_step = record_step, dem = dem, friction = friction,
                            start_h = start_h, start_y = start_y, rain = rain, inflow = inflow, bctype = bctype,
                            bcval = bcval, infiltration = infiltration, losses = losses, prefix = prefix, values = values,
                            stats_file = stats_file, hmin = hmin, slmax = slmax, cfl = cfl, theta = theta, vrouting = vrouting,
                            dtmax = dtmax, dtinf = dtinf)

"""
9. End the current grass session
"""

grutl.end_GRASS_sesion(user)

"""
X. Run the simulation in the terminal:

    $ itzi run output_itzi_file

XI. Export the rasters of interest to GeoTiff file so that they can be analysed separately using:
    'export_tif.py'
"""

print('\nGRASS has been set correctly. Ready to run itzi\n')
