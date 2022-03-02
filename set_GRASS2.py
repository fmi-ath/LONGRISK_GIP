import os
import subprocess
from pathlib import Path
import sys
# import some convenient GRASS GIS Python API parts
from grass.script import core as gcore
import grass.script as gscript
import grass.script.setup as gsetup
# import grass python libraries
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import temporal as t

import modules.GRASS_utils as grutl
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
    for value in ('temporary', 'processed', 'itzi_output', 'grass_db'):
        p = get_path_for(value, config)
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)

def main(config_filename):
    #* ---
    #* We read the input information from the ini file and then the magic happens
    #* ---

    config = create_config_dictionary_from_config_file(config_filename)
    # parse existing file
    grass_info = config['grass_info']
    grass_time = config['grass_time']
    grass_input = config['grass_input']
    grass_output = config['grass_output']
    grass_statistics = config['grass_statistics']
    grass_options = config['grass_options']

    #* Create a rc file for grass if it doesn't exist yet
    # If GISRC is set, assume valid path. Else default to grass's default: $HOME/.grass7/rc
    gisrc = os.environ.get('GISRC', None)
    if gisrc is None or not gisrc:
        raise RuntimeError('Environment variable GISRC is not set! Cannot run grass.')

    gisrc = Path(gisrc).expanduser().resolve()

    if not gisrc.exists():
        print(f'GISRC file {gisrc} does not exist! Creating one based on config file:')
        rcstr = (f"GISDBASE: {grass_info['grass_db']}\n"
                f"LOCATION_NAME: {grass_info['location']}\n"
                f"MAPSET: {grass_info['mapset']}\n"
                "GUI: text\n")
        print(f'\n{rcstr}')
        with open(gisrc, 'w', encoding='utf-8') as f:
            f.write(rcstr)

    #* ---
    #* 1. We define the paths in which we would like to work out the simulation
    #*     - grass_db: str. Path of working location
    #*     - location: str. Name of the location (Ex. Helsinki)
    #*     - mapset: str. User that will be working on that location
    #*
    #* 2. We initiate the grass session with 'initiate_GRASS_session'. If sessions does
    #*    not exists, it creates it. If it exists, it opens it and load the files
    #* ---

    grassdata_path = get_path_for('grass_db',config)
    grass_info['grassdata_path'] = grassdata_path
    location = grass_info['location']
    mapset = grass_info['mapset']
    CRS = grass_info['CRS']

    mapset_path = Path(os.path.join(grassdata_path, location, mapset))

    if mapset_path.exists():
        subprocess.call(["rm", "-r", str(mapset_path)])

    user = grutl.initiate_GRASS_session(str(grassdata_path), location, mapset, CRS_code=CRS)

    #* ---
    #* 3. Add all the relevant rasters for the simulation.
    #*    For conviniency, the name of the rasters in GRASS will be the same as the
    #*    .tif file names but without the extension (i.e. input = Helsinki_cropped.tif,
    #*    output = Helsinki_cropped)
    #* ---

    rasters_path = Path(grassdata_path)

    g.list(flags = 'p', type = 'raster')

    grutl.import_multiple_raster_files(rasters_path, search_criteria = '*.tif')

    # Set the region to match the raster of interest
    g.region(raster = f'DEM_cropped@{mapset}')

    # We would like to mask data that falls outside the boundaries for the simulation
    r.mask(raster = f'DEM_cropped@{mapset}')

    # If you would like to check the imported files
    g.list(flags = 'p', type = 'raster')

    #* ---
    #* 5. Create time and space dataset with rain data required for the simulation.
    #* ---

    constant_rain_used = config['rain']['constant']
    if constant_rain_used:
        stds = 'constant_rain.tif'
        gcore.run_command('r.in.gdal', input = grassdata_path / stds, output = stds)
        stds = Path(stds).stem
    else:
        stds = 'rain_minutely'

        t.create(output=stds, semantictype='mean', title='Rain Rate', description='Rain rate data for itzi')

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

        start_time = grass_time['start_time']
        increment_number = int(grass_time['increment_number'])
        increment_unit = grass_time['increment_unit']

        grutl.create_rain_raster_text_file(rain_path, rain_txt_file, search_criteria='*.tif',
                                        start_time=start_time, increment_number=increment_number,
                                        increment_unit=increment_unit)

        t.register(type = 'raster', input = stds, file = str(rain_txt_file))

    #* ---
    #* 8. We now create the izti configuration file to run the simulation. Here we set all the
    #*    parameters we need and call the respective rasters and dataset we would like to use for the
    #*    simulation. All the parameters must be strings.
    #* ---

    itzi_output_path = get_path_for('itzi_output', config)
    output_itzi_file = itzi_output_path / 'itzi_config_file.ini'
    grass_input['dem'] = grass_input['dem'] or 'DEM_cropped'
    grass_input['friction'] = grass_input['friction'] or 'friction'
    # tähän ehto jos constant tehdään niin kuin tehtiin start_h.tif:n kanssa ja muuten otetaan grass tietokanta.
    grass_input['rain'] = stds

    start_y = grass_input['start_y']
    if not start_y:
        grass_input['start_y'] = ''
    inflow = grass_input['inflow']
    if not inflow:
        grass_input['inflow'] = ''
    bctype = grass_input['bctype']
    if not bctype:
        grass_input['bctype'] = ''
    bcval = grass_input['bcval']
    if not bcval:
        grass_input['bcval'] = ''
    
    infiltration = grass_input['infiltration']
    if not infiltration:

        if config['rain']['infiltration_rate']:

            inf_stds = 'infiltration_minutely'

            t.create(output=inf_stds, semantictype='mean', title='Infiltration Rate',
                    description='Infiltration rate data for itzi')

            infiltration_path = rasters_path / 'infiltration'

            grutl.import_multiple_raster_files(infiltration_path, search_criteria = '*.tif')

            infiltration_txt_file = infiltration_path / 'infiltration_registering_data.txt'
            start_time = grass_time['start_time']
            increment_number = int(grass_time['increment_number'])
            increment_unit = grass_time['increment_unit']

            grutl.create_rain_raster_text_file(infiltration_path, infiltration_txt_file,
                                            search_criteria='*.tif', start_time=start_time,
                                            increment_number=increment_number,
                                            increment_unit=increment_unit)

            t.register(type = 'raster', input = inf_stds, file = str(infiltration_txt_file))

            grass_input['infiltration'] = inf_stds

        else:
            grass_input['infiltration'] = 'infiltration_0'

    # Add the start_h file
    start_h = grass_input['start_h']
    if not start_h:
        start_h = ''
    else:
        # Vaikuttaako tässä se, että alkukartta on cropattu?
        gcore.run_command('r.in.gdal', input = grassdata_path / start_h, output = start_h)
        start_h = Path(start_h).stem
    grass_input['start_h'] = start_h

    grass_input['losses'] = grass_input['losses'] or 'losses'
    grass_output['prefix'] = grass_output['prefix'] or f'{mapset}_itzi'
    grass_statistics['stats_file'] = grass_statistics['stats_file'] or itzi_output_path / f'{mapset}_itzi.csv'

    grutl.create_itzi_config_file(output_itzi_file, config)
     
    #* ---
    #* 9. End the current grass session
    #* ---

    grutl.end_GRASS_session(user)

    #* ---
    #* X. Run the simulation in the terminal:
    #*
    #*     $ itzi run output_itzi_file
    #*
    #* XI. Export the rasters of interest to GeoTiff file so that they can be analysed separately using:
    #*     'export_tif.py'
    #* ---

    print('\nGRASS has been set correctly. Ready to run itzi\n')

if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)