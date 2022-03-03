import os, subprocess, sys
import modules.GRASS_utils as grutl
import modules.common as cm
from pathlib import Path
from grass.script import core as gcore
from grass.pygrass.modules.shortcuts import general as g, raster as r, temporal as t


def set_up_start_h_file_if_specified(config):
    start_h = config['grass_input']['start_h']
    if not start_h:
        return ''
    gcore.run_command('r.in.gdal', input = Path(config['grass_info']['grass_db']) / start_h, output = start_h)
    return Path(start_h).stem

def create_gisrc_file_for_grass_if_needed(grass_info):
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

def set_up_mapset_path(grass_info):
    grassdata_path = Path(grass_info['grass_db'])
    mapset_path = Path(os.path.join(grassdata_path, grass_info['location'], grass_info['mapset']))
    if mapset_path.exists():
        subprocess.call(["rm", "-r", str(mapset_path)])

def add_all_relevant_rasters_to_grass(config):
    mapset = config['grass_info']['mapset'] 
    g.list(flags = 'p', type = 'raster')
    grutl.import_multiple_raster_files(Path(config['grass_info']['grass_db']))
    # Set the region to match the raster of interest
    g.region(raster = f'DEM_cropped@{mapset}')
    # We would like to mask data that falls outside the boundaries for the simulation
    r.mask(raster = f'DEM_cropped@{mapset}')
    # If you would like to check the imported files
    g.list(flags = 'p', type = 'raster')
    config['grass_input']['start_h'] = set_up_start_h_file_if_specified(config)

def set_up_infiltration(config):
    if not config['grass_input']['infiltration']:
        config['grass_input']['infiltration'] = 'infiltration_0'
        if config['rain']['infiltration_rate']:
            inf_stds = 'infiltration_minutely'
            t.create(output=inf_stds, semantictype='mean', title='Infiltration Rate', description='Infiltration rate data for itzi')
            infiltration_path = Path(config['grass_info']['grass_db']) / 'infiltration'
            grutl.import_multiple_raster_files(infiltration_path)
            infiltration_txt_file = infiltration_path / 'infiltration_registering_data.txt'
            grutl.create_rain_raster_text_file(infiltration_path, infiltration_txt_file, config['grass_time'])
            t.register(type = 'raster', input = inf_stds, file = str(infiltration_txt_file))
            config['grass_input']['infiltration'] = inf_stds

def set_up_rain(config):
    if config['rain']['constant']:
        stds = 'constant_rain.tif'
        gcore.run_command('r.in.gdal', input = Path(config['grass_info']['grass_db']) / stds, output = stds)
        return Path(stds).stem
    
    stds = 'rain_minutely'
    t.create(output=stds, semantictype='mean', title='Rain Rate', description='Rain rate data for itzi')

    #* 4.1. Import the minutely recorded radar rain events.
    rain_path = Path(config['grass_info']['grass_db']) / 'rain'
    grutl.import_multiple_raster_files(rain_path)
    g.list(flags = 'p', type = 'raster')

    #* 4.2. Register the minutely rain radar data in the time and space dataset created previously.
    #*    As there are many files and the only way grass allows registering many files at once is with a
    #*    .txt file indicating the name of the files, we will first create the file and then register.
    rain_txt_file = rain_path / 'rain_registering_data.txt'
    grutl.create_rain_raster_text_file(rain_path, rain_txt_file, config['grass_time'])
    t.register(type = 'raster', input = stds, file = str(rain_txt_file))
    return stds

def main(config_filename):
    config = cm.create_config_dictionary_from_config_file(config_filename)
    create_gisrc_file_for_grass_if_needed(config['grass_info'])

    #* 1. We define the paths in which we would like to work out the simulation
    #*     - grass_db: str. Path of working location
    #*     - location: str. Name of the location (Ex. Helsinki)
    #*     - mapset: str. User that will be working on that location
    #* 2. We initiate the grass session with 'initiate_GRASS_session'. If sessions does
    #*    not exists, it creates it. If it exists, it opens it and load the files
    set_up_mapset_path(config['grass_info'])
    user = grutl.initiate_GRASS_session(config['grass_info'])

    #* 3. Add all the relevant rasters for the simulation.
    #*    For conviniency, the name of the rasters in GRASS will be the same as the
    #*    .tif file names but without the extension (i.e. input = Helsinki_cropped.tif,
    #*    output = Helsinki_cropped)
    add_all_relevant_rasters_to_grass(config)

    #* 4. Create time and space dataset with rain data required for the simulation.
    config['grass_input']['rain'] = set_up_rain(config)

    #* 5. Here infiltration is set up for grass
    set_up_infiltration(config)
    
    #* 6. We now create the itzi configuration file to run the simulation. Here we set all the
    #*    parameters we need and call the respective rasters and dataset we would like to use for the
    #*    simulation.
    #*    Then the GRASS session is finished
    grutl.create_itzi_config_file(config)
    grutl.end_GRASS_session(user)
    print('\nGRASS has been set correctly. Ready to run itzi\n')


if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters_helsinki.yml'")) from e
    main(config_file_name)