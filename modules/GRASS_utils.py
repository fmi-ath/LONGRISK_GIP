# pylint: disable=import-error  # GRASS is installed inside container, not in athras -> pylint doesn't find it
import configparser
import os
import glob
from pathlib import Path
from datetime import datetime, timedelta

from grass_session import Session
from grass.script import core as gcore
import grass.script as gscript
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import temporal as t

def initiate_GRASS_session(grass_info):
    """Initiates a GRASS session with the given parameters.

        Parameters
        ----------
        mygisdb : str
            Path to the GIS database.
        mylocation : str
            Name of the location for the current session.
        mymapset : str
            Name of the mapset for the current session.
        CRS_code : int | str
            Optional CRS for the current session.
    """
    mygisdb = grass_info['grass_db']
    mylocation = grass_info['location']
    mymapset = grass_info['mapset']
    CRS_code = grass_info['CRS'] or 4326

    crs_as_epsg = f'EPSG:{CRS_code}'

    os.environ.update(dict(GRASS_COMPRESS_NULLS = '1', GRASS_COMPRESSOR = 'ZSTD'))

    # create a PERMANENT mapset; create a Session instance
    PERMANENT = Session()
    # hint: EPSG code lookup: https://epsg.io
    PERMANENT.open(gisdb = mygisdb, location = mylocation, create_opts = crs_as_epsg)

    if mymapset == 'PERMANENT':

        return PERMANENT

    PERMANENT.close()

    user = Session()
    user.open(gisdb = mygisdb, location = mylocation, mapset = mymapset, create_opts = crs_as_epsg)

    return user

def end_GRASS_session(session_name):
    """Ends the current GRASS session."""
    session_name.close()

def check_remove_existing_files(maps_remove = True, datasets_remove = True):
    """
    Based on:
    https://grass.osgeo.org/grass79/manuals/libpython/temporal_framework.html

        Check if there are already files in the current session and removes them.

        Parameters
        ----------
        maps_remove: bool
            Boolean indicating if maps (if any) are to be removed
        datasets_remove: bool
            Boolean indicating if datasets (if any) are to be removed
    """
    print('\n Files to be deleted: \n\n')

    g.list(flags = 'p', type = 'all')
    t.list()

    if maps_remove:
        g.remove(type = 'all', pattern = '*', flags = 'fb')

    if datasets_remove:
        t.remove(flags = 'rf')

def import_multiple_raster_files(path, search_criteria = '*.tif'):
    """Imports many rasters to the current GRASS session for a given search criteria.

        Parameters
        ----------
        path : str or Path
            Folder where to search for the files
        search_criteria : list or str
            List of criteria for searhing the rasters to be imported.
    """
    files_path = os.path.join(path, search_criteria)

    for name in glob.iglob(files_path):
        gcore.run_command('r.in.gdal', input = name, output = Path(name).stem)

def create_rain_raster_text_file(rain_raster_path, output_file, grass_time, search_criteria = '*.tif'):
    """Creates the text file needed for registering the rain rasters in the created space and time
    dataset for the simulation

        Parameters
        ----------
        rain_raster_path : str
            Path to the rain rasters.
        output_file : str
            Name of the rain rasters names file.
        search_criteria : list or str
            List of criteria for searhing the rasters to be imported.
        start_time: str or datetime object
            Both with format 'yyyy-mm-dd HH:MM' or just 'yyyy-mm-dd'
        increment_number: int
            Time increment of the rain files
        increment_unit: str
            Specifies if increment_number is seconds, minutes, hours, days, months or years

        Returns
        -------
        out : file
            Text file with rain raster names and time intervals if indicated
    """
    start_time = grass_time['start_time']
    increment_number = grass_time['increment_number']
    increment_unit = grass_time['increment_unit']

    files_path = os.path.join(rain_raster_path, search_criteria)

    file_list = glob.glob(files_path)

    n_files = len(file_list)

    assert n_files != 0, 'There are no files in the path given'

    f = open(output_file, "w", encoding='utf-8')

    if start_time is None:
        for filename in file_list:
            _, tail = os.path.split(filename)

            new_raster = os.path.splitext(tail)[0]

            f.write(new_raster)
            f.write('\n')

    else:
        if not isinstance(start_time, datetime):
            start_time = datetime.strptime(start_time, "%Y-%m-%d %H:%M")
            print(f'The start_time was converted to datetime: {start_time:%Y-%m-%d %H:%M}')

        delta = timedelta(**{increment_unit: increment_number})

        initial_time = start_time - delta
        current_time = start_time

        for filename in sorted(file_list):
            _, tail = os.path.split(filename)

            new_raster = os.path.splitext(tail)[0]

            f.write(f'{new_raster}|{initial_time:%Y-%m-%d %H:%M}|{current_time:%Y-%m-%d %H:%M}\n')

            initial_time = current_time
            current_time += delta

    f.close()

def create_itzi_config_file(config_data: dict):
    """Creates the .ini file needed for running itzi simulation according to parameters specified
    by user.

    For this function, all parameters MUST be given as strings.
    See itzi documentation for info about the parameters: itzi.readthedocs.io/

    Parameters
    ----------
    config_data : dict
        config data from the config file. The data might have some modification effects given by set_GRASS.
    Returns
    -------
    out : file
        Text file with parameters for itzi simnulation
    """
    itzi_output_path = Path(config_data['folders']['itzi_output_files'])
    output_file = itzi_output_path / 'itzi_config_file.ini'
    
    grass_input = config_data['grass_input']
    input_kws = {
        'dem': grass_input['dem'] or 'DEM_cropped',
        'friction': grass_input['friction'] or 'friction',
        'start_h': grass_input['start_h'],
        'start_y': grass_input['start_y'] or '',
        'rain': grass_input['rain'],
        'inflow': grass_input['inflow'] or '',
        'bctype': grass_input['bctype'] or '',
        'bcval': grass_input['bcval'] or '',
        'infiltration': grass_input['infiltration'],
        'effective_porosity': None,
        'capillary_pressure': None,
        'hydraulic_conductivity': None,
        'losses': grass_input['losses'] or 'losses',
    }

    mapset = config_data['grass_info']['mapset']
    grass_kws = {
        "grass_bin": config_data['grass_info']['grass_bin_path'],
        "grassdata": config_data['grass_info']['grass_db'],
        "location": config_data['grass_info']['location'],
        "mapset": mapset,
        "region": None,
        "mask": None,
    }

    config = configparser.ConfigParser()

    def _insert_in_loop(section_name, items):
        config.add_section(section_name)
        for key, value in items.items():
            if value:
                config.set(section_name, key, str(value))

    #? TIME SECTION
    secname = 'time'
    config.add_section(secname)

    grass_time = config_data['grass_time']
    if grass_time['start_time'] and grass_time['end_time']:
        config.set(secname, 'start_time', str(grass_time['start_time']))
        config.set(secname, 'end_time', str(grass_time['end_time']))
    elif grass_time['start_time'] and grass_time['duration']:
        config.set(secname, 'start_time', str(grass_time['start_time']))
        config.set(secname, 'duration', str(grass_time['duration']))
    elif grass_time['duration']:
        config.set('time', 'duration', str(grass_time['duration']))
    else:
        time_error_message = ('Possible combinations for time:\n'
                              '- start_time and end_time\n'
                              '- start_time and duration\n'
                              '- duration only\n')
        raise ValueError(time_error_message)

    config.set(secname, 'record_step', str(grass_time['record_step']))

    #? INPUT SECTION
    _insert_in_loop('input', input_kws)

    #? OUTPUT SECTION
    secname = 'output'
    config.add_section(secname)
    config.set(secname, 'prefix', str(config_data['grass_output']['prefix'] or f'{mapset}_itzi'))
    values = config_data['grass_output']['values'] or []
    config.set(secname, 'values', ', '.join(values))

    #? STATISTICS SECTION
    secname = 'statistics'
    config.add_section(secname)
    config.set(secname, 'stats_file', str(config_data['grass_statistics']['stats_file'] or itzi_output_path / f'{mapset}_itzi.csv'))

    #? OPTIONS SECTION
    _insert_in_loop('options', config_data['grass_options'])

    #? DRAINAGE SECTION
    _insert_in_loop('drainage', config_data['drainage'])

    #? GRASS SECTION
    _insert_in_loop('grass', grass_kws)

    # Write config to file
    with open(output_file, "w", encoding='utf-8') as f:
        config.write(f, space_around_delimiters=True)

def GRASS_export_rasters(output_path, mapset, search_criteria = '*'):
    """
    Based on:
    https://baharmon.github.io/python-in-grass

        Exports the itzi simulation rasters to .tiff files

        Parameters
        ----------
        output_path : str
            Path where the rasters will be located.
        mapset : str
            Name of the mapset for the current session.
        search_criteria : list or str
            List of criteria for searhing the rasters to be exported.

        Returns
        -------
        out : file
            Tiff files from itzi simulation
    """
    raster_list = gscript.list_grouped('rast', pattern = search_criteria)[mapset]

    for raster in raster_list:
        new_path = os.path.join(output_path, f'{raster}.tif')
        gcore.run_command('r.out.gdal', input=raster, output=new_path, format='GTiff')