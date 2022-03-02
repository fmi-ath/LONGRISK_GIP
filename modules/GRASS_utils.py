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
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import temporal as t

def initiate_GRASS_session(mygisdb, mylocation, mymapset, CRS_code = 4326):

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

    crs_as_epsg = f'EPSG:{CRS_code}'

    os.environ.update(dict(GRASS_COMPRESS_NULLS = '1', GRASS_COMPRESSOR = 'ZSTD'))

    # create a PERMANENT mapset; create a Session instance
    PERMANENT = Session()
    # hint: EPSG code lookup: https://epsg.io
    PERMANENT.open(gisdb = mygisdb, location = mylocation,
                   create_opts = crs_as_epsg)

    if mymapset == 'PERMANENT':

        return PERMANENT

    PERMANENT.close()

    user = Session()
    user.open(gisdb = mygisdb, location = mylocation, mapset = mymapset,
                   create_opts = crs_as_epsg)

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
        #g.remove(type = 'all', pattern = '*', flags = 'f')

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


def create_rain_raster_text_file(rain_raster_path, output_file, search_criteria = '*.tif',
                                 start_time = None, increment_number = None, increment_unit = None):

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

    files_path = os.path.join(rain_raster_path, search_criteria)

    file_list = glob.glob(files_path)

    n_files = len(file_list)

    assert n_files != 0, 'There are no files in the path given'

    f = open(output_file, "w", encoding='utf-8')

    if start_time is None:

        for filename in file_list:

            _, tail = os.path.split(filename)

            #print(tail)
            new_raster = os.path.splitext(tail)[0]

            #print('hour: %d, minute: %d' % (hour, minute))
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

            #print(tail)
            new_raster = os.path.splitext(tail)[0]

            #print('hour: %d, minute: %d' % (hour, minute))
            f.write(f'{new_raster}|{initial_time:%Y-%m-%d %H:%M}|{current_time:%Y-%m-%d %H:%M}\n')

            initial_time = current_time
            current_time += delta

    f.close()


def create_itzi_config_file(output_file, config_data: dict): 


    """
    # grass
    grass_bin=None, grassdata=None, location=None, mapset=None,
    region=None, mask=None"""
    """Creates the .ini file needed for running itzi simulation according to parameters specified
    by user.

    For this function, all parameters MUST be given as strings.
    See itzi documentation for info about the parameters: itzi.readthedocs.io/

    Parameters
    ----------
    output_file : str
        Name of the output .ini file.

    Returns
    -------
    out : file
        Text file with parameters for itzi simnulation
    """
    grass_info = config_data['grass_info']
    grass_time = config_data['grass_time']
    grass_input = config_data['grass_input']
    grass_output = config_data['grass_output']
    grass_options = config_data['grass_options']
    drainage_kws = config_data['drainage']
    grass_statistics = config_data['grass_statistics']

    if (grass_time['record_step'] is None) or (grass_input['dem'] is None) or (grass_input['friction'] is None):
        raise ValueError('record_step, dem and friction are mandatory arguments for the simulation')
    
    # These values are none as they lacked a value in the former version.
    effective_porosity=None
    capillary_pressure=None
    hydraulic_conductivity=None 
    
    input_kws = {
        'dem': grass_input['dem'],
        'friction': grass_input['friction'],
        'start_h': grass_input['start_h'],
        'start_y': grass_input['start_y'],
        'rain': grass_input['rain'],
        'inflow': grass_input['inflow'],
        'bctype': grass_input['bctype'],
        'bcval': grass_input['bcval'],
        'infiltration': grass_input['infiltration'],
        'effective_porosity': effective_porosity,
        'capillary_pressure': capillary_pressure,
        'hydraulic_conductivity': hydraulic_conductivity,
        'losses': grass_input['losses'],
    }    
    region=None
    mask=None
    grass_kws = {
        "grass_bin": grass_info['grass_bin_path'],
        "grassdata": grass_info['grassdata_path'],
        "location": grass_info['location'],
        "mapset": grass_info['mapset'],
        "region": region,
        "mask": mask,
    }

    config = configparser.ConfigParser()

    def _insert_in_loop(section_name, items):
        config.add_section(section_name)
        for key, value in items.items():
            if value is not None:
                config.set(section_name, key, str(value))

    #? TIME SECTION
    secname = 'time'
    config.add_section(secname)

    # start_time -- Start of the simulation. Format yyyy-mm-dd HH:MM
    # end_time -- End of the simulation. Format yyyy-mm-dd HH:MM
    # duration -- Duration of the simulation. Format HH:MM:SS
    # record_step -- Duration between two records. Format HH:MM:SS
    if (grass_time['start_time'] is not None) and (grass_time['end_time'] is not None):
        config.set(secname, 'start_time', str(grass_time['start_time']))
        config.set(secname, 'end_time', str(grass_time['end_time']))
    elif (grass_time['start_time'] is not None) and (grass_time['duration'] is not None):
        config.set(secname, 'start_time', str(grass_time['start_time']))
        config.set(secname, 'duration', str(grass_time['duration']))
    elif grass_time['duration'] is not None:
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
    config.set(secname, 'prefix', str(grass_output['prefix']))
    values = grass_output['values'] if grass_output['values'] is not None else []
    config.set(secname, 'values', ', '.join(values))

    #? STATISTICS SECTION
    secname = 'statistics'
    config.add_section(secname)
    config.set(secname, 'stats_file', str(grass_statistics['stats_file']))

    #? OPTIONS SECTION
    _insert_in_loop('options', grass_options)

    #? DRAINAGE SECTION
    _insert_in_loop('drainage', drainage_kws)

    #? GRASS SECTION
    _insert_in_loop('grass', grass_kws)

    # Write config to file
    with open(output_file, "w", encoding='utf-8') as f:
        config.write(f, space_around_delimiters=True)

#def GRASS_export_rasters(raster_path, output_path, search_criteria = '*'):
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
