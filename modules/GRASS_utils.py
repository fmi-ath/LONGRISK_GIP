# pylint: disable=import-error  # GRASS is installed inside container, not in athras -> pylint doesn't find it
import configparser
import os
import glob
from pathlib import Path
from datetime import datetime, timedelta
from re import sub

import numpy as np
from grass_session import Session
from grass.script import core as gcore
import grass.script as gscript
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import temporal as t

def initiate_GRASS_sesion(mygisdb, mylocation, mymapset, CRS_code = 4326):

    """Initiates a GRASS session with the given parameters.

        Parameters
        ----------
        mygisdb : str
            Path to the GIS database.
        mylocation : str
            Name of the location for the current session.
        mymapset : str
            Name of the mapset for the current session.
        CRS_code : int
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

def end_GRASS_sesion(session_name):

    """Ends the current GRASS session.

        Parameters
        ----------
        mymapset : str
            Name of the mapset for the current session.
    """

    session_name.close()

    return

def check_remove_existing_files(mapset, search_criteria = '*', maps_remove = True, datasets_remove = True):

    """
    Based on:
    https://grass.osgeo.org/grass79/manuals/libpython/temporal_framework.html

        Check if there are already files in the current session and removes them.

        Parameters
        ----------
        mapset : str
            Name of the mapset for the current session.
        search_criteria : list or str
            List of criteria for searhing the maps or datasets to be removed.
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

    f = open(output_file, "w", encoding='utf-8')

    files_path = os.path.join(rain_raster_path, search_criteria)

    file_list = glob.glob(files_path)

    n_files = len(file_list)

    assert n_files != 0, 'There are no files in the path given'

    if start_time is None:

        for file in glob.iglob(files_path):

            _, tail = os.path.split(file)

            #print(tail)
            new_raster = os.path.splitext(tail)[0]

            #print('hour: %d, minute: %d' % (hour, minute))
            f.write(new_raster + '\n')

    else:

        if not isinstance(start_time, datetime):

            start_time = datetime.strptime(start_time, "%Y-%m-%d %H:%M")
            print('The start_time was converted to datetime: ', start_time.strftime("%Y-%m-%d %H:%M"))

        delta = timedelta(**{increment_unit: increment_number})

        initial_time = start_time - delta
        current_time = start_time

        for file in sorted(glob.iglob(files_path)):

            _, tail = os.path.split(file)

            #print(tail)
            new_raster = os.path.splitext(tail)[0]

            #print('hour: %d, minute: %d' % (hour, minute))
            f.write(f'{new_raster}|{initial_time:%Y-%m-%d %H:%M}|{current_time:%Y-%m-%d %H:%M}\n')

            initial_time = current_time
            current_time += delta

    f.close()

    return

def create_itzi_config_file(output_file, record_step=None, dem=None, friction=None,
                            # input
                            start_h=None, start_y=None, rain=None,
                            inflow=None, bctype=None, bcval=None, infiltration=None,
                            effective_porosity=None, capillary_pressure=None,
                            hydraulic_conductivity=None, losses=None, #dem, friction
                            # output
                            prefix=None, values=None,
                            # time
                            start_time=None, end_time=None, duration=None, #record_step
                            # statistics
                            stats_file=None,
                            # options
                            hmin=None, slmax=None, cfl=None, theta=None,
                            vrouting=None, dtmax=None, dtinf=None,
                            # drainage
                            swmm_inp=None, output=None,
                            orifice_coeff=None, free_weir_coeff=None, submerged_weir_coeff=None,
                            # grass
                            grass_bin=None, grassdata=None, location=None, mapset=None,
                            region=None, mask=None,
                            ):

    """Creates the .ini file needed for running itzi simulation according to parameters specified by user.
        For this function, all parameters MUST be given as strings.
        See itzi documentation for info about the parameters: itzi.readthedocs.io/

        Parameters
        ----------
        output_file : str
            Name of the output .ini file.
        record_step :
        dem :
        friction :

        Returns
        -------
        out : file
            Text file with parameters for itzi simnulation
    """

    if (record_step is None) or (dem is None) or (friction is None):
        raise ValueError('record_step, dem and friction are mandatory arguments for the simulation')

    input_kws = {
        'dem': dem,
        'friction': friction,
        'start_h': start_h,
        'start_y': start_y,
        'rain': rain,
        'inflow': inflow,
        'bctype': bctype,
        'bcval': bcval,
        'infiltration': infiltration,
        'effective_porosity': effective_porosity,
        'capillary_pressure': capillary_pressure,
        'hydraulic_conductivity': hydraulic_conductivity,
        'losses': losses,
    }

    options_kws = {
        'hmin': hmin,
        'slmax': slmax,
        'cfl': cfl,
        'theta': theta,
        'vrouting': vrouting,
        'dtmax': dtmax,
        'dtinf': dtinf,
    }

    drainage_kws = {
        "swmm_inp": swmm_inp,
        "output": output,
        "orifice_coeff": orifice_coeff,
        "free_weir_coeff": free_weir_coeff,
        "submerged_weir_coeff": submerged_weir_coeff,
    }

    grass_kws = {
        "grass_bin": grass_bin,
        "grassdata": grassdata,
        "location": location,
        "mapset": mapset,
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
    if (start_time is not None) and (end_time is not None):
        config.set(secname, 'start_time', str(start_time))
        config.set(secname, 'end_time', str(end_time))
    elif (start_time is not None) and (duration is not None):
        config.set(secname, 'start_time', str(start_time))
        config.set(secname, 'duration', str(duration))
    elif duration is not None:
        config.set('time', 'duration', str(duration))
    else:
        time_error_message = ('Possible combinations for time:\n'
                              '- start_time and end_time\n'
                              '- start_time and duration\n'
                              '- duration only\n')
        raise ValueError(time_error_message)

    config.set(secname, 'record_step', str(record_step))

    #? INPUT SECTION
    _insert_in_loop('input', input_kws)

    #? OUTPUT SECTION
    secname = 'output'
    config.add_section(secname)
    config.set(secname, 'prefix', str(prefix))
    values = values if values is not None else []
    config.set(secname, 'values', ', '.join(values))

    #? STATISTICS SECTION
    secname = 'statistics'
    config.add_section(secname)
    config.set(secname, 'stats_file', str(stats_file))

    #? OPTIONS SECTION
    _insert_in_loop('options', options_kws)

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
