# A file for refactoring and testable code
import sys
import os
import glob
import numpy as np
import rasterio as rio
import modules.water_depth_map as wdm
from pathlib import Path
import modules.Landcover as LC
import modules.raster_utils as utl
import yaml

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

def read_dem_file_array(input_filename):
    with rio.open(input_filename, 'r') as dem_data:
        return dem_data.read(1)

def read_dem_file_metadata(input_filename):
    with rio.open(input_filename, 'r') as dem_data:
        return dem_data.meta

def create_start_water_depth_file(input_file_path, output_file_path, config):
    water_depth_method_exists_in_config = (config['water']['method'] != None)
    
    #if water_depth_method_exists_in_config:
        #dem, profile = wdm.single_dem(input_file_path)
        # TÄMÄ EI TOIMI YML:lle
        #watermap = wdm.main(dem, profile)
        #wdm.write_watermap_to_output(watermap, output_file_path, profile)

def add_missing_subfolders(config_dict):
    store_folder = get_path_for('processed', config_dict)  # dump final modified files here
    print('Checking that output folders exist...', end=' ')
    ensure_folders_exist(config_dict)
    
    # subfolders for this script
    (get_path_for('temporary', config_dict) / 'rain').mkdir(exist_ok=True)
    (get_path_for('grass_db', config_dict) / 'rain').mkdir(exist_ok=True)
    (get_path_for('grass_db', config_dict) / 'infiltration').mkdir(exist_ok=True)
    print('Done.')
    return 1

def create_config_dictionary_from_yml_file(config_filename):
    with open(config_filename) as config_yml:
        parsed_yaml_file = yaml.load(config_yml, Loader=yaml.FullLoader)
    return parsed_yaml_file

def merge_if_needed(config):
    # Tämä on huono ohjelma koska tekee mergen ja palaluttaa tiedoston, eli kaksi asiaa.
    merge_cfg = config['merging']
    merge_is_used = merge_cfg['DEM_merge_boolean']
    if merge_is_used:
        rasters_folder_path = get_path_for('DEM', config)
        merged_file_path = os.path.join(get_path_for('temporary', config), 'DEM_merged.tif')
        merge_search_criteria = merge_cfg['DEM_merge_search_criteria']
        print(merge_search_criteria)
        utl.raster_merge(rasters_folder_path, merged_file_path, search_criteria=merge_search_criteria)
        raster_to_crop_path = merged_file_path
    else:
        raster_to_crop_path = config['cropping']['DEM_file_path_to_crop']
    print("tassa:")
    print(raster_to_crop_path)
    return raster_to_crop_path

def set_up_config_for_crop(config):
    config['cropping']['crop_vector_file'] = config['cropping']['shapefile_file_path']
    config['cropping']['crop_polygon_crs'] = config['cropping']['polygon_crs']
    x_coords = config['cropping']['x_polygon_coords']
    y_coords = config['cropping']['y_polygon_coords']
    
    config['cropping']['crop_polygon_coords'] = None
    if x_coords and y_coords:
        config['cropping']['crop_polygon_coords'] = list(zip(x_coords, y_coords))
    
    config['cropping']['null_value'] = config['cropping']['mask_value'] if config['cropping']['mask_value'] else -9999.

    config['cropping']['dem_cropped_file_path'] = os.path.join(get_path_for('grass_db', config), 'DEM_cropped.tif')
    
    config['cropping']['crop_search_criteria'] = config['cropping']['DEM_crop_search_criteria']
    
    return config

def crop_merged_rasters(raster_to_crop_path, config):
    #* ---
    #* Cropping utilities
    #* ---
    config['cropping']['crop_is_used'] = config['cropping']['DEM_crop_boolean']
    config = set_up_config_for_crop(config)
    if config['cropping']['crop_is_used']:
        utl.raster_crop(raster_to_crop_path,
                        cropped_file_path=config['cropping']['dem_cropped_file_path'],
                        search_criteria=config['cropping']['crop_search_criteria'],
                        vector_file=config['cropping']['crop_vector_file'],
                        polygon_coords=config['cropping']['crop_polygon_coords'],
                        polygon_crs=config['cropping']['crop_polygon_crs'],
                        mask_value=config['cropping']['null_value'],
                        mask_all_touched=True
                    )
    return config
    # If you would like to visualize it
    #rio.plot.show(rasterio.open(cropped_file_path), cmap='terrain')

def set_up_rain_rasters(config):
    #* 3.1. Here a map describing the water depth in the initial timestep is created so that
    # it corresponds to the cropped DEM above or just the dem if it is not cropped.
    output_file_path = os.path.join(get_path_for('grass_db', config), 'start_h.tif')
    input_file_for_depth_and_constant_rain = config['cropping']['dem_cropped_file_path'] if config['cropping']['crop_is_used'] else (config['merging']['merged_file_path'] if config['merging']['DEM_merge_boolean'] else config['cropping']['raster_to_crop_path'])
    create_start_water_depth_file(input_file_for_depth_and_constant_rain, output_file_path, config), 

    #* 4. Now we proceed to extract the rain rasters
    rain_cfg = config['rain']
    constant_rain_is_used = rain_cfg['constant']
    if constant_rain_is_used:
        # create constant rain file
        intensity_constant = rain_cfg['intensity_constant']
        rain_output_file_path = os.path.join(get_path_for('grass_db', config), 'constant_rain.tif')
        rain_dem, rain_profile = wdm.single_dem(input_file_for_depth_and_constant_rain)
        wdm.write_watermap_to_output(np.ones_like(rain_dem) * intensity_constant, rain_output_file_path, rain_profile)

    else:
        ascii_files_path = rain_cfg['ascii_files_path']
        GTiff_files_path = rain_cfg['GTiff_files_path']

        try:
            xmin = rain_cfg['x_min']
            ymax = rain_cfg['y_max']
            EPSG_code = rain_cfg['EPSG_code']
        except ValueError as e:
            raise ValueError('xmin, ymax and EPSG_code are mandatory') from e

        if GTiff_files_path:
            rain_crop_search_criteria = config['rain']['rain_crop_search_criteria']
        elif ascii_files_path:
            print('\nRain data is in ascii format and its being converted to Tif format...\n')

            GTiff_files_path = os.path.join(get_path_for('temporary', config), 'rain')

            utl.ascii_to_geotiff(ascii_files_path, GTiff_files_path, xmin, ymax,
                                search_criteria=rain_cfg['ascii_search_criteria'], CRS_code=EPSG_code,
                                xmax=rain_cfg['x_max'], xrotation=rain_cfg['x_rotation'], xres=rain_cfg['x_res'],
                                yrotation=rain_cfg['y_rotation'], ymin=rain_cfg['y_min'], yres=rain_cfg['y_res'])
            print('\nDONE\n')
        else:
            raise ValueError('No rain data path provided')

        if rain_cfg['rain_relocation']: # BUG: Any string is True, not just "True"
            print('\nRain data is being relocated...\n')
            utl.rain_relocation(GTiff_files_path, xmin, ymax, rain_cfg['X_target'], rain_cfg['Y_target'], rain_cfg['X_radar_rain'],
                                rain_cfg['Y_radar_rain'], search_criteria=rain_cfg['relocation_search_criteria'], CRS_code=EPSG_code,
                                xmax=rain_cfg['x_max'], xrotation=rain_cfg['x_rotation'], xres=rain_cfg['x_res'],
                                yrotation=rain_cfg['y_rotation'], ymin=rain_cfg['y_min'], yres=rain_cfg['y_res'])

            rain_crop_search_criteria = '*_relocated.tif'

        print('\nDONE\n')

        #* 5. Let's check if the rain rasters have same projection as raster of interest
        reference_file_path = os.path.join(get_path_for('grass_db', config), 'DEM_cropped.tif')

        if Path(GTiff_files_path).suffix == '':
            source_fname = next(glob.iglob(f"{GTiff_files_path}/*.tif"))  # we need only one file
            reproj = utl.raster_check_projection(source_fname, ref_fp = reference_file_path)
        else:
            reproj = utl.raster_check_projection(GTiff_files_path, ref_fp = reference_file_path)

        if reproj:
            root_reproj_rain_file = os.path.join(get_path_for('temporary', config), 'rain')

            print('\nRain file(s) being reprojected to match DEMs projection...\n')
            utl.raster_reproject(GTiff_files_path, reproj_fp=root_reproj_rain_file,
                                ref_fp=reference_file_path, search_criteria=rain_crop_search_criteria)

            rain_crop_search_criteria = "*_reproj.tif"

        if rain_cfg['rain_crop_boolean']:
            print('\nRain files are being cropped...\n')
            rain_cropped_file_path = os.path.join(get_path_for('grass_db', config), 'rain')
            utl.raster_crop(GTiff_files_path, cropped_file_path=rain_cropped_file_path,
                            search_criteria=rain_crop_search_criteria, vector_file=config['crop_vector_file'],
                            polygon_coords=config['crop_polygon_coords'], polygon_crs=config['crop_polygon_crs'],
                            mask_value=config['null_value'], mask_all_touched=True)

            print('\nDONE\n')

def set_up_landcover(config):
    merge_cfg = config['merging']
    if merge_cfg['landcover_merge_boolean']:
        landcover_folder_path = merge_cfg['landcover_file_path_to_merge']
        merged_file_path = os.path.join(get_path_for('temporary', config), 'landcover_merged.tif')
        merge_search_criteria  = merge_cfg['landcover_merge_search_criteria']
        utl.raster_merge(landcover_folder_path, merged_file_path, search_criteria=merge_search_criteria)
        landcover_to_crop_path = merged_file_path
    else:
        landcover_to_crop_path = config['cropping']['landcover_file_path_to_crop']

    # Let's check that the projection of the landcover map matches the projection of DEM map
    root_reproj_landcover_file = os.path.join(get_path_for('temporary', config), 'landcover_reproj.tif')
    reproj = utl.raster_check_projection(landcover_to_crop_path, ref_fp = config['dem_cropped_file_path'])

    if reproj:
        print('\nReprojection needed for Landcover file to match DEM projection\n')
        utl.raster_reproject(landcover_to_crop_path, reproj_fp=root_reproj_landcover_file,
                            ref_fp=config['dem_cropped_file_path'])
        landcover_to_crop_path = root_reproj_landcover_file
    config['cropping']['landcover_to_crop_path'] = landcover_to_crop_path

def crop_landcover(config):
    if config['cropping']['landcover_crop_boolean']:
        config['landcover_cropped_file_path'] = os.path.join(get_path_for('grass_db', config), 'landcover_cropped.tif')
        config['crop_search_criteria'] = config['cropping']['landcover_crop_search_criteria']

        utl.raster_crop(config['cropping']['landcover_to_crop_path'], cropped_file_path=config['landcover_cropped_file_path'],
                        search_criteria=config['crop_search_criteria'], vector_file=config['crop_vector_file'],
                        polygon_coords=config['crop_polygon_coords'], polygon_crs=config['crop_polygon_crs'],
                        mask_value=config['null_value'], mask_all_touched=True)

def initialize_landcover(config):
    landcover = LC.Landcover(config['landcover_cropped_file_path'])
    landcover.get_friction()
    landcover.get_losses()
    rain_cfg = config['rain']
    merge_cfg = config['merging']
    if config['imperviousness']['imperviousness_file_boolean']:
        if merge_cfg['imperviousness_merge_boolean']:
            merged_file_path = os.path.join(get_path_for('temporary', config), 'imperviousness_merged.tif')
            utl.raster_merge(merge_cfg['imperviousness_file_path_to_merge'], merged_file_path,
                            search_criteria=merge_cfg['imperviousness_merge_search_criteria'])
            imperviousness_file_path = merged_file_path
        else:
            imperviousness_file_path = config['cropping']['imperviousness_file_path_to_crop']

        root_reproj_imperviousness_file = os.path.join(get_path_for('temporary', config), 'imperviousness_reproj.tif')
        reproj = utl.raster_check_projection(imperviousness_file_path, ref_fp=config['dem_cropped_file_path'])

        if reproj:

            print('\nImperviousness file is being reprojected to match DEMs projection...\n')
            utl.raster_reproject(imperviousness_file_path, reproj_fp=root_reproj_imperviousness_file,
                                ref_fp=config['dem_cropped_file_path'])
            imperviousness_file_path = root_reproj_imperviousness_file

        if config['cropping']['imperviousness_crop_boolean']:

            imperviousness_cropped_file_path = os.path.join(get_path_for('grass_db', config),
                                                            'imperviousness_cropped.tif')
            crop_search_criteria = config['cropping']['imperviousness_crop_search_criteria']

            utl.raster_crop(imperviousness_file_path,cropped_file_path=imperviousness_cropped_file_path,
                            search_criteria=crop_search_criteria, vector_file=config['crop_vector_file'],
                            polygon_coords=config['crop_polygon_coords'], polygon_crs=config['crop_polygon_crs'],
                            mask_value=config['null_value'], mask_all_touched=True)

            imperviousness_file_path = imperviousness_cropped_file_path

            utl.set_raster_resolution(imperviousness_file_path, config['landcover_cropped_file_path'],
                                    binary=False, mask_value=241)

        if rain_cfg['infiltration_rate']:
            infiltration_output_path = os.path.join(get_path_for('grass_db', config), 'infiltration')

        else:
            infiltration_output_path = get_path_for('grass_db', config)

        landcover.get_infiltration(imperviousness_file_path, os.path.join(get_path_for('grass_db', config), 'rain'),
                                output_folder=infiltration_output_path,
                                infiltration_rate=rain_cfg['infiltration_rate'])


def main(config_file_name):
    #* ---
    #* We read the input information from the ini file and then the magic happens
    #* ---
    
    # instantiate
    config = create_config_dictionary_from_yml_file(config_file_name)

    # Path setup
    add_missing_subfolders(config)

    #* 1. Downdload DEMs of interest:
    #*     - e.g. 2x2 DEM files from https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en
    #* 2. If the region of interest lays in the intersection of different DEMs, we can merge them into a
    #*    single GeoTiff file:
    config['cropping']['raster_to_crop_path'] = merge_if_needed(config)
    
    #* 3. We would probably like to do the analysis over a specific region and not the entire map.
    #*    Therefore, we proceed to crop the file according to a given polygon cooredinates or vector
    #*    shapefile.
    #*
    #*    If vector file not given, then both polygon coordinates and crs MUST be given
    config = crop_merged_rasters(config['cropping']['raster_to_crop_path'], config)

    set_up_rain_rasters(config)

    #* 6. We now turn to landuse information to define: friction, losses and infiltration.
    #*    Landuse information for Finland downloaded from:
    #*    https://www.avoindata.fi/data/fi/dataset/corine-maanpeite-2018
    #*    The file provides the different classes of the terrain with 20x20 m resolution
    #*
    #*    Let's first check whether if the landcover map has the same projection as our DEM raster
    #*    *As I already ran it once, reprojection was needed and I will skip it as it takes some minutes
    #*     to reproject.
    #*    But if its for the first time, do it!
    set_up_landcover(config)
    
    #* 7. Let's crop the landcover file before extracting the information we need
    crop_landcover(config)

    #* 8. Now we intialize the Landcover class and extract friction, losses and infiltration
    initialize_landcover(config)
    print('\nDONE with DEM, LANDCOVER and RAIN rasters\n')


if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)