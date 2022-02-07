# A file for refactoring and testable code
import sys
import os
import glob
import numpy as np
import rasterio as rio
import modules.water_depth_map as wdm
from pathlib import Path
from modules import common
import modules.Landcover as LC
import modules.raster_utils as utl
from configparser import ConfigParser


def read_dem_file_array(input_filename):
    with rio.open(input_filename, 'r') as dem_data:
        return dem_data.read(1)

def read_dem_file_metadata(input_filename):
    with rio.open(input_filename, 'r') as dem_data:
        return dem_data.meta

def create_start_water_depth_file(input_file_path, output_file_path, config):
    
    water_depth_method_exists_in_config = ((config['water_method']) != '')
    
    if water_depth_method_exists_in_config:
        dem, profile = wdm.single_dem(input_file_path)
        watermap = wdm.main(dem, profile)
        wdm.write_watermap_to_output(watermap, output_file_path, profile)

def add_missing_subfolders(config_dict):
    store_folder = common.get_path_for('processed', config_dict)  # dump final modified files here
    print('Checking that output folders exist...', end=' ')
    common.ensure_folders_exist(config_dict)
    
    # subfolders for this script
    (common.get_path_for('temporary', config_dict) / 'rain').mkdir(exist_ok=True)
    (common.get_path_for('grass_db', config_dict) / 'rain').mkdir(exist_ok=True)
    (common.get_path_for('grass_db', config_dict) / 'infiltration').mkdir(exist_ok=True)
    print('Done.')
    return 1

def create_config_dictionary_from_config_file(config_filename):
    config_dict = ConfigParser()
    config_dict.read(config_filename)
    return config_dict


def main(config_file_name):

    
    #* ---
    #* We read the input information from the ini file and then the magic happens
    #* ---
    

    # instantiate
    config = create_config_dictionary_from_config_file(config_file_name)

    # Path setup
    add_missing_subfolders(config)

    #* ---
    #* 1. Downdload DEMs of interest:
    #*     - e.g. 2x2 DEM files from https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en
    #* 2. If the region of interest lays in the intersection of different DEMs, we can merge them into a
    #*    single GeoTiff file:
    #* ---
    merge_is_used = config.getboolean('merging', 'DEM_merge_boolean')
    if merge_is_used:

        rasters_folder_path = common.get_path_for('DEM', config)
        merged_file_path = os.path.join(common.get_path_for('temporary', config), 'DEM_merged.tif')

        dsc = config.get('merging', 'DEM_merge_search_criteria')
        merge_search_criteria = dsc.split()

        utl.raster_merge(rasters_folder_path, merged_file_path, search_criteria=merge_search_criteria)

        raster_to_crop_path = merged_file_path

    else:

        raster_to_crop_path = config.get('cropping', 'DEM_file_path_to_crop')

    # If you would like to visualize it
    #rio.plot.show(rasterio.open(merged_file_path), cmap='terrain')

    #* ---
    #* 3. We would probably like to do the analysis over a specific region and not the entire map.
    #*    Therefore, we proceed to crop the file according to a given polygon cooredinates or vector
    #*    shapefile.
    #*
    #*    If vector file not given, then both polygon coordinates and crs MUST be given
    #* ---
    #* ---
    #* Cropping utilities
    #* ---

    crop_vector_file = config.get('cropping', 'shapefile_file_path')
    crop_polygon_crs = config.get('cropping', 'polygon_crs')

    if not crop_vector_file:
        crop_vector_file = None

    if not crop_polygon_crs:
        crop_polygon_crs = None

    x_coords = config.get('cropping', 'x_polygon_coords', fallback='')
    y_coords = config.get('cropping', 'y_polygon_coords', fallback='')
    if x_coords and y_coords:

        try:
            lx = list(map(float, x_coords.split()))
        except ValueError as e:
            raise ValueError('Polygon coords should be float or integer') from e
        try:
            ly = list(map(float, y_coords.split()))
        except ValueError as e:
            raise ValueError('Polygon coords should be float or integer') from e

        crop_polygon_coords = list(zip(lx, ly))

    else:
        crop_polygon_coords = None

    null_value = config.get('cropping', 'mask_value')
    null_value = float(null_value) if null_value else -9999.

    
    
    crop_is_used = config.getboolean('cropping', 'DEM_crop_boolean')
    if crop_is_used:

        dem_cropped_file_path = os.path.join(common.get_path_for('grass_db', config), 'DEM_cropped.tif')
        crop_search_criteria = config.get('cropping', 'DEM_crop_search_criteria')

        utl.raster_crop(raster_to_crop_path,
                        cropped_file_path=dem_cropped_file_path,
                        search_criteria=crop_search_criteria,
                        vector_file=crop_vector_file,
                        polygon_coords=crop_polygon_coords,
                        polygon_crs=crop_polygon_crs,
                        mask_value=null_value,
                        mask_all_touched=True
                    )

    # If you would like to visualize it
    #rio.plot.show(rasterio.open(cropped_file_path), cmap='terrain')

    #* ---
    #* 3.1. Here a map describing the water depth in the initial timestep is created so that 
    # it corresponds to the cropped DEM above or just the dem if it is not cropped.
    #* ---


    output_file_path = os.path.join(common.get_path_for('grass_db', config), 'start_h.tif')
    input_file_for_depth_and_constant_rain = dem_cropped_file_path if crop_is_used else (merged_file_path if merge_is_used else raster_to_crop_path)
    create_start_water_depth_file(input_file_for_depth_and_constant_rain, output_file_path), 

    #* ---
    #* 4. Now we proceed to extract the rain rasters
    #* ---

    constant_rain_is_used = config.getboolean('rain', 'constant')
    if constant_rain_is_used:
        # create constant rain file
        intensity_constant = config.getfloat('rain', 'intensity_constant')
        rain_output_file_path = os.path.join(common.get_path_for('grass_db', config), 'constant_rain.tif')
        rain_dem, rain_profile = wdm.single_dem(input_file_for_depth_and_constant_rain)
        wdm.write_watermap_to_output(np.ones_like(rain_dem) * intensity_constant, rain_output_file_path, rain_profile)

    else:
        ascii_files_path = config.get('rain', 'ascii_files_path')
        GTiff_files_path = config.get('rain', 'GTiff_files_path')

        try:
            xmin = config.getfloat('rain', 'x_min')
            ymax = config.getfloat('rain', 'y_max')
            EPSG_code = config.getint('rain', 'EPSG_code')
        except ValueError as e:
            raise ValueError('xmin, ymax and EPSG_code are mandatory') from e

        xmax = config.get('rain', 'x_max')
        ymin = config.get('rain', 'y_min')
        xres = config.get('rain', 'x_res')
        yres = config.get('rain', 'y_res')
        xrotation = config.get('rain', 'x_rotation')
        yrotation = config.get('rain', 'y_rotation')

        xmax = float(xmax) if xmax else None
        ymin = float(ymin) if ymin else None
        xres = float(xres) if xres else None
        yres = float(yres) if yres else None
        xrotation = float(xrotation) if xrotation else None
        yrotation = float(yrotation) if yrotation else None

        if GTiff_files_path:
            rain_crop_search_criteria = config.get('rain', 'rain_crop_search_criteria')

        elif ascii_files_path:

            print('\nRain data is in ascii format and its being converted to Tif format...\n')

            GTiff_files_path = os.path.join(common.get_path_for('temporary', config), 'rain')
            ascii_search_criteria = config.get('rain', 'ascii_search_criteria')

            utl.ascii_to_geotiff(ascii_files_path, GTiff_files_path, xmin, ymax,
                                search_criteria=ascii_search_criteria, CRS_code=EPSG_code,
                                xmax=xmax, xrotation=xrotation, xres=xres,
                                yrotation=yrotation, ymin=ymin, yres=yres)

            print('\nDONE\n')

        else:

            raise ValueError('No rain data path provided')

        if config.get('rain', 'rain_relocation'): # BUG: Any string is True, not just "True"

            print('\nRain data is being relocated...\n')

            X_target = config.getfloat('rain', 'X_target')
            Y_target = config.getfloat('rain', 'Y_target')
            X_radar_rain = config.getfloat('rain', 'X_radar_rain')
            Y_radar_rain = config.getfloat('rain', 'Y_radar_rain')
            relocation_search_criteria = config.get('rain', 'relocation_search_criteria')

            utl.rain_relocation(GTiff_files_path, xmin, ymax, X_target, Y_target, X_radar_rain,
                                Y_radar_rain, search_criteria=relocation_search_criteria,CRS_code=EPSG_code,
                                xmax=xmax, xrotation=xrotation, xres=xres,
                                yrotation=yrotation, ymin=ymin, yres=yres)

            rain_crop_search_criteria = '*_relocated.tif'

        print('\nDONE\n')

        #* ---
        #* 5. Let's check if the rain rasters have same projection as raster of interest
        #* ---

        reference_file_path = os.path.join(common.get_path_for('grass_db', config), 'DEM_cropped.tif')

        if Path(GTiff_files_path).suffix == '':
            source_fname = next(glob.iglob(f"{GTiff_files_path}/*.tif"))  # we need only one file
            reproj = utl.raster_check_projection(source_fname, ref_fp = reference_file_path)

        else:

            reproj = utl.raster_check_projection(GTiff_files_path, ref_fp = reference_file_path)

        if reproj:

            root_reproj_rain_file = os.path.join(common.get_path_for('temporary', config), 'rain')

            print('\nRain file(s) being reprojected to match DEMs projection...\n')
            utl.raster_reproject(GTiff_files_path, reproj_fp=root_reproj_rain_file,
                                ref_fp=reference_file_path, search_criteria=rain_crop_search_criteria)

            rain_crop_search_criteria = "*_reproj.tif"

        if config.getboolean('rain', 'rain_crop_boolean'):

            print('\nRain files are being cropped...\n')

            rain_cropped_file_path = os.path.join(common.get_path_for('grass_db', config), 'rain')

            utl.raster_crop(GTiff_files_path, cropped_file_path=rain_cropped_file_path,
                            search_criteria=rain_crop_search_criteria, vector_file=crop_vector_file,
                            polygon_coords=crop_polygon_coords, polygon_crs=crop_polygon_crs,
                            mask_value=null_value, mask_all_touched=True)

            print('\nDONE\n')

    #* ---
    #* 6. We now turn to landuse information to define: friction, losses and infiltration.
    #*    Landuse information for Finland downloaded from:
    #*    https://www.avoindata.fi/data/fi/dataset/corine-maanpeite-2018
    #*    The file provides the different classes of the terrain with 20x20 m resolution
    #*
    #*    Let's first check whether if the landcover map has the same projection as our DEM raster
    #*    *As I already ran it once, reprojection was needed and I will skip it as it takes some minutes
    #*     to reproject.
    #*    But if its for the first time, do it!
    #* ---

    if config.getboolean('merging', 'landcover_merge_boolean'):

        landcover_folder_path = config.get('merging', 'landcover_file_path_to_merge')
        merged_file_path = os.path.join(common.get_path_for('temporary', config), 'landcover_merged.tif')

        lsc = config.get('merging', 'landcover_merge_search_criteria')

        merge_search_criteria = lsc.split()

        utl.raster_merge(landcover_folder_path, merged_file_path, search_criteria=merge_search_criteria)

        landcover_to_crop_path = merged_file_path

    else:

        landcover_to_crop_path = config.get('cropping', 'landcover_file_path_to_crop')

    # Let's check that the projection of the landcover map matches the projection of DEM map

    root_reproj_landcover_file = os.path.join(common.get_path_for('temporary', config), 'landcover_reproj.tif')

    reproj = utl.raster_check_projection(landcover_to_crop_path, ref_fp = dem_cropped_file_path)

    if reproj:

        print('\nReprojection needed for Landcover file to match DEM projection\n')
        utl.raster_reproject(landcover_to_crop_path, reproj_fp=root_reproj_landcover_file,
                            ref_fp=dem_cropped_file_path)
        landcover_to_crop_path = root_reproj_landcover_file

    #* ---
    #* 7. Let's crop the landcover file before extracting the information we need
    #* ---

    if config.getboolean('cropping', 'landcover_crop_boolean'):
        landcover_cropped_file_path = os.path.join(common.get_path_for('grass_db', config), 'landcover_cropped.tif')
        crop_search_criteria = config.get('cropping', 'landcover_crop_search_criteria')

        utl.raster_crop(landcover_to_crop_path, cropped_file_path=landcover_cropped_file_path,
                        search_criteria=crop_search_criteria, vector_file=crop_vector_file,
                        polygon_coords=crop_polygon_coords, polygon_crs=crop_polygon_crs,
                        mask_value=null_value, mask_all_touched=True)

    #* ---
    #* 8. Now we intialize the Landcover class and extract friction, losses and infiltration
    #* ---

    landcover = LC.Landcover(landcover_cropped_file_path)

    landcover.get_friction()

    landcover.get_losses()

    if config.getboolean('imperviousness', 'imperviousness_file_boolean'):

        if config.getboolean('merging', 'imperviousness_merge_boolean'):

            imperviousness_folder_path = config.get('merging', 'imperviousness_file_path_to_merge')
            merged_file_path = os.path.join(common.get_path_for('temporary', config), 'imperviousness_merged.tif')

            isc = config.get('merging', 'imperviousness_merge_search_criteria')

            merge_search_criteria = isc.split()

            utl.raster_merge(imperviousness_folder_path, merged_file_path,
                            search_criteria=merge_search_criteria)

            imperviousness_file_path = merged_file_path

        else:

            imperviousness_file_path = config.get('cropping', 'imperviousness_file_path_to_crop')

        root_reproj_imperviousness_file = os.path.join(common.get_path_for('temporary', config), 'imperviousness_reproj.tif')

        reproj = utl.raster_check_projection(imperviousness_file_path, ref_fp=dem_cropped_file_path)

        if reproj:

            print('\nImperviousness file is being reprojected to match DEMs projection...\n')
            utl.raster_reproject(imperviousness_file_path, reproj_fp=root_reproj_imperviousness_file,
                                ref_fp=dem_cropped_file_path)
            imperviousness_file_path = root_reproj_imperviousness_file

        if config.getboolean('cropping', 'imperviousness_crop_boolean'):

            imperviousness_cropped_file_path = os.path.join(common.get_path_for('grass_db', config),
                                                            'imperviousness_cropped.tif')
            crop_search_criteria = config.get('cropping', 'imperviousness_crop_search_criteria')

            utl.raster_crop(imperviousness_file_path,cropped_file_path=imperviousness_cropped_file_path,
                            search_criteria=crop_search_criteria, vector_file=crop_vector_file,
                            polygon_coords=crop_polygon_coords, polygon_crs=crop_polygon_crs,
                            mask_value=null_value, mask_all_touched=True)

            imperviousness_file_path = imperviousness_cropped_file_path

            utl.set_raster_resolution(imperviousness_file_path, landcover_cropped_file_path,
                                    binary=False, mask_value=241)

        rain_file_path = os.path.join(common.get_path_for('grass_db', config), 'rain')

        if config.getboolean('rain', 'infiltration_rate'):

            infiltration_output_path = os.path.join(common.get_path_for('grass_db', config), 'infiltration')

        else:

            infiltration_output_path = common.get_path_for('grass_db', config)

        landcover.get_infiltration(imperviousness_file_path, rain_file_path,
                                output_folder=infiltration_output_path,
                                infiltration_rate=config.getboolean('rain', 'infiltration_rate'))

    print('\nDONE with DEM, LANDCOVER and RAIN rasters\n')


if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)