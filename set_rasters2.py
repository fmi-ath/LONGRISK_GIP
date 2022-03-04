import sys, os, glob
import numpy as np
import modules.water_depth_map as wdm
from pathlib import Path
import modules.Landcover as LC, modules.raster_utils as utl, modules.common as cm


def create_start_water_depth_file(input_file_path, output_file_path, config):
    if config['water']['method']:
        dem, profile = wdm.single_dem(input_file_path)
        watermap = wdm.main(dem, profile)
        wdm.write_watermap_to_output(watermap, output_file_path, profile)

def add_missing_subfolders(config_dict):
    print('Checking that output folders exist...', end=' ')
    cm.ensure_folders_exist(config_dict)
    # subfolders for this script
    (Path(config_dict['folders']['temporary_files']) / 'rain').mkdir(exist_ok=True)
    (Path(config_dict['grass_info']['grass_db']) / 'rain').mkdir(exist_ok=True)
    (Path(config_dict['grass_info']['grass_db']) / 'infiltration').mkdir(exist_ok=True)
    print('Done.')

def main(config_file_name):
    config = cm.create_config_dictionary_from_config_file(config_file_name)
    add_missing_subfolders(config)
    #* 1. Downdload DEMs of interest:
    #*     - e.g. 2x2 DEM files from https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en
    #* 2. If the region of interest lays in the intersection of different DEMs, we can merge them into a
    #*    single GeoTiff file:
    merge_is_used = config['merging']['DEM_merge_boolean']
    if merge_is_used:
        rasters_folder_path = Path(config['folders']['DEM_files'])
        merged_file_path = os.path.join(Path(config['folders']['temporary_files']), 'DEM_merged.tif')
        utl.raster_merge(rasters_folder_path, merged_file_path, search_criteria=config['merging']['DEM_merge_search_criteria'])
        raster_to_crop_path = merged_file_path
    else:
        raster_to_crop_path = config['cropping']['DEM_file_path_to_crop']
    #* 3. We would probably like to do the analysis over a specific region and not the entire map.
    #*    Therefore, we proceed to crop the file according to a given polygon cooredinates or vector
    #*    shapefile.
    #*
    #*    If vector file not given, then both polygon coordinates and crs MUST be given
    crop_is_used = config['cropping']['DEM_crop_boolean']
    if crop_is_used:
        dem_cropped_file_path = os.path.join(config['grass_info']['grass_db'], 'DEM_cropped.tif')
        utl.raster_crop(raster_to_crop_path, config['cropping'],cropped_file_path=dem_cropped_file_path,
                        search_criteria=config['cropping']['DEM_crop_search_criteria'])
    #* 3.1. Here a map describing the water depth in the initial timestep is created so that 
    # it corresponds to the cropped DEM above or just the dem if it is not cropped.
    input_file_for_depth_and_constant_rain = dem_cropped_file_path if crop_is_used else (merged_file_path if merge_is_used else raster_to_crop_path)
    create_start_water_depth_file(input_file_for_depth_and_constant_rain, os.path.join(config['grass_info']['grass_db'], 'start_h.tif'), config)
    #* 4. Now we proceed to extract the rain rasters
    if config['rain']['constant']:
        # create constant rain file
        intensity_constant = config['rain']['intensity_constant']
        rain_output_file_path = os.path.join(config['grass_info']['grass_db'], 'constant_rain.tif')
        rain_dem, rain_profile = wdm.single_dem(input_file_for_depth_and_constant_rain)
        wdm.write_watermap_to_output(np.ones_like(rain_dem) * intensity_constant, rain_output_file_path, rain_profile)
    else:
        GTiff_files_path = config['rain']['GTiff_files_path']

        if GTiff_files_path:
            rain_crop_search_criteria = config['rain']['rain_crop_search_criteria']
        elif config['rain']['ascii_files_path']:
            print('\nRain data is in ascii format and its being converted to Tif format...\n')
            GTiff_files_path = config['folders']['temporary_files']+'/rain'
            utl.ascii_to_geotiff(GTiff_files_path, config['rain'])
            print('\nDONE\n')
        else:
            raise ValueError('No rain data path provided')

        if config['rain']['rain_relocation']:
            print('\nRain data is being relocated...\n')
            utl.rain_relocation(GTiff_files_path, config['rain'])
            rain_crop_search_criteria = '*_relocated.tif'

        print('\nDONE\n')
        #* 5. Let's check if the rain rasters have same projection as raster of interest
        reference_file_path = os.path.join(config['grass_info']['grass_db'], 'DEM_cropped.tif')
        if Path(GTiff_files_path).suffix == '':
            source_fname = next(glob.iglob(f"{GTiff_files_path}/*.tif"))  # we need only one file
            reproj = utl.raster_check_projection(source_fname, ref_fp = reference_file_path)
        else:
            reproj = utl.raster_check_projection(GTiff_files_path, ref_fp = reference_file_path)

        if reproj:
            root_reproj_rain_file = os.path.join(config['folders']['temporary_files'], 'rain')
            print('\nRain file(s) being reprojected to match DEMs projection...\n')
            utl.raster_reproject(GTiff_files_path, reproj_fp=root_reproj_rain_file,
                                ref_fp=reference_file_path, search_criteria=rain_crop_search_criteria)
            rain_crop_search_criteria = "*_reproj.tif"

        if config['rain']['rain_crop_boolean']:
            print('\nRain files are being cropped...\n')
            rain_cropped_file_path = os.path.join(config['grass_info']['grass_db'], 'rain')
            utl.raster_crop(GTiff_files_path, config['cropping'], cropped_file_path=rain_cropped_file_path,
                            search_criteria=rain_crop_search_criteria)
            print('\nDONE\n')
    #* 6. We now turn to landuse information to define: friction, losses and infiltration.
    #*    Landuse information for Finland downloaded from:
    #*    https://www.avoindata.fi/data/fi/dataset/corine-maanpeite-2018
    #*    The file provides the different classes of the terrain with 20x20 m resolution
    #*
    #*    Let's first check whether if the landcover map has the same projection as our DEM raster
    #*    *As I already ran it once, reprojection was needed and I will skip it as it takes some minutes
    #*     to reproject.
    #*    But if its for the first time, do it!
    if config['merging']['landcover_merge_boolean']:
        landcover_folder_path = config['merging']['landcover_file_path_to_merge']
        merged_file_path = os.path.join(config['folders']['temporary_files'], 'landcover_merged.tif')
        utl.raster_merge(landcover_folder_path, merged_file_path, search_criteria=config['merging']['landcover_merge_search_criteria'])
        landcover_to_crop_path = merged_file_path
    else:
        landcover_to_crop_path = config['cropping']['landcover_file_path_to_crop']
    # Let's check that the projection of the landcover map matches the projection of DEM map
    reproj = utl.raster_check_projection(landcover_to_crop_path, ref_fp = dem_cropped_file_path)
    if reproj:
        print('\nReprojection needed for Landcover file to match DEM projection\n')
        utl.raster_reproject(landcover_to_crop_path, reproj_fp=os.path.join(config['folders']['temporary_files'], 'landcover_reproj.tif'),
                            ref_fp=dem_cropped_file_path)
        landcover_to_crop_path = os.path.join(config['folders']['temporary_files'], 'landcover_reproj.tif')
    #* 7. Let's crop the landcover file before extracting the information we need
    if config['cropping']['landcover_crop_boolean']:
        landcover_cropped_file_path = os.path.join(config['grass_info']['grass_db'], 'landcover_cropped.tif')
        utl.raster_crop(landcover_to_crop_path, config['cropping'], cropped_file_path=landcover_cropped_file_path,
                        search_criteria=config['cropping']['landcover_crop_search_criteria'])
    #* 8. Now we intialize the Landcover class and extract friction, losses and infiltration
    landcover = LC.Landcover(landcover_cropped_file_path)
    landcover.get_friction()
    landcover.get_losses()
    
    if config['imperviousness']['imperviousness_file_boolean']:
        if config['merging']['imperviousness_merge_boolean']:
            imperviousness_folder_path = config['merging']['imperviousness_file_path_to_merge']
            merged_file_path = os.path.join(config['folders']['temporary_files'], 'imperviousness_merged.tif')
            utl.raster_merge(imperviousness_folder_path, merged_file_path, search_criteria=config['merging']['imperviousness_merge_search_criteria'])
            imperviousness_file_path = merged_file_path
        else:
            imperviousness_file_path = config['cropping']['imperviousness_file_path_to_crop']

        root_reproj_imperviousness_file = os.path.join(config['folders']['temporary_files'], 'imperviousness_reproj.tif')
        reproj = utl.raster_check_projection(imperviousness_file_path, ref_fp=dem_cropped_file_path)
        if reproj:
            print('\nImperviousness file is being reprojected to match DEMs projection...\n')
            utl.raster_reproject(imperviousness_file_path, reproj_fp=root_reproj_imperviousness_file, ref_fp=dem_cropped_file_path)
            imperviousness_file_path = root_reproj_imperviousness_file

        if config['cropping']['imperviousness_crop_boolean']:
            imperviousness_cropped_file_path = os.path.join(config['grass_info']['grass_db'], 'imperviousness_cropped.tif')
            utl.raster_crop(imperviousness_file_path,config['cropping'],cropped_file_path=imperviousness_cropped_file_path,
                            search_criteria=config['cropping']['imperviousness_crop_search_criteria'])

            imperviousness_file_path = imperviousness_cropped_file_path

            utl.set_raster_resolution(imperviousness_file_path, landcover_cropped_file_path, binary=False, mask_value=241)

        if config['rain']['infiltration_rate']:
            infiltration_output_path = os.path.join(config['grass_info']['grass_db'], 'infiltration')
        else:
            infiltration_output_path = config['grass_info']['grass_db']

        landcover.get_infiltration(imperviousness_file_path, os.path.join(config['grass_info']['grass_db'], 'rain'),
                                output_folder=infiltration_output_path, infiltration_rate=config['rain']['infiltration_rate'])

    print('\nDONE with DEM, LANDCOVER and RAIN rasters\n')


if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)