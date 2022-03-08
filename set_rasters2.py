import sys, glob
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
    (Path(config_dict['folders']['temporary_files']) / 'rain').mkdir(exist_ok=True)
    (Path(config_dict['grass_info']['grass_db']) / 'rain').mkdir(exist_ok=True)
    (Path(config_dict['grass_info']['grass_db']) / 'infiltration').mkdir(exist_ok=True)
    print('Done.')

def merge_DEMs(config):
    merge_is_used = config['merging']['DEM_merge_boolean']
    raster_to_crop_path = config['cropping']['DEM_file_path_to_crop']
    if merge_is_used:
        merged_file_path = config['folders']['temporary_files'] + '/DEM_merged.tif'
        utl.raster_merge(config['folders']['DEM_files'],
                         merged_file_path, 
                         search_criteria=config['merging']['DEM_merge_search_criteria'])
        raster_to_crop_path = merged_file_path
    
    return merge_is_used,merged_file_path,raster_to_crop_path

def crop_DEMs(config, raster_to_crop_path):
    crop_is_used = config['cropping']['DEM_crop_boolean']
    if crop_is_used:
        dem_cropped_file_path = config['grass_info']['grass_db'] + '/DEM_cropped.tif'
        utl.raster_crop(raster_to_crop_path, 
                        config['cropping'],
                        cropped_fp=dem_cropped_file_path,
                        search_criteria=config['cropping']['DEM_crop_search_criteria'])
                        
    return crop_is_used,dem_cropped_file_path

def set_up_rain_rasters(config):
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
    return GTiff_files_path,rain_crop_search_criteria

def reproject_rain_rasters(config, GTiff_fp, rain_crop_sc):
    ref_fp = config['grass_info']['grass_db'] + '/DEM_cropped.tif'
    if Path(GTiff_fp).suffix == '':
        source_fname = next(glob.iglob(f"{GTiff_fp}/*.tif"))  # we need only one file
        reproj = utl.raster_check_projection(source_fname, ref_fp = ref_fp)
    else:
        reproj = utl.raster_check_projection(GTiff_fp, ref_fp = ref_fp)

    if reproj:
        root_reproj_rain_file = config['folders']['temporary_files'] + '/rain'
        print('\nRain file(s) being reprojected to match DEMs projection...\n')
        utl.raster_reproject(GTiff_fp, 
                             reproj_fp=root_reproj_rain_file,
                             ref_fp=ref_fp, 
                             search_criteria=rain_crop_sc)
        rain_crop_sc = "*_reproj.tif"
    return rain_crop_sc

def crop_rain_rasters(config, GTiff_fp, rain_crop_sc):
    if config['rain']['rain_crop_boolean']:
        print('\nRain files are being cropped...\n')
        rain_crop_fp = config['grass_info']['grass_db'] + '/rain'
        utl.raster_crop(GTiff_fp, config['cropping'], 
                        cropped_fp=rain_crop_fp,
                        search_criteria=rain_crop_sc)
        print('\nDONE\n')

def extract_rain_rasters(config, input_file):
    if config['rain']['constant']:
        intensity_c = config['rain']['intensity_constant']
        rain_output_fp = config['grass_info']['grass_db'] + '/constant_rain.tif'
        rain_dem, rain_profile = wdm.single_dem(input_file)
        wdm.write_watermap_to_output(np.ones_like(rain_dem) * intensity_c, rain_output_fp, rain_profile)
    else:
        GTiff_fp, rain_crop_sc = set_up_rain_rasters(config)
        #* 5. Let's check if the rain rasters have same projection as raster of interest
        rain_crop_sc = reproject_rain_rasters(config, GTiff_fp, rain_crop_sc)
        crop_rain_rasters(config, GTiff_fp, rain_crop_sc)

def main(config_file_name):
    config = cm.create_config_dictionary_from_config_file(config_file_name)
    add_missing_subfolders(config)
    #* 1. Downdload DEMs of interest:
    #*     - e.g. 2x2 DEM files from https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en
    #* 2. If the region of interest lays in the intersection of different DEMs, we can merge them into a
    #*    single GeoTiff file:
    merge_is_used, merge_fp, raster_to_crop_path = merge_DEMs(config)
    #* 3. We would probably like to do the analysis over a specific region and not the entire map.
    #*    Therefore, we proceed to crop the file according to a given polygon cooredinates or vector
    #*    shapefile.
    #*    If vector file not given, then both polygon coordinates and crs MUST be given
    crop_is_used, dem_cropped_fp = crop_DEMs(config, raster_to_crop_path)
    #* 3.1. Here a map describing the water depth in the initial timestep is created so that 
    # it corresponds to the cropped DEM above or just the dem if it is not cropped.
    depth_rain_file = dem_cropped_fp if crop_is_used else (merge_fp if merge_is_used else raster_to_crop_path)
    create_start_water_depth_file(depth_rain_file, config['grass_info']['grass_db'] + '/start_h.tif', config)
    #* 4. Now we proceed to extract the rain rasters
    extract_rain_rasters(config, depth_rain_file)
    #* 6. We now turn to landuse information to define: friction, losses and infiltration.
    #*    Landuse information for Finland downloaded from:
    #*    https://www.avoindata.fi/data/fi/dataset/corine-maanpeite-2018
    #*    The file provides the different classes of the terrain with 20x20 m resolution
    #*
    #*    Let's first check whether if the landcover map has the same projection as our DEM raster
    #*    *As I already ran it once, reprojection was needed and I will skip it as it takes some minutes
    #*     to reproject.
    #*    But if its for the first time, do it!
    landcover_to_crop_path = config['cropping']['landcover_file_path_to_crop']
    if config['merging']['landcover_merge_boolean']:
        merged_file_path = config['folders']['temporary_files'] + '/landcover_merged.tif'
        utl.raster_merge(config['merging']['landcover_file_path_to_merge'], 
                         merged_file_path, 
                         search_criteria=config['merging']['landcover_merge_search_criteria'])
        landcover_to_crop_path = merged_file_path
    # Let's check that the projection of the landcover map matches the projection of DEM map
    reproj = utl.raster_check_projection(landcover_to_crop_path, ref_fp = dem_cropped_fp)
    if reproj:
        print('\nReprojection needed for Landcover file to match DEM projection\n')
        utl.raster_reproject(landcover_to_crop_path, 
                             reproj_fp=config['folders']['temporary_files'] + '/landcover_reproj.tif',
                             ref_fp=dem_cropped_fp)
        landcover_to_crop_path = config['folders']['temporary_files'] + '/landcover_reproj.tif'
    #* 7. Let's crop the landcover file before extracting the information we need
    if config['cropping']['landcover_crop_boolean']:
        landcover_crop_fp = config['grass_info']['grass_db'] + '/landcover_cropped.tif'
        utl.raster_crop(landcover_to_crop_path, 
                        config['cropping'], 
                        cropped_fp=landcover_crop_fp,
                        search_criteria=config['cropping']['landcover_crop_search_criteria'])
    #* 8. Now we intialize the Landcover class and extract friction, losses and infiltration
    landcover = LC.Landcover(landcover_crop_fp)
    landcover.get_friction()
    landcover.get_losses()
    
    if config['imperviousness']['imperviousness_file_boolean']:
        imperviousness_fp = config['cropping']['imperviousness_file_path_to_crop']
        if config['merging']['imperviousness_merge_boolean']:
            imperviousness_fp = config['folders']['temporary_files'] + '/imperviousness_merged.tif'
            utl.raster_merge(config['merging']['imperviousness_file_path_to_merge'],
                             imperviousness_fp, 
                             search_criteria=config['merging']['imperviousness_merge_search_criteria'])

        reproj = utl.raster_check_projection(imperviousness_fp, ref_fp=dem_cropped_fp)
        if reproj:
            print('\nImperviousness file is being reprojected to match DEMs projection...\n')
            utl.raster_reproject(imperviousness_fp, 
                                 reproj_fp=config['folders']['temporary_files'] + '/imperviousness_reproj.tif',
                                 ref_fp=dem_cropped_fp)
            imperviousness_fp = config['folders']['temporary_files'] + '/imperviousness_reproj.tif'

        if config['cropping']['imperviousness_crop_boolean']:
            imperviousness_cropped_file_path = config['grass_info']['grass_db'] + '/imperviousness_cropped.tif'
            utl.raster_crop(imperviousness_fp,
                            config['cropping'],
                            cropped_fp=imperviousness_cropped_file_path,
                            search_criteria=config['cropping']['imperviousness_crop_search_criteria'])

            imperviousness_fp = imperviousness_cropped_file_path

            utl.set_raster_resolution(imperviousness_fp, landcover_crop_fp)
        
        infiltration_output_path = config['grass_info']['grass_db']
        if config['rain']['infiltration_rate']: 
            infiltration_output_path += '/infiltration'

        landcover.get_infiltration(imperviousness_fp, 
                                   config['grass_info']['grass_db'] + '/rain',
                                   output_folder=infiltration_output_path, 
                                   infiltration_rate=config['rain']['infiltration_rate'])

    print('\nDONE with DEM, LANDCOVER and RAIN rasters\n')


if __name__ == '__main__':
    try:
        config_file_name = sys.argv[1]
    except IndexError as e:
        raise RuntimeError(("No config file given! Please provide filename, e.g. "
                        "'GRASS_itzi_parameters.ini'")) from e
    main(config_file_name)