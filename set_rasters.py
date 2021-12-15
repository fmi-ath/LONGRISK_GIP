import os
import glob
from pathlib import Path

import rasterio as rio

import modules.Landcover as LC
import modules.raster_utils as utl
from modules import common

#* ---
#* We read the input information from the ini file and then the magic happens
#* ---

# instantiate
config = common.CONFIG

# Path setup
temp_folder = common.get_path_for('temporary')
store_folder = common.get_path_for('processed')  # dump final modified files here
grassdata_folder = common.get_path_for('grass_db')

print('Checking that output folders exist...', end=' ')
common.ensure_folders_exist()
# subfolders for this script
(temp_folder / 'rain').mkdir(exist_ok=True)
(grassdata_folder / 'rain').mkdir(exist_ok=True)
(grassdata_folder / 'infiltration').mkdir(exist_ok=True)
print('Done.')

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

#* ---
#* 1. Downdload DEMs of interest:
#*     - e.g. 2x2 DEM files from https://tiedostopalvelu.maanmittauslaitos.fi/tp/kartta?lang=en
#* 2. If the region of interest lays in the intersection of different DEMs, we can merge them into a
#*    single GeoTiff file:
#* ---

if config.getboolean('merging', 'DEM_merge_boolean'):

    rasters_folder_path = common.get_path_for('DEM')
    merged_file_path = os.path.join(temp_folder, 'DEM_merged.tif')

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

if config.getboolean('cropping', 'DEM_crop_boolean'):

    dem_cropped_file_path = os.path.join(grassdata_folder, 'DEM_cropped.tif')
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
#* 4. Now we proceed to extract the rain rasters
#* ---

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

    GTiff_files_path = os.path.join(temp_folder, 'rain')
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

reference_file_path = os.path.join(grassdata_folder, 'DEM_cropped.tif')

if Path(GTiff_files_path).suffix == '':
    source_fname = next(glob.iglob(f"{GTiff_files_path}/*.tif"))  # we need only one file
    reproj = utl.raster_check_projection(source_fname, ref_fp = reference_file_path)

else:

    reproj = utl.raster_check_projection(GTiff_files_path, ref_fp = reference_file_path)

if reproj:

    root_reproj_rain_file = os.path.join(temp_folder, 'rain')

    print('\nRain file(s) being reprojected to match DEMs projection...\n')
    utl.raster_reproject(GTiff_files_path, reproj_fp=root_reproj_rain_file,
                         ref_fp=reference_file_path, search_criteria=rain_crop_search_criteria)

    rain_crop_search_criteria = "*_reproj.tif"

if config.getboolean('rain', 'rain_crop_boolean'):

    print('\nRain files are being cropped...\n')

    rain_cropped_file_path = os.path.join(grassdata_folder, 'rain')

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
    merged_file_path = os.path.join(temp_folder, 'landcover_merged.tif')

    lsc = config.get('merging', 'landcover_merge_search_criteria')

    merge_search_criteria = lsc.split()

    utl.raster_merge(landcover_folder_path, merged_file_path, search_criteria=merge_search_criteria)

    landcover_to_crop_path = merged_file_path

else:

    landcover_to_crop_path = config.get('cropping', 'landcover_file_path_to_crop')

# Let's check that the projection of the landcover map matches the projection of DEM map

root_reproj_landcover_file = os.path.join(temp_folder, 'landcover_reproj.tif')

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
    landcover_cropped_file_path = os.path.join(grassdata_folder, 'landcover_cropped.tif')
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
        merged_file_path = os.path.join(temp_folder, 'imperviousness_merged.tif')

        isc = config.get('merging', 'imperviousness_merge_search_criteria')

        merge_search_criteria = isc.split()

        utl.raster_merge(imperviousness_folder_path, merged_file_path,
                         search_criteria=merge_search_criteria)

        imperviousness_file_path = merged_file_path

    else:

        imperviousness_file_path = config.get('cropping', 'imperviousness_file_path_to_crop')

    root_reproj_imperviousness_file = os.path.join(temp_folder, 'imperviousness_reproj.tif')

    reproj = utl.raster_check_projection(imperviousness_file_path, ref_fp=dem_cropped_file_path)

    if reproj:

        print('\nImperviousness file is being reprojected to match DEMs projection...\n')
        utl.raster_reproject(imperviousness_file_path, reproj_fp=root_reproj_imperviousness_file,
                             ref_fp=dem_cropped_file_path)
        imperviousness_file_path = root_reproj_imperviousness_file

    if config.getboolean('cropping', 'imperviousness_crop_boolean'):

        imperviousness_cropped_file_path = os.path.join(grassdata_folder,
                                                        'imperviousness_cropped.tif')
        crop_search_criteria = config.get('cropping', 'imperviousness_crop_search_criteria')

        utl.raster_crop(imperviousness_file_path,cropped_file_path=imperviousness_cropped_file_path,
                        search_criteria=crop_search_criteria, vector_file=crop_vector_file,
                        polygon_coords=crop_polygon_coords, polygon_crs=crop_polygon_crs,
                        mask_value=null_value, mask_all_touched=True)

        imperviousness_file_path = imperviousness_cropped_file_path

        utl.set_raster_resolution(imperviousness_file_path, landcover_cropped_file_path,
                                  binary=False, mask_value=241)

    rain_file_path = os.path.join(grassdata_folder, 'rain')

    if config.getboolean('rain', 'infiltration_rate'):

        infiltration_output_path = os.path.join(grassdata_folder, 'infiltration')

    else:

        infiltration_output_path = grassdata_folder

    landcover.get_infiltration(imperviousness_file_path, rain_file_path,
                               output_folder=infiltration_output_path,
                               infiltration_rate=config.getboolean('rain', 'infiltration_rate'))

print('\nDONE with DEM, LANDCOVER and RAIN rasters\n')
