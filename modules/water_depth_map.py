"""Create a water depth map for custom scenarios"""
from configparser import ConfigParser
from math import inf
from os import write

import numpy as np
import rasterio as rio
import rasterio.merge
import cv2 as cv


def single_dem(input_file):
    """Read single DEM file"""
    with rio.open(input_file, 'r') as src:
        profile = src.profile
        data = src.read(1)
    return data, profile

def merged_dem(input_rasters):
    """Read multiple DEMs and merge them into one"""
    array, affine = rasterio.merge.merge(input_rasters)
    array = array[0, :, :] # We are interested in only the first band

    rows, cols = array.shape

    # L4133A.tif's profile, but updated with the shape and affine of the merged raster
    profile = {
        'driver': 'GTiff',
        'dtype': 'float32',
        'nodata': -9999.0,
        'width': cols,
        'height': rows,
        'count': 1,
        'crs': "EPSG:3067",
        'transform': affine,
        'blockxsize': 256,
        'blockysize': 256,
        'tiled': True,
        'compress': 'lzw',
        'interleave': 'band',
    }

    return array, profile

def read_water_ini(configfile='./GRASS_itzi_parameters.ini'):
    """Read configuration file"""
    config = ConfigParser()
    config.read(configfile)
    return config

# WATER LEVEL GENERATION METHODS
# 1. Add constant water everywhere
# 2. Add water everywhere under certain height
    # 2.1 Relative to reference height (N2000)
    # 2.2 Relative to reference pixel's height
# 3. Add water with flood-fill algorithm from a point
# E. Add water to only certain landcover tiles (e.g. sea or lakes) with above methods

def generate_constant_water_depth_map(dem, depth):
    """Generate a uniform water depth map"""
    return np.ones_like(dem) * depth

def generate_height_restricted_water_depth_map(dem, depth, elevation):
    # Copy raster values
    modified = dem.copy()
    # Water surface level will be at reference height + added water depth
    surface_level = elevation + depth
    # Set values < reference to reference level
    modified[dem < surface_level] = surface_level
    # Calculate delta between modified and original levels
    # as that is how much water each pixel requires to have the water surface at
    # the defined level
    delta = modified - dem
    return delta

def generate_water_depth_map_using_floodfill(dem, coordinate, addon=0.005):
    elevation_at_seed=get_pixel_elevation(dem, coordinate)
    dem_floodable_points_by_elevation = np.where(dem <= elevation_at_seed, 255, 0).astype(np.uint8)
    possibly_flooded_points_mask = np.zeros(np.asarray(dem_floodable_points_by_elevation.shape)+2, dtype=np.uint8)
    
    # Here the possibly_flooded_points_mask is created by the flood fill algorithn
    cv.floodFill(dem_floodable_points_by_elevation, possibly_flooded_points_mask, coordinate, 255, flags=8)
    
    # The size of the mask is fitted to the size of the dem array
    mask_normalized = possibly_flooded_points_mask[1:-1, 1:-1]

    # Only remains to compute the amount of water that certain point should obtain.
    watermap = dem.copy()
    watermap[mask_normalized==1]=elevation_at_seed+addon
    watermap = watermap - dem
    watermap[mask_normalized!=1]=0
    
    return watermap


def get_pixel_elevation(dem, coordinate):
    # The coordinate might not be the same as array pointer.
    return dem[coordinate]

def main(dem, profile):
    # Get configuration
    config = read_water_ini()  # ConfigParser object

    config = config['water']  # SectionProxy (dict-like)

    # Algorithm
    method = config.get('method')
    if not method:
        raise ValueError('Missing method parameter!')

    depth = config.get('depth', fallback='0')
    depth = float(depth) if depth else 0

    if method == 'constant':
        # Early exit for constant depth watermap
        return generate_constant_water_depth_map(dem, depth)

    # Reference elevation for other methods
    x = config.get('x')
    y = config.get('y')
    x = int(x) if x else None
    y = int(y) if y else None
    
    if x and y:
        elevation = get_pixel_elevation(dem, (x,y))
    else:
        elevation = config.get('elevation', fallback='0')
        elevation = float(elevation) if elevation else 0

    # Other methods
    if method == 'limited':
        watermap = generate_height_restricted_water_depth_map(dem, depth, elevation)
    elif method == 'floodfill':
        if not x and y:
            raise ValueError(f'Floodfill requires seed coordinates x and y')
        addon=config.get('addon')
        addon=float(addon) if addon else 0

        watermap = generate_water_depth_map_using_floodfill(dem,(x,y),addon)
    else:
        raise ValueError(f'Unknown method {method}')

    return watermap

def write_watermap_to_output(watermap, output_file, profile) -> None:
    """Store watermap to disk or to GRASS GIS database as in TIF-format"""
    with rio.open(output_file, 'w', **profile) as dst:
        dst.write(watermap.astype(rio.uint8), 1)


if __name__ == '__main__':
    # MML DEMs have sea water levels != 0, and they may be different between images
    # -> artifical differences between image borders!
    FILE_PATH = __file__[:(__file__.rfind('/')+1)]
    # DEM file(s) for area of interest
    BASE_RASTERS = [
        FILE_PATH + '../data/DEM_2m_HKI/L4133A.tif',
        FILE_PATH + '../data/DEM_2m_HKI/L4133B.tif',
        FILE_PATH + '../data/DEM_2m_HKI/L4133C.tif',
        FILE_PATH + '../data/DEM_2m_HKI/L4133D.tif',
    ]

    # Read base DEM
    # dem, profile = (np.ones((5,5)), {})
    dem, profile = single_dem(FILE_PATH + '../../out/grassdata/DEM_cropped.tif')
    for i in range(len(dem)):
        print(dem[30,i]-dem[15,216],i)
    
    print(dem[15,216])
    print(dem[15,217])
    print(dem[15,215])
    print(dem[1,1])

    print(dem[15,216]-dem[15,216])
    print(dem[15,217]-dem[15,216])
    print(dem[15,215]-dem[15,216])