"""Common utility functions. This file does not import any other modules from
this folder.
"""
from pathlib import Path
import yaml
from osgeo import gdal
import numpy as np


def create_config_dictionary_from_config_file(config_filename):
    with open(config_filename) as file:
        config = yaml.full_load(file)
    return config

def ensure_folders_exist(config) -> None:
    """Create output folders if they don't exist"""
    for value in (('folders','temporary_files'), ('folders','processed_files'), ('folders','itzi_output_files'), ('grass_info','grass_db')):
        p = Path(config[value[0]][value[1]])
        if not p.suffix:
            p.mkdir(parents=True, exist_ok=True)

def save_GTiff_raster(projection_wkt, geotransform, array, 
                      GTiff_destination_file):
    nrows, ncols = np.shape(array)
    driver = gdal.GetDriverByName('GTiff')
    output_raster = driver.Create(GTiff_destination_file, ncols, nrows, 1 ,gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    output_raster.SetProjection(projection_wkt)
    output_raster.GetRasterBand(1).WriteArray(array)
    output_raster.FlushCache()
    output_raster = None