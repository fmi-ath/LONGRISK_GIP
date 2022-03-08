import os
import glob
from pathlib import Path

from osgeo import gdal, osr
import numpy as np
from modules.common import save_GTiff_raster

class Landcover:
    """Landcover object"""
    def __init__(self, root_landcover_file):
        """Initialize a Landcover object.

        Parameters
        ----------
        root_landcover_file : str
            Root path of the landcover tif file.
        ----------
        """
        dataset = gdal.Open(root_landcover_file, gdal.GA_ReadOnly)
        proj = osr.SpatialReference(wkt = dataset.GetProjection())
        band = dataset.GetRasterBand(1) # Note GetRasterBand() takes band no. starting from 1 not 0
        arr = band.ReadAsArray()

        self._root = Path(root_landcover_file).resolve()
        self._arr = arr
        self._geotransform = dataset.GetGeoTransform()
        self._crs_code = int(proj.GetAttrValue('AUTHORITY',1))
        self._crs = f'EPSG:{self._crs_code}'
        self._crs_wkt = proj.ExportToWkt()

        band = None
        dataset = None

    def get_friction(self):
        """Calculate friction at each pixel according to Landcover's file category.

        Store results to 'friction.tif'.
        """
        arr = self._arr
        new_arr = np.zeros(arr.shape)
        # Landcover category and its friction value
        singles = (
            (14, 0.04),
            (15, 0.08),
            (16, 0.03),
        )

        # Multiple subsequent categories have same value (inclusive range)
        ranges = (
            (1, 4, 0.08),
            (5, 7, 0.02),
            (8, 9, 0.05),
            (10, 13, 0.08),
            (17, 21, 0.04),
            (22, 29, 0.1),
            (30, 31, 0.04),
            (32, 36, 0.05),
            (37, 39, 0.04),
            (40, 45, 0.06),
            (46, 48, 0.03),
        )

        for value, new_value in singles:
            new_arr[arr == value] = new_value

        for start, stop, value in ranges:
            new_arr[(arr >= start) & (arr <= stop)] = value

        output_file = os.path.join(self._root.parent, 'friction.tif')
        save_GTiff_raster(self._crs_wkt, self._geotransform, new_arr, output_file)

    def get_losses(self):
        """Calculate losses at each pixel according to Landcover's file category.
        Store results to 'losses.tif'.
        """
        arr = self._arr
        # Only classes from 1 to 7 and 11 have value 10, others have value of 0
        new_arr = np.where(((arr >= 1) & (arr <= 7)) + (arr == 11), 10, 0)
        output_file = os.path.join(self._root.parent, 'losses.tif')
        save_GTiff_raster(self._crs_wkt, self._geotransform, new_arr, output_file)

    def get_infiltration(self, imperviousness_fp, rain_fp, output_folder,
                         infiltration_rate=True):
        """Calculates friction at each pixel according to Landcover's file category.

        Args:
            imperviousness_file_path (str): Raster containing imperviousness values
            rain_file_path (str): Folder or file containing rainfall values
            output_folder (str): Folder to store output files
            infiltration_rate (bool, optional): If True, calculate average infiltration from
                                                rainfall data. Defaults to True.

        Laskennassa käytetty valuntakerroin eri Corine2012-maankäyttöluokille
        *The runoff factor used in the calculation for different Corine2012 landuse categories*

        Runoff factor is affected by imperviousness of the surface

        Imperviousness taken from:
        https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness/status-maps/imperviousness-density-2018?tab=download

        Imperviousness tif has values from 0 to 100 indicating the percentage of imperviousness

        Infiltration would be 1 - runoff coefficient
        """
        arr = self._arr

        dataset = gdal.Open(imperviousness_fp, gdal.GA_ReadOnly)
        imperviousness_arr = dataset.GetRasterBand(1).ReadAsArray()
        dataset = None

        infiltration_arr = np.zeros(arr.shape)

        # This is about two orders of magnitude faster than the previous double loop
        infiltration_arr[(arr >= 1) & (arr <= 7)] = np.clip(imperviousness_arr[(arr >= 1) & (arr <= 7)], 65, 95) / 100  # pylint: disable=line-too-long
        infiltration_arr[(arr >= 8) & (arr <= 10)] = 0.05
        infiltration_arr[arr==11] = np.clip(imperviousness_arr[arr==11], 65, 95) / 100
        infiltration_arr[(arr >= 12) & (arr <= 15)] = 0.2
        infiltration_arr[(arr >= 16) & (arr <= 20)] = 0.2
        infiltration_arr[(arr >= 21) & (arr <= 36)] = 0.1
        infiltration_arr[(arr >= 37) & (arr <= 45)] = 0.05
        infiltration_arr[(arr >= 46) & (arr <= 48)] = 1

        infiltration_arr = 1 - infiltration_arr # Infiltration coefficients from runoff coefficients
        infiltration_arr = np.round(infiltration_arr, 2)  # Round for prettier output

        if not infiltration_rate:
            save_GTiff_raster(self._crs_wkt, self._geotransform, infiltration_arr, os.path.join(output_folder, 'infiltration.tif'))
            return

        if not Path(rain_fp).suffix:
            rain_files = sorted(glob.glob(os.path.join(rain_fp, '*.tif')))
        else:
            rain_files = sorted(glob.glob(rain_fp))

        for i, t in enumerate(rain_files):
            rain_dataset = gdal.Open(t, gdal.GA_ReadOnly)
            rain_arr = rain_dataset.GetRasterBand(1).ReadAsArray()
            rain_dataset = None

            av_rain = np.mean(rain_arr[rain_arr >= 0])

            infil_rate = infiltration_arr * av_rain

            output_file = os.path.join(output_folder, f'infiltration_{i:0>4}.tif')
            save_GTiff_raster(self._crs_wkt, self._geotransform, infil_rate, output_file)