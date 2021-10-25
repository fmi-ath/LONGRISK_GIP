import os
import glob
from pathlib import Path

from osgeo import gdal, osr
import numpy as np

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

        dataset = gdal.Open(root_landcover_file, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
        proj = osr.SpatialReference(wkt = dataset.GetProjection())
        geotransform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        arr = band.ReadAsArray()

        self._root = Path(root_landcover_file).resolve() # Resolve to remove relative paths
        self._arr = arr
        self._shape = arr.shape
        self._geotransform = geotransform
        self._crs_code = int(proj.GetAttrValue('AUTHORITY',1)) # it returns a string consisting of the EPSG number
        self._crs = f'EPSG:{self._crs_code}'

    def get_friction(self):
        """Calculate friction at each pixel according to Landcover's file category.

        Store results to 'friction.tif'.
        """

        arr = self._arr

        new_arr = np.zeros(self._shape)

        # TODO: See if this double-loop can be rewritten using np.nditer function or indexing tricks
        for i in range(self._shape[0]):

            for j in range(self._shape[1]):

                element = arr[i, j]

                if 1 <= element <= 4:
                    friction_value = 0.08

                elif 5 <= element <= 7:
                    friction_value = 0.02

                elif 8 <= element <= 9:
                    friction_value = 0.05

                elif 10 <= element <= 13:
                    friction_value = 0.08

                elif element == 14:
                    friction_value = 0.04

                elif element == 15:
                    friction_value = 0.08

                elif element == 16:
                    friction_value = 0.03

                elif 17 <= element <= 21:
                    friction_value = 0.04

                elif 22 <= element <= 29:
                    friction_value = 0.1

                elif 30 <= element <= 31:
                    friction_value = 0.04

                elif 32 <= element <= 36:
                    friction_value = 0.05

                elif 37 <= element <= 39:
                    friction_value = 0.04

                elif 40 <= element <= 45:
                    friction_value = 0.06

                elif 46 <= element <= 48:
                    friction_value = 0.03

                else:
                    continue  # Skip unknown values

                new_arr[i, j] = friction_value

        #geotransform = dataset.GetGeoTransform()
        output_file = os.path.join(self._root.parent, 'friction.tif')
        self._write_raster_file(output_file, new_arr)


    def get_losses(self):
        """Calculate losses at each pixel according to Landcover's file category.

        Store results to 'losses.tif'.
        """

        arr = self._arr

        new_arr = np.where(((arr >= 1) & (arr <= 7)) + (arr == 11), 10, 0)

        output_file = os.path.join(self._root.parent, 'losses.tif')
        self._write_raster_file(output_file, new_arr)


    def get_infiltration(self, imperviousness_file_path, rain_file_path, output_folder, infiltration_rate = True):
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

        dataset = gdal.Open(imperviousness_file_path, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
        band = dataset.GetRasterBand(1)
        imperviousness_arr = band.ReadAsArray()

        infiltration_arr = np.zeros(self._shape)

        # TODO: I wonder if this is faster than the previous double loop? Here we have ~18 numpy loops...
        infiltration_arr[(arr >= 1) & (arr <= 7)] = np.clip(imperviousness_arr[(arr >= 1) & (arr <= 7)], 65, 95) / 100
        infiltration_arr[(arr >= 8) & (arr <= 10)] = 0.05
        infiltration_arr[arr==11] = np.clip(imperviousness_arr[arr==11], 65, 95) / 100
        infiltration_arr[(arr >= 12) & (arr <= 15)] = 0.2
        infiltration_arr[(arr >= 16) & (arr <= 20)] = 0.2
        infiltration_arr[(arr >= 21) & (arr <= 36)] = 0.1
        infiltration_arr[(arr >= 37) & (arr <= 45)] = 0.05
        infiltration_arr[(arr >= 46) & (arr <= 48)] = 0.05

        infiltration_arr = 1 - infiltration_arr # Infiltration coefficients from runoff coefficients

        if not infiltration_rate:
            output_file = os.path.join(output_folder, 'infiltration.tif')
            self._write_raster_file(output_file, infiltration_arr)

            return

        if not Path(rain_file_path).suffix:
            rain_files = sorted(glob.glob(os.path.join(rain_file_path, '*.tif')))
        else:
            rain_files = sorted(glob.glob(rain_file_path))

        for i, t in enumerate(rain_files):
            rain_dataset = gdal.Open(t, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
            rainband = rain_dataset.GetRasterBand(1)
            rain_arr = rainband.ReadAsArray()

            av_rain = np.mean(rain_arr[rain_arr >= 0])

            infil_rate = infiltration_arr * av_rain

            #file_name = os.path.splitext(os.path.split(t)[1])[0]
            #file_time = file_name[4:8]

            output_file = os.path.join(output_folder, f'infiltration_{i:0>4}.tif')
            self._write_raster_file(output_file, infil_rate)


    def _write_raster_file(self, filename, array, dtype=gdal.GDT_Float32):
        """Convenience function to write a calculated raster to a GeoTiff file.

        Args:
            filename (str): Path to the output file
            array (numpy.array): Raster image to-be-written
            dtype (gdal constant, optional): GDAL Data type for the file. Defaults to gdal.GDT_Float32.
        """
        driver = gdal.GetDriverByName('GTiff')

        raster = driver.Create(filename, self._shape[1], self._shape[0], 1, dtype)
        raster.SetGeoTransform(self._geotransform)  # Specify its coordinates

        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(self._crs_code)

        raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system to the file
        raster.GetRasterBand(1).WriteArray(array)   # Writes my array to the raster

        raster.FlushCache()
        raster = None
