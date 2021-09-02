import os
import glob
from pathlib import Path

from osgeo import gdal, osr
import numpy as np
import modules.utils as utl

class Landcover:

    def __init__(self, root_landcover_file):
        """Initialize a Landcover object.

        Parameters
        ----------
        root_radar_file : str
            Root path of the landcover tif file.
        ----------
        """

        dataset = gdal.Open(root_landcover_file, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
        proj = osr.SpatialReference(wkt = dataset.GetProjection())
        geotransform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        arr = band.ReadAsArray()

        self._root_ = Path(root_landcover_file).resolve() # Resolve to remove relative paths
        self._arr_ = arr
        self._shape_ = arr.shape
        self._geotransform_ = geotransform
        self._crs_code_ = int(proj.GetAttrValue('AUTHORITY',1)) # it returns a string consisting of the EPSG number
        self._crs_ = f'EPSG:{self._crs_code_}'

    def get_friction(self):

        """Calcualtes friction at each pixel according to Landcover's file category.

            Parameters
            ----------


            Returns
            -------
            out : file
                Friction file
        """

        arr = self._arr_

        new_arr = np.zeros(self._shape_)

        for i in range(self._shape_[0]):

            for j in range(self._shape_[1]):

                if (arr[i, j] >= 1) & (arr[i, j] <= 4):

                    new_arr[i, j] = 0.08

                elif (arr[i, j] >= 5) & (arr[i, j] <= 7):

                    new_arr[i, j] = 0.02

                elif (arr[i, j] >= 8) & (arr[i, j] <= 9):

                    new_arr[i, j] = 0.05

                elif (arr[i, j] >= 10) & (arr[i, j] <= 13):

                    new_arr[i, j] = 0.08

                elif arr[i, j] == 14:

                    new_arr[i, j] = 0.04

                elif arr[i, j] == 15:

                    new_arr[i, j] = 0.08

                elif arr[i, j] == 16:

                    new_arr[i, j] = 0.03

                elif (arr[i, j] >= 17) & (arr[i, j] <= 21):

                    new_arr[i, j] = 0.04

                elif (arr[i, j] >= 22) & (arr[i, j] <= 29):

                    new_arr[i, j] = 0.1

                elif (arr[i, j] >= 30) & (arr[i, j] <= 31):

                    new_arr[i, j] = 0.04

                elif (arr[i, j] >= 32) & (arr[i, j] <= 36):

                    new_arr[i, j] = 0.05

                elif (arr[i, j] >= 37) & (arr[i, j] <= 39):

                    new_arr[i, j] = 0.04

                elif (arr[i, j] >= 40) & (arr[i, j] <= 45):

                    new_arr[i, j] = 0.06

                elif (arr[i, j] >= 46) & (arr[i, j] <= 48):

                    new_arr[i, j] = 0.03

        #geotransform = dataset.GetGeoTransform()
        output_file = os.path.join(os.path.split(self._root_)[0], 'friction.tif')

        output_raster = gdal.GetDriverByName('GTiff').Create(output_file, self._shape_[1], self._shape_[0], 1 ,gdal.GDT_Float32)  # Open the file
        output_raster.SetGeoTransform(self._geotransform_)  # Specify its coordinates
        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(self._crs_code_)

        output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                           # to the file
        output_raster.GetRasterBand(1).WriteArray(new_arr)   # Writes my array to the raster

        output_raster.FlushCache()

        return

    def get_losses(self):

        """Calcualtes losses at each pixel according to Landcover's file category.

            Parameters
            ----------


            Returns
            -------
            out : file
                Losses file
        """

        arr = self._arr_

        new_arr = np.where(((arr >= 1) & (arr <= 7)) + (arr == 11), 10, 0)

        output_file = os.path.join(os.path.split(self._root_)[0], 'losses.tif')

        #geotransform = dataset.GetGeoTransform()
        output_raster = gdal.GetDriverByName('GTiff').Create(output_file, self._shape_[1], self._shape_[0], 1 ,gdal.GDT_Float32)  # Open the file
        output_raster.SetGeoTransform(self._geotransform_)  # Specify its coordinates
        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(self._crs_code_)                     # This one specifies WGS84 lat long.
                                                     # Anyone know how to specify the
                                                     # IAU2000:49900 Mars encoding?
        output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                           # to the file
        output_raster.GetRasterBand(1).WriteArray(new_arr)   # Writes my array to the raster

        output_raster.FlushCache()

        return

    def get_infiltration(self, imperviousness_file_path, rain_file_path, output_folder, infiltration_rate = True):

        """Calcualtes friction at each pixel according to Landcover's file category.
            Laskennassa käytetty valuntakerroin eri Corine2012-maankäyttöluokille
            *The runoff factor used in the calculation for different Corine2012 landuse categories*

            Runoff factor is affected by imperviousness of the surface

            Imperviousness taken from:
            https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness/status-maps/imperviousness-density-2018?tab=download

            Imperviousness tif has values from 0 to 100 indicating the percentage of imperviousness

            Infiltration would be 1 - runoff coefficient

            Parameters
            ----------


            Returns
            -------
            out : file
                Friction file
        """

        arr = self._arr_

        dataset = gdal.Open(imperviousness_file_path, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
        proj = osr.SpatialReference(wkt = dataset.GetProjection())
        CRS_code = proj.GetAttrValue('AUTHORITY',1)

        geotransform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        imperviousness_arr = band.ReadAsArray()

        infiltration_arr = np.zeros(self._shape_)

        for i in range(self._shape_[0]):

            for j in range(self._shape_[1]):

                if (arr[i, j] >= 1) & (arr[i, j] <= 7):

                    if imperviousness_arr[i, j] < 65:

                        coef = 0.65

                    elif imperviousness_arr[i, j] > 95:

                        coef = 0.95

                    else:

                        coef = imperviousness_arr[i, j] / 100

                    infiltration_arr[i, j] = coef

                elif (arr[i, j] >= 8) & (arr[i, j] <= 10):

                    infiltration_arr[i, j] = 0.05

                elif arr[i, j] == 11:

                    if imperviousness_arr[i, j] < 65:

                        coef = 0.65

                    elif imperviousness_arr[i, j] > 95:

                        coef = 0.95

                    else:

                        coef = imperviousness_arr[i, j] / 100

                    infiltration_arr[i, j] = coef

                elif (arr[i, j] >= 12) & (arr[i, j] <= 15):

                    infiltration_arr[i, j] = 0.2

                elif (arr[i, j] >= 16) & (arr[i, j] <= 20):

                    infiltration_arr[i, j] = 0.2

                elif (arr[i, j] >= 21) & (arr[i, j] <= 36):

                    infiltration_arr[i, j] = 0.1

                elif (arr[i, j] >= 37) & (arr[i, j] <= 45):

                    infiltration_arr[i, j] = 0.05

                elif (arr[i, j] >= 46) & (arr[i, j] <= 48):

                    infiltration_arr[i, j] = 1

        infiltration_arr = 1 - infiltration_arr # Infiltration coefficients from runoff coefficients

        if infiltration_rate:

            if not Path(rain_file_path).suffix:

                rain_files = sorted(glob.glob(os.path.join(rain_file_path, '*.tif')))

            else:

                rain_files = sorted(glob.glob(rain_file_path))

        else:

            output_file = os.path.join(output_folder, 'infiltration.tif')

            output_raster = gdal.GetDriverByName('GTiff').Create(output_file, self._shape_[1], self._shape_[0], 1 ,gdal.GDT_Float32)  # Open the file
            output_raster.SetGeoTransform(self._geotransform_)  # Specify its coordinates
            srs = osr.SpatialReference()                 # Establish its coordinate encoding
            srs.ImportFromEPSG(self._crs_code_)                     # This one specifies WGS84 lat long.
                                                         # Anyone know how to specify the
                                                         # IAU2000:49900 Mars encoding?
            output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                               # to the file
            output_raster.GetRasterBand(1).WriteArray(infiltration_arr)   # Writes my array to the raster

            output_raster.FlushCache()

            return

        i = 0

        for t in rain_files:

            rain_dataset = gdal.Open(t, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
            rainband = rain_dataset.GetRasterBand(1)
            rain_arr = rainband.ReadAsArray()

            av_rain = np.mean(rain_arr[rain_arr >= 0])

            infil_rate = infiltration_arr * av_rain

            #file_name = os.path.splitext(os.path.split(t)[1])[0]
            #file_time = file_name[4:8]

            #geotransform = dataset.GetGeoTransform()
            output_file = os.path.join(output_folder, 'infiltration_' + str(i).zfill(4) + '.tif')

            output_raster = gdal.GetDriverByName('GTiff').Create(output_file, self._shape_[1], self._shape_[0], 1 ,gdal.GDT_Float32)  # Open the file
            output_raster.SetGeoTransform(self._geotransform_)  # Specify its coordinates
            srs = osr.SpatialReference()                 # Establish its coordinate encoding
            srs.ImportFromEPSG(self._crs_code_)                     # This one specifies WGS84 lat long.
                                                         # Anyone know how to specify the
                                                         # IAU2000:49900 Mars encoding?
            output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                               # to the file
            output_raster.GetRasterBand(1).WriteArray(infil_rate)   # Writes my array to the raster

            output_raster.FlushCache()

            i += 1

        return
