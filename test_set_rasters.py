import os
import unittest
import numpy as np
import set_rasters2 as sr
import rasterio as rio
import rasterio.crs as crs
from rasterio.transform import from_origin
# These tests have to be ran inside the Docker container 
# where the app runs. Run run_docker.sh to achieve this.

class test_set_rasters(unittest.TestCase):
    setup_done = False
    config_filename = "test_GRASS_itzi_parameters.ini"
    test_config_dict = dict()
    dem_location = ""
    mock_raster = np.array([[]])

    def setUp(self):
        if not test_set_rasters.setup_done:
            self.create_mock_DEM_file()
            test_set_rasters.test_config_dict = self.get_test_config_dict()
            test_set_rasters.setup_done = True
        
    def create_mock_DEM_file(self):
        test_set_rasters.mock_raster = np.random.randint(5, size = (21, 10))
        test_set_rasters.dem_location = "mock_DEM.tif"

        mock_DEM_data = rio.open(test_set_rasters.dem_location, 'w', driver='GTiff',
                                    height=self.mock_raster.shape[0], width=self.mock_raster.shape[1],
                                    count=1, dtype='float32',
                                    nodata=-9999.0, crs=crs.CRS.from_epsg(3067),
                                    transform=rio.Affine(2.0, 0.0, 362000.0,0.0, -2.0, 6672000.0))

        mock_DEM_data.write(test_set_rasters.mock_raster, 1)
        mock_DEM_data.close()

    def get_test_config_dict(self):
        return dict()

    def test_get_a_dictionary_from_create_config_dictionary(self):
        config_dict = sr.create_config_dictionary_from_config_file(self.config_filename)
        self.assertIsInstance(type(config_dict), type(dict))

    def test_get_a_dictionary_of_folders_from_created_config_dict(self):
        config_dict = sr.create_config_dictionary_from_config_file(self.config_filename)
        folders_dict = config_dict['folders']
        self.assertIsInstance(type(folders_dict), type(dict))

    def test_get_dem_files_folder_from_config(self):
        config_dict = sr.create_config_dictionary_from_config_file(self.config_filename)
        self.assertEqual(config_dict['folders']['DEM_files'], "/data/DEM_2m_HKI")

    def assert_equal_two_rasters(self, test_raster, true_raster):
        self.assertEqual(len(test_raster), len(true_raster))
        length_of_first = len(true_raster[0])
        
        for row in range(len(test_raster)):
            self.assertEqual(len(test_raster[row]), length_of_first)
            self.assertEqual(len(true_raster[row]), length_of_first)
            
            for column in range(len(test_raster[row])):
                self.assertEqual(test_raster[row][column], true_raster[row][column])


if __name__ == '__main__':
    unittest.main()