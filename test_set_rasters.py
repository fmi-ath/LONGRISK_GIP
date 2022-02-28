import os
import unittest
import numpy as np
import set_rasters2 as sr
import rasterio as rio
import rasterio.crs as crs
from rasterio.transform import from_origin
import osgeo
# import yaml
from modules import common
from pathlib import Path

# These tests have to be ran inside the Docker container 
# where the app runs. Run run_docker.sh to achieve this.

class test_set_rasters(unittest.TestCase):
    setup_done = False
    config_filename_yml = "test_GRASS_itzi_parameters.yml"
    test_config_dict = dict()
    dem_locations = []
    mock_rasters = []

    def setUp(self):
        if not test_set_rasters.setup_done:
            self.create_mock_DEM_files()
            test_set_rasters.setup_done = True
        
    def create_mock_DEM_files(self):
        test_set_rasters.mock_rasters = [np.random.randint(5, size = (21, 10)), np.random.randint(5, size = (21, 10))]
        test_set_rasters.dem_locations = ["mock_DEM.tif", "mock_2_DEM.tif"]        

        for i in [0,1]:
            mock_DEM_data = rio.open(test_set_rasters.dem_locations[i], 'w', driver='GTiff',
                                    height=self.mock_rasters[i].shape[0], width=self.mock_rasters[i].shape[1],
                                    count=1, dtype='float32',
                                    nodata=-9999.0, crs=crs.CRS.from_epsg(3067),
                                    transform=rio.Affine(2.0, 0.0, 362000.0,0.0, -2.0, 6672000.0))

            mock_DEM_data.write(test_set_rasters.mock_rasters[i], 1)
            mock_DEM_data.close()

    def test_get_a_dictionary_from_create_config_dictionary_yml_file(self):
        config_dict =sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        self.assertIsInstance(config_dict, type(dict()))

    def test_get_a_dictionary_of_folders_from_created_yml_dict(self):
        config_dict = sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        folders_dict = config_dict['folders']
        self.assertIsInstance(type(folders_dict), type(dict))

    def test_get_dem_files_folder_from_yml_config(self):
        config_dict = sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        self.assertEqual(config_dict['folders']['DEM_files'], "/data/DEM_2m_HKI")

    def assert_equal_two_rasters(self, test_raster, true_raster):
        self.assertEqual(len(test_raster), len(true_raster))
        length_of_first = len(true_raster[0])
        
        for row in range(len(test_raster)):
            self.assertEqual(len(test_raster[row]), length_of_first)
            self.assertEqual(len(true_raster[row]), length_of_first)
            
            for column in range(len(test_raster[row])):
                self.assertEqual(test_raster[row][column], true_raster[row][column])

    def test_open_dem_file_array(self):
        dem_data_array = sr.read_dem_file_array(test_set_rasters.dem_locations[0])
        self.assert_equal_two_rasters(dem_data_array, self.mock_rasters[0])

    def test_open_dem_file_metadata(self):
        correct_dem_metadata = {'driver': 'GTiff', 'dtype': 'float32', 'nodata': -9999.0,
                                'width': 10, 'height': 21, 'count': 1, 'crs': crs.CRS.from_epsg(3067),
                                'transform': rio.Affine(2.0, 0.0, 362000.0,0.0, -2.0, 6672000.0)}
        dem_metadata = sr.read_dem_file_metadata(test_set_rasters.dem_locations[0])
        self.assertEqual(dem_metadata, correct_dem_metadata)

    def test_add_missing_subfolders(self):
        config = sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        self.assertEqual(sr.add_missing_subfolders(config), 1)

    def test_merge_if_needed_gives_correct_path_True(self):
        config = sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        
        config['merging']['DEM_merge_boolean'] = True
        config['merging']['DEM_file_path_to_merge'] = os.getcwd()
        config['folders']['temporary_files'] = os.getcwd()
        config['folders']['DEM_files'] = os.getcwd()
        config['merging']['DEM_merge_search_criteria'] = "mock_*"
        
        true_True_path = os.path.join(sr.get_path_for('temporary', config), 'DEM_merged.tif')
        
        output_True_path = sr.merge_if_needed(config)
        self.assertEqual(output_True_path, true_True_path)

    def test_merge_if_needed_gives_correct_path_False(self):
        config = sr.create_config_dictionary_from_yml_file(self.config_filename_yml)
        config['merging']['DEM_merge_boolean'] = False
        config['cropping']['DEM_file_path_to_crop'] = "hello"
        output_False_path = sr.merge_if_needed(config)
        self.assertEqual(output_False_path, "hello")



if __name__ == '__main__':
    unittest.main()