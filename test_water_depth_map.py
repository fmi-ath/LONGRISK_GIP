import unittest 
import rasterio as rio
import modules.water_depth_map as wdm
import numpy as np
import cv2 as cv

class test_water_depth_map(unittest.TestCase):
    def test_get_single_dem_gives_same_file(self):
        test_file_path = __file__[:(__file__.rfind('/')+1)]+"../data/DEM_2m_HKI/L4131C.tif"
        
        with rio.open(test_file_path, 'r') as src:
            profile_true = src.profile
            data_true = src.read(1)

        data_test, profile_test = wdm.single_dem(test_file_path)
        for key in profile_test.keys():
            self.assertEquals(profile_test[key], profile_true[key])
        # MOCKAA
        #for row in range(len(data_test)):
        #    for col in range(len(data_test[0])):
        #        self.assertEquals(data_test[row][col],data_true[row][col])

    def compare_dem(self, dem1, dem2):
        for row in range(len(dem1)): 
            for col in range(len(dem1[0])):
                self.assertEqual(dem1[row,col], dem2[row,col])

    def test_FloodFill_1x1_watermap(self):
        dem = np.array([[0]])
        watermap = wdm.generate_water_depth_map_using_floodfill(dem, (0,0), addon=0)
        self.assertEquals(watermap, np.array([[0]]))

    def test_FloodFill_2x2_watermap(self):
        dem = np.array([[1,0],[1,0]])
        result = np.array([[0,1],[0,1]])
        watermap = wdm.generate_water_depth_map_using_floodfill(dem, (0,0), addon=0)
        self.compare_dem(result, watermap)
        
    def test_FloodFill_3x4_watermap(self):
        dem = np.array([[1,2,3,4],[4,3,2,1],[1,2,3,4]])
        result = np.array([[2,1,0,0],[0,0,1,2],[2,1,0,0]])
        watermap = wdm.generate_water_depth_map_using_floodfill(dem, (1,1), addon=0)
        self.compare_dem(result, watermap)
    
    def test_FloodFill_5x5_watermap(self):
        dem = np.array([[0.001, 0.001, 0.001, 0.001, 0.001],[0.003, 0.003, 0.003, 0.003, 0.003],[0.0, 0.001, 0.002, 0.001, 0.0],[0.003, 0.003, 0.003, 0.003, 0.003],[0.001, 0.001, 0.001, 0.001, 0.001]])
        result = np.array([[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0.002, 0.001, 0, 0.001, 0.002],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0]])
        watermap = wdm.generate_water_depth_map_using_floodfill(dem, (2,2), addon=0)
        self.compare_dem(result, watermap)

    def test_diagonal_3x3(self):
        dem = np.array([[0,2,0],[2,1,2],[0,2,0]])
        result = np.array([[1,0,1],[0,0,0],[1,0,1]])
        watermap = wdm.generate_water_depth_map_using_floodfill(dem, (1,1), addon=0)
        print(watermap)
        self.compare_dem(result, watermap)



if __name__ == '__main__':
    unittest.main()
        