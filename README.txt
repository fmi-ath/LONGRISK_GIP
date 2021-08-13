Author: Miguel Aldana
Date: July 30th 2021


This is module is meant for running itzï simulations in GRASS without opening GRASS explicitly through Python.

For running the simulation, the user has to fill in the file 'GRASS_itzi_parameters.ini', according to the files available
for the simulation and parameters needed for it.

The following packages have to be installed in order to perform correctly:

- Grass
- itzi
- rasterio
- osgeo
- shapely
- geopandas
- fiona
- pycrs
- pyproj

Apart from filling the file 'GRASS_itzi_parameters.ini', the user does not have to specify anything else, just download the neccesary DEM files and indicate the corresponding folder of their locations. The module runs the simulation over a region according to shapefile or polygon specified by the user in the same file.

After filling the file, run the command in the terminal:

  $ make all
  
After the process is completed, in the current folder the user will find a folder called 'GRASS_itzi', which includes the final rasters used in 'GRASS_itzi/grassdata' and itzi outputs in 'GRASS_itzi/itzi_files'. 
