all: run_raster run_grass run_itzi run_export clean_infil

generate: 
	@mkdir GRASS_itzi GRASS_itzi/temp GRASS_itzi/temp/rain GRASS_itzi/itzi_files GRASS_itzi/grassdata GRASS_itzi/grassdata/rain GRASS_itzi/grassdata/infiltration
	
run_raster:
	@echo "\n--> Setting DEM, landcover (friction, losses, infiltration) and rain rasters <--\n"
	@python3 -W ignore set_rasters2.py GRASS_itzi_parameters_helsinki.yml
	
run_rain:
	@echo "\n--> Setting rain rasters <--\n"
	@python3 -W ignore set_rain_rasters.py GRASS_itzi_parameters.ini
	
run_grass:
	@echo "\n--> Setting GRASS environment <--\n"
	@python3 -W ignore set_GRASS2.py GRASS_itzi_parameters_helsinki.yml
	
run_itzi:
	@echo "\n--> Running itzi... <--\n"
	@cd ../../../tmp/gip/itzi
	@itzi run GRASS_itzi/itzi_config_file.ini
	@cd ../../../usr/src/app
	
run_export:
	@echo "\n--> Exporting itzi output files... <--\n"
	@python3 -W ignore export_itzi_output2.py GRASS_itzi_parameters.ini

clean_temp: 
	@rm -r GRASS_itzi/temp
	
clean_infil: 
	@rmdir GRASS_itzi/grassdata/infiltration
	
clean_all:
	@rm -r GRASS_itzi
