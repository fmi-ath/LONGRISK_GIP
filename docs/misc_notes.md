---
title: Miscellaneous LONGRISK_GIP notes
author: Petteri Karsisto
year: 2021
---

# Miscellaneous notes

A collection of miscellaneous notes.

## Configuration files

- `GRASS_itzi_parameters.ini`: This file includes all configurable parameters. The settings are
    described in the file, there is no companion documentation for it. The sections starting with
    `grass_` are used in generating the Itzi configuration file.
    - Note that the makefile is a bit outdated, it assumes that the configuration file is in the
        same folder as the makefile, and it's named `GRASS_itzi_parameters.ini`. Run the Python
        code manually if you want to specify which configuration file to use.
- `itzi_config_file.ini`: This file contains the Itzi configuration. It is generated from the above
    configuration file and will be located in `itzi_output_files` directory (as defined in
    `GRASS_itzi_parameters.ini`). See Itzi documentation for parameter explanations.
     (see also `modules.GRASS_utils.create_itzi_config_file()`)

## Python

Order of execution (based on `make all`)
1. `set_rasters.py`
    - depends on `modules/Landcover.py` and `modules/utils.py`
2. `set_GRASS.py`
    - depends on `modules/GRASS_utils.py`
3. `export_itzi_output`

`set_GRASS.py` and `modules/GRASS_utils.py` depend on GRASS
