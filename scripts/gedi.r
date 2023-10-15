# GEDI on gee


packages <- c("rgee", "reticulate", "tidyverse")
lapply(packages, require, character.only=TRUE) # returns true if the package is loaded

rgee_environment_dir <- "/usr/bin/python3"

reticulate::use_python(rgee_environment_dir, required=T)
rgee::ee_install_set_pyenv(
py_path = rgee_environment_dir,
py_env = "Python3"
)

Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

rgee::ee_Initialize(drive = T)

gedi <- ee$FeatureCollection('LARSE/GEDI/GEDI04_A_002/GEDI04_A_2022157233128_O19728_03_T11129_02_003_01_V002')

ee_print(gedi)
Map$setCenter(lon = -94.77616, lat = 38.9587, zoom = 14)
Map$addLayer(eeObject = gedi, name = "gedi")

# feature collection showing the footprints

