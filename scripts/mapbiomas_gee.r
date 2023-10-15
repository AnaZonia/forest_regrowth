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

ages <- ee$Image("users/celsohlsj/public/secondary_vegetation_age_collection71_v5")

ages_2020 <- ee$Image("users/celsohlsj/public/secondary_vegetation_age_collection71_v5")$select('classification_2020')
# get/check the metadata
ee_print(ages_2020)

# centering the plot using Map$setCenter
Map$setCenter(lon = -48, lat = -16, zoom = 4)
# for visualization we use Map$addLayer
vis_ages = list(
  min = 1,
  max = 30)
Map$addLayer(eeObject = ages_2020, visParams = vis_ages, name = "ages_2020")


# Image: MapBiomas land cover map ---------------------------------

# mapbiomas dataset
land_cover = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")
ee_print(land_cover)
# take a look at the different bands

# load year 2020 land cover classification 
land_cover_2020 = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")$
  select("classification_2020")
ee_print(land_cover_2020)

# color palette
# official: https://mapbiomas-br-site.s3.amazonaws.com/downloads/Colecction%206/Cod_Class_legenda_Col6_MapBiomas_BR.pdf
idx_mapbiomas = c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 13, 14, 15, 18, 19, 39, 20, 40, 41,  36,  46,  47,  48,  9, 
                  21,  22,  23,  24,  30,  25,  26,  33,  31,  27) 
hex_mapbiomas = c("#129912", "#006400", "#00ff00", "#687537", "#6b9932", "#BBFCAC", "#45C2A5", "#B8AF4F", "#968c46",
                  "#665a3a", "#f1c232", "#FFFFB2", "#FFD966", "#E974ED", "#D5A6BD", "#e075ad", "#C27BA0", "#982c9e", "#e787f8", "#f3b4f1",
                  "#cca0d4", "#d082de", "#cd49e4", "#ad4413", "#fff3bf", "#EA9999", "#DD7E6B", "#aa0000", "#af2a2a", "#ff3d3d", "#0000FF",
                  "#0000FF", "#02106f", "#D5D5E5")
map_biomas_palette = rep("#FFFFFF", 49)
map_biomas_palette[idx_mapbiomas] = hex_mapbiomas
vis_mapbiomas = list(
  min = 1,
  max = 49,
  palette = map_biomas_palette
)
Map$addLayer(eeObject = land_cover_2020, visParams = vis_mapbiomas, 'mapbiomas_2020')

## You can then plot your own shapefiles on top of this view to visualize land cover

# now lets get the 1985 map
land_cover_1985 = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")$
  select("classification_1985")
ee_print(land_cover_1985)

# and plot both 2010 and 2020 to see the changes
Map$addLayer(eeObject = land_cover_1985, visParams = vis_mapbiomas, 'mapbiomas_1985') +
  Map$addLayer(eeObject = land_cover_2020, visParams = vis_mapbiomas, 'mapbiomas_2020')

