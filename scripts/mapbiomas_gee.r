####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Ana Avila - Dec 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

packages <- c("rgee", "reticulate", "tidyverse", "sf", "mapview")
lapply(packages, require, character.only=TRUE) # returns true if the package is loaded

rgee_environment_dir <- "/usr/bin/python3"

reticulate::use_python(rgee_environment_dir, required=T)
rgee::ee_install_set_pyenv(
py_path = rgee_environment_dir,
py_env = "Python3"
)
 
Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

rgee::ee_Initialize() #use drive = TRUE to access things in your google drive

####################################################################
# -------------------------   Switches    ------------------------ #
# land_use = FALSE
# climate = FALSE
# mature = FALSE

####################################################################

####################################################################
# -------------------   Age and Biomass    ----------------------- #
####################################################################

# Load images and feature collection
age <- ee$Image('users/celsohlsj/public/secondary_vegetation_age_collection71_v5')$
                select('classification_2020')
biomass <- ee$Image('projects/ee-ana-zonia/assets/biomass_2020')
sd <- ee$Image('projects/ee-ana-zonia/assets/biomass_sd_2020')
cwd <- ee$Image('projects/ee-ana-zonia/assets/cwd_chave')
amazon_biome <- ee$FeatureCollection('projects/ee-ana-zonia/assets/amazon_biome_border')
ecoregions <- ee$FeatureCollection('RESOLVE/ECOREGIONS/2017')

# Clip images to Amazon biome
age <- age$clip(amazon_biome$geometry())$updateMask(age$gt(0)) # keep only pixels with ages > 0 (secondary forests)
biomass <- biomass$clip(amazon_biome$geometry())
sd <- sd$clip(amazon_biome$geometry())
cwd <- cwd$clip(amazon_biome$geometry())
ecoregions <- ecoregions$filterBounds(amazon_biome$geometry())
ecoregions_list <- ecoregions$toList(ecoregions$size())
ee_print(ecoregions_list$getInfo())

# Reproject images to 10m
biomass_10m <- biomass$reproject(crs = age$projection(), scale = 10)
cwd_10m <- cwd$reproject(crs = age$projection(), scale = 10)
sd_10m <- sd$reproject(crs = age$projection(), scale = 10)

# Reaggregate to 30m (mean value)
aggregated_biomass <- biomass_10m$reduceResolution(reducer = ee$Reducer$mean())$reproject(crs = age$projection())
aggregated_cwd <- biomass_10m$reduceResolution(reducer = ee$Reducer$mean())$reproject(crs = age$projection())
aggregated_sd <- biomass_10m$reduceResolution(reducer = ee$Reducer$mean())$reproject(crs = age$projection())

# Mask only to regions with age greater than zero (secondary forests)
aggregated_biomass <- aggregated_biomass$updateMask(age)$toInt16()
aggregated_sd <- aggregated_sd$updateMask(age)$toInt16()
aggregated_cwd <- aggregated_cwd$updateMask(age)$toInt16()

ecoreg <- ecoregions_list.get(2)
print(ecoreg)
ee_print(ecoreg$getInfo('ECO_NAME'))

# export_data <- function(img, prefix) {
#   # Export data for each ecoregion
#   for (i in 1:31) {
#     ecoreg <- ecoregions_list$get(i)
#     img_clipped <- img$clip(ecoreg)
#     projection <- img_clipped$projection()$getInfo()

#     ee$batch$Export$image$toDrive(
#     image = img_clipped,
#     description = paste0(prefix, '_', ecoreg$get('ECO_NAME')),
#     crs = img$projection()$crs(),
#     crsTransform = img$projection()$getInfo()$transform,
#     region = ecoregion$geometry(),
#     maxPixels = 4e10
#     )
#   }
# }

export_data(age, 'age')

####################################################################
# -----------------------   Land Use    -------------------------- #
####################################################################

# # mapbiomas dataset
# land_cover = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")$select('classification_1985')
# ee_print(land_cover)
# # take a look at the different bands
# land_cover_2020 = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")$
#   select("classification_2020")

# # color palette
# # official: https://mapbiomas-br-site.s3.amazonaws.com/downloads/Colecction%206/Cod_Class_legenda_Col6_MapBiomas_BR.pdf
# idx_mapbiomas = c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 13, 14, 15, 18, 19, 39, 20, 40, 41,  36,  46,  47,  48,  9, 
#                   21,  22,  23,  24,  30,  25,  26,  33,  31,  27) 
# hex_mapbiomas = c("#129912", "#006400", "#00ff00", "#687537", "#6b9932", "#BBFCAC", "#45C2A5", "#B8AF4F", "#968c46",
#                   "#665a3a", "#f1c232", "#FFFFB2", "#FFD966", "#E974ED", "#D5A6BD", "#e075ad", "#C27BA0", "#982c9e", "#e787f8", "#f3b4f1",
#                   "#cca0d4", "#d082de", "#cd49e4", "#ad4413", "#fff3bf", "#EA9999", "#DD7E6B", "#aa0000", "#af2a2a", "#ff3d3d", "#0000FF",
#                   "#0000FF", "#02106f", "#D5D5E5")
# map_biomas_palette = rep("#FFFFFF", 49)
# map_biomas_palette[idx_mapbiomas] = hex_mapbiomas
# vis_mapbiomas = list(
#   min = 1,
#   max = 49,
#   palette = map_biomas_palette
# )

# Map$addLayer(eeObject = land_cover_2020, visParams = vis_mapbiomas, 'mapbiomas_2020')

# ## You can then plot your own shapefiles on top of this view to visualize land cover

# # now lets get the 1985 map
# land_cover_1985 = ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")$
#   select("classification_1985")
# ee_print(land_cover_1985)

# # and plot both 2010 and 2020 to see the changes
# Map$addLayer(eeObject = land_cover_1985, visParams = vis_mapbiomas, 'mapbiomas_1985') +
#   Map$addLayer(eeObject = land_cover_2020, visParams = vis_mapbiomas, 'mapbiomas_2020')

# asset = ee$Image("projects/ee-ana-zonia/assets/biomass_mask_2018")
# ee_print(asset)
# Map$addLayer(asset, visParams =)

####################################################################
# ------------------------    Fire    ---------------------------- #
####################################################################


# num_fires <- sum(fire_brick_masked)

# # find when was last fire observed (how many years before observed regrowth)
# fire_brick_flipped <- fire_brick_masked [[ c(rev(order(names(fire_brick_masked)))) ]]
# last_fire <- which.lyr(fire_brick_flipped == 1)


####################################################################
# ---------------------   Mature Forest    ----------------------- #
####################################################################

# import mapbiomas forest cover script
# select for the pixels between 200 and 300
# mature[mature > 300 | mature < 200] <- NA # only leave behind mature values
# mature_mask <- mask(mature, app(mature, fun = sum))
# consider only pixels that have been mature the whole time

# get biomass of surrounding mature forests
# with adjustable buffer size



####################################################################
# -------------------------   Climate    ------------------------- #
####################################################################

# cwd <- ee$Image('projects/ee-ana-zonia/assets/cwd_brazilian_amazon_chave')
# terraclimate <- ee$Image("IDAHO_EPSCOR/TERRACLIMATE")
