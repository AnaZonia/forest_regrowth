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
age <- ee$Image('users/celsohlsj/public/secondary_vegetation_age_collection71_v5')$select('classification_2020')
biomass <- ee$Image('projects/ee-ana-zonia/assets/biomass_2020')
sd <- ee$Image('projects/ee-ana-zonia/assets/biomass_sd_2020')
cwd <- ee$Image('projects/ee-ana-zonia/assets/cwd_chave')
amazon_biome <- ee$FeatureCollection('projects/ee-ana-zonia/assets/amazon_biome_border')

# Define visualization parameters
age_viz <- list(min = 0, max = 33, palette = c('00FFFF', '0000FF'))
agbd_viz <- list(min = 0, max = 415, palette = c('00FFFF', '0000FF'))
cwd_viz <- list(min = -831.75, max = 0, palette = c('00FFFF', '0000FF'))

Map$addLayer(age, age_viz, 'age')
Map$addLayer(biomass, agbd_viz, 'agbd')
Map$addLayer(cwd, cwd_viz, 'cwd')

# Clip images to Amazon biome
age <- age$clip(amazon_biome$geometry())$updateMask(age$gt(0))$rename('age') # keep only pixels with ages > 0 (secondary forests)
biomass <- biomass$clip(amazon_biome$geometry())
sd <- sd$clip(amazon_biome$geometry())
cwd <- cwd$clip(amazon_biome$geometry())

# Reproject images to 10m
ages_10m <- age$reproject(crs = age$projection()$crs(), scale = 10)
biomass_10m <- biomass$reproject(crs = age$projection(), scale = 10)
cwd_10m <- cwd$reproject(crs = age$projection(), scale = 10)
sd_10m <- sd$reproject(crs = age$projection(), scale = 10)

# Mask everything to only encompass secondary forest patches
biomass_10m <- biomass_10m$updateMask(ages_10m)
ages_10m <- ages_10m$updateMask(biomass_10m)
sd_10m <- sd_10m$updateMask(ages_10m)
cwd_10m <- cwd_10m$updateMask(ages_10m)

# Reaggregate to 30m (mean value)
aggregated_biomass <- biomass_10m$reduceResolution(reducer = ee$Reducer$mean(), maxPixels = 1024)$reproject(crs = age$projection()$crs(), scale = 30)
aggregated_cwd <- cwd_10m$reduceResolution(reducer = ee$Reducer$mean(), maxPixels = 1024)$reproject(crs = age$projection()$crs(), scale = 30)
aggregated_sd <- sd_10m$reduceResolution(reducer = ee$Reducer$mean(), maxPixels = 1024)$reproject(crs = age$projection()$crs(), scale = 30)


Map$addLayer(aggregated_cwd, cwd_viz, 'cwd')
Map$addLayer(aggregated_biomass, age_viz, 'ages_2020')

Map$addLayer(amazon_imgcol$select('cwd'), cwd_viz, 'cwd')
Map$addLayer(amazon_imgcol$select('age'), age_viz, 'ages_2020')
Map$addLayer(amazon_imgcol$select('agbd'), agbd_viz, 'agbd')


# Create latitude and longitude bands
# proj <- ages_10m$projection()
# latlon <- ee$Image$pixelLonLat()$reproject(proj)
# latlon <- latlon$updateMask(ages_10m)
# lat <- latlon$select("latitude")
# lon <- latlon$select("longitude")

# Create secondary_amazon image
secondary_amazon <- ee$Image(ages_10m)$addBands(lat)$addBands(lon)$addBands(biomass_10m$rename('biomass'))$addBands(sd_10m$rename('sd'))$addBands(cwd_10m$rename('cwd'))

# Split image for exporting
# Construct grid and intersect with country polygon
grid <- ages_10m$geometry()$coveringGrid(ages_10m$geometry()$projection(), 500000)

for (f in 1:grid$size()$getInfo()) {
  print(f)
  second_feature <- ee$Feature(grid$toList(grid$size()$getInfo())$get(f))
  second_tiles <- secondary_amazon$clip(second_feature)
  sampl_clip <- second_tiles$stratifiedSample(100)
  # export all features in a FeatureCollection as one file
  task <- ee$batch$Export$table(sampl_clip, paste0('sampl_clip', f), list(fileFormat = 'CSV'))
  task$start()
  ee_monitoring()
}

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
