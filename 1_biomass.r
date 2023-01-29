####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Downloading and processing non-GEDI biomass data.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# imports and makes dataframes from santoro and potapov raw data
####################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("grimbough/rhdf5")
library(rhdf5) # for handling raw GEDI data
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R) # for extracting raw GEDI data
library(terra)
library(ncdf4)
library(tidyverse)
library(rgdal)

setwd("/home/aavila/forest_regrowth")
regrowth_mask <- rast('0000000000-0000095232_mask.tif')
coords <- c(-0.5377764, -3.2823093, -48.32644, -43.99998) #ext(regrowth_mask)

# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
# more instructions on https://github.com/VangiElia/GEDI4R
# a package built specially for processing GEDI4A biomass data.

####################################################################
##########              BODY                 #######################
####################################################################


outdir <- tempdir()
#173 output files
GEDI_download <- l4_download(
coords[1], coords[2], coords[3], coords[4], # ul_lat,lr_lat,ul_lon,lr_lon
outdir = outdir,
from = "2020-01-01",
to = "2020-12-31",
just_path = F)

#jpeg("rplot.jpg", width = 350, height = "350")
GEDI_download <- paste0('./GEDI_raw/', list.files('./GEDI_raw/', pattern = '.h5'))
l4 <- l4_getmulti(GEDI_download,just_colnames = F)
l4 <- subset(l4, l4_quality_flag == 1)
l4 <- subset(l4, degrade_flag == 0)

#coords <- c(max(age$lat), min(age$lat), min(age$lon), max(age$lon)) #coordinates for a test region in the northeastern Amazon
clipped <- l4_clip(l4,c(coords[3], coords[2], coords[4], coords[1]))
clipped <- subset(clipped, agbd < 600)
jpeg("rplot_mat.jpg")
l4_plotagb(clipped,n=100,h=c(100,100))
dev.off()





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  SANTORO ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Santoro et al 2018 data -> GlobBiomass ESA (Mg/ha)
# 100m resolution, 2010

biomass1 = raster("./santoro/N00W060_agb.tif")
biomass2 = raster("./santoro/N00W100_agb.tif")
biomass3 = raster("./santoro/N40W060_agb.tif")

biomass = merge(biomass1, biomass2, biomass3)

## crop and mask
r2 <- crop(biomass, extent(BRA))
r3 <- mask(r2, BRA) #save this somewhere
e <- extent(-48.32644, -43.99998, -3.2823093, -0.5377764) # for testing purposes, the coordinates of region 0000000000.0000095232
r4 <- crop(biomass1, e)

bm_test <- getValues(r4)
bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
bm_test <- na.omit(bm_test)
bm_test[,c("x","y")] <- xyFromCell(r4, bm_test$cell)

#biomass = biomass[ymin < biomass$y & biomass$y < ymax & xmin < biomass$x & biomass$x < xmax,]

biomass = cbind(biomass, LongLatToUTM(biomass$x, biomass$y))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  POTAPOV ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Potapov et al 2020 -> GLAD Forest Canopy Height (m)
# 30m resolution, 2019

biomass = raster("Forest_height_2019_SAM.tif")

## crop and mask
r2 <- crop(biomass, extent(BRA))
r3 <- mask(r2, BRA)

#writeRaster(r3, "Forest_height_2019_Brazil.tif")

e <- extent(xmin, xmax, ymin, ymax)

biomass = raster("Forest_height_2019_Brazil.tif")
biomass_cropped <- crop(biomass,e)


#biomass_df <- as.data.frame(biomass_cropped, xy=T, na.rm = TRUE)
bm_test <- values(biomass_cropped)

bm_tst_complete <- na.omit(bm_test)

bm_test <- getValues(biomass_cropped)
bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
bm_test <- na.omit(bm_test)
bm_test[,c("x","y")] <- xyFromCell(biomass_cropped, bm_test$cell)

biomass = cbind(biomass_with_data, LongLatToUTM(biomass_with_data$x, biomass_with_data$y))