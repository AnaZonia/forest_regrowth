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

files <- list.files(path = './data/secondary_forest_age_v2_2018', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
secondary_ages <- lapply(files, rast)

biomass1 = rast("./data/santoro/N00W050_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif")
biomass2 = rast("./data/santoro/N00W060_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif")
biomass3 = rast("./data/santoro/N00W070_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif")

biomass <- merge(biomass1, biomass2, biomass3)

#biomass = rast("./data/santoro/S10W050_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif")

tst <- crop(biomass, secondary_ages[[2]])
ages <- crop(secondary_ages[[2]], tst)

plot(ages)

new_raster <- resample(biomass, secondary_ages[[2]], 'near')

ages[ages==0] = NA
global(ages, 'isNA')
ncell(ages)


extrapolation <- c(secondary_ages[[2]], biomass)
names(extrapolation) <- c('age', 'agbd')

means <- aggregate(agbd ~ age, extrapolation, median)
colnames(means) <- c('age', 'agbd')
means

plot(means$age, means$agbd)

sd <- aggregate(data$agbd, by = list(data$age), FUN = sd)
colnames(sd) <- c('age', 'sd')
data <- merge(data, means, by = 'age')
data <- merge(data, sd, by = 'age')
data[abs(data$agbd-data$mean) < 0.25*data$sd, ]
