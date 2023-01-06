####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Downloading and processing GEDI4A data.
# Ana Avila - Jan 2023
####################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("grimbough/rhdf5")
library(rhdf5) # for handling raw GEDI data
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R) # for extracting raw GEDI data
library(raster) # Might not need this one
library(ncdf4)
library(tidyverse)
library(rgdal)
library(hdf5r)

setwd("/home/aavila/forest_regrowth")
coords <- c(-0.5377764, -3.2823093, -48.32644, -43.99998) #coordinates for a test region in the northeastern Amazon

# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
# NASA's h5 files are downloaded in Polaris under "/home/aavila/forest_regrowth/GEDI_raw".
#If you'd like to download them again, it is fast.

####################################################################
##########              ATTEMPT 1                 ################## 
####################################################################
# downloads GEDI4A product from cloud using GEDI4R package
# converts hdf5 into netcdf files for easier handling (thanks to Dat Nguyen for help with netcdf)
# makes biomass dataframe from raw product with lat, lon, UTM coords, and AGBD (aboveground biomass)
# more instructions on https://github.com/VangiElia/GEDI4R   ->   a package built specially for processing GEDI4A biomass data.

############ GEDI 4A documentation states that: ##############
# https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density_V2_1.html#acqmatmethods
# The GEDI instrument consists of three lasers producing a total of eight beam ground transects,
# which instantaneously sample eight ~25 m footprints spaced approximately every 60 m along-track.
# The GEDI beam transects are spaced approximately 600 m apart on the Earth's surface in the
# cross-track direction, for an across-track width of ~4.2 km.
#       - So, each beam gives values for a few footprints.

outdir <- ("./GEDI_raw/")
GEDI_download = l4_download(
  coords[1], coords[2], coords[3], coords[4], # ul_lat,lr_lat,ul_lon,lr_lon
  outdir = outdir,
  from = "2020-01-01",
  to = "2020-12-31",
  just_path = F )

nc_data <- c(paste0("./GEDI_raw/", list.files("./GEDI_raw/", pattern = '.h5')))
test_convert <- function(nc_file) tryCatch(nc_open(nc_file), error = function(e) e)
nc_data2 <- lapply(nc_data, test_convert)

df_from_nc <- function(nc_file){
    # Get variable names
    #nc_file <- crop(nc_file, shape)
    nmv = names(nc_file$var)
    # Which ones are agbd
    nmv.values = nmv[grepl("/agbd$", nmv)]
    # Which ones are lat/lon
    nmv.lat = nmv[grepl("/lat_lowestmode$", nmv)]
    nmv.lon = nmv[grepl("/lon_lowestmode$", nmv)]

    # Get data for the first one, might want to iterate over 1:length(nmv.values), or beams or whatever
    # -9999 seems to be their way of doing NA
    df1 = data.frame(agbd = ncvar_get(nc_file, nmv.values[1]),
        lat = ncvar_get(nc_file, nmv.lat[1]),
        lon = ncvar_get(nc_file, nmv.lon[1]))
    df1 <- subset(df1, agbd > 0)

    for (i in 2:8){   # because there's eight beams, we merge them all together.
        df2 = data.frame(agbd = ncvar_get(nc_file, nmv.values[i]),
            lat = ncvar_get(nc_file, nmv.lat[i]),
            lon = ncvar_get(nc_file, nmv.lon[i]))
        df2 <- subset(df2, agbd > 0)
        df1 <- rbind(df1, df2)
    }

    df3 <- subset(df1, lat < coords[1] & lat > coords[2] & lon > coords[3] & lon < coords[4])
    head(df3)
}

df_list <- lapply(nc_data2, df_from_nc)
df_unified <- bind_rows(df_list)

hist(df_unified$agbd,breaks <- 5000)

#looks like a strange negative-exponential. Absolutely no clue as to why.
# NOTE: when looking at nmv, there are many variables stored in this file.
# maybe there is another action that needs to be done with one of the other variables to get the right value?

####################################################################
##########              ATTEMPT 2                 ################## 
####################################################################

# Using GEDI4R package: (guideline in https://github.com/VangiElia/GEDI4R)

dataname <- l4_get(nc_data[[1]],just_colnames = T)  # the following error is returned:
# gediL4_path has lenght==1,l4_get will be used in single thread mode
# Error in H5File.open(filename, mode, file_create_pl, file_access_pl) : 
#   HDF5-API Errors:
#     error #000: ../../../src/H5F.c in H5Fopen(): line 509: unable to open file
#         class: HDF5
#         major: File accessibilty
#         minor: Unable to open file

#     error #001: ../../../src/H5Fint.c in H5F__open(): line 1400: unable to open file
#         class: HDF5
#         major: File accessibilty
#         minor: Unable to open file

#     error #002: ../../../src/H5Fint.c in H5F_open(): line 1738: file close degree doesn't match
#         class: HDF5
#         major: File accessibilty
#         minor: Unable to initialize object

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Looking into the function l4_getmulti, I find that it uses the following snippet to read code:

level4a_h5 <- hdf5r::H5File$new(nc_data[[1]], mode = "r") #mode="r" to open in reading mode
# which is where the mentioned error happens. Somehow, hdf5r can't open this file.

# However, when I use an alternative package, this command runs just fine:
level4a_h5 <- rhdf5::h5ls(nc_data[[1]])
# the file is there, and it's readable. However,
level4a_h5 <- rhdf5::H5Fopen(nc_data[[1]])
# Error in rhdf5::H5Fopen(nc_data[[1]]) : 
#   HDF5. File accessibility. Unable to open file.

# This, however, works fine:
level4a_h5 <- rhdf5::h5read(nc_data[[1]], "BEAM1011/agbd")

# Bit 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#To make sure that l4_download isn't the issue,
#(the files aren't corrupted from the download function):

man_dwnld <- "/home/aavila/forest_regrowth/GEDI_manual_dwnld/GEDI04_A_2020019064709_O06246_04_T00785_02_002_02_V002.h5"

gediL4_path <- l4_get(man_dwnld)

hist(gediL4_path$agbd,breaks <- 5000)

level4a_h5 <- hdf5r::H5File$new(man_dwnld, mode = "r") #mode="r" to open in reading mode

# However, when I use an alternative package, this command runs just fine:
level4a_h5 <- rhdf5::h5ls(man_dwnld)
# the file is there, and it's readable. However,
level4a_h5 <- rhdf5::H5Fopen(nc_data[[1]])
# Error in rhdf5::H5Fopen(nc_data[[1]]) : 
#   HDF5. File accessibility. Unable to open file.

# This, however, works fine:
level4a_h5 <- rhdf5::h5read(nc_data[[1]], "BEAM1011/agbd")
