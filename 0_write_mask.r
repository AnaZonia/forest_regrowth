####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Writing raster masks selecting for only secondary forest patches at 2019.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAPBIOMAS raw data comes divided in 12 different regions of the Amazon, each with their own identifier (a location number).
# This script builds rasters from raw materials containing only regrowth events.
# This allows to subset the data to contain only pixels that show regrowth history
# The output masks are input into 1_processing_pred_var.R
####################################################################


# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)
path <- './mapbiomas/regrowth_raw'
files <- list.files(path)
#newname <- sub('none-','', files) ## making names easier to read, standardizing the names
#file.rename(file.path(path,files), file.path(path, newname)) ## renaming it.
locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)
#locations <- locations[1:length(locations)]

# create a raster stack with one layer per year, for this location.
# this stack will be merged to make the regrowth-only mask.

# here we go location by location:
  # import raw data
  # make a mask with only areas that show regrowth - this mask will be used to reduce the size of files we're handling
  # save these rasters in the ./regrowth_rasters directory

#for (i in 1:length(locations)){}
# since each mask takes 6h and the ssh connection tends to break after too long of inactivity, we start with an example location:
  location <- '0000000000-0000095232'

  files_tmp <- list.files(path = './mapbiomas/regrowth_raw', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
  files_tmp <- sort(files_tmp) #making sure it's ordered by year

  regrowth_list <- c()
  for (i in 1:length(files_tmp)){
    regrowth_list <- c(regrowth_list, raster(files_tmp[i]))
    print(i)
  }

  # obtain the raster file for all years within that location.
  # to make processing lighter, subset only the pixels that have shown regrowth history.
  # here, we are (1) making a mask registering all regrowth moments and (2) subsetting rasters based on that mask.
  # a regrowth moment is flagged with the value "503", therefore:
  for (i in 1:32){
    print(i)
    print(Sys.time())
    regrowth_list[[i]][regrowth_list[[i]]!=503] <- NA # only leave behind values not including the regrowth moment
    writeRaster(regrowth_list[[i]], file.path(paste0('./mapbiomas/regrowth_rasters/', location, '_', c(1987+i), ".tif"))) # save rasters with the year on the folder created for each location.
  }

  # each mask takes about 6 hours to be processed.
  # these files are being stored in the machine to avoid losing all progress in case of a disconnection of the ssh,
  # as well as avoid having to redo everything in case there is an error that could be fixed later by editing the files.

  filepaths <- paste0("/home/aavila/Documents/forest_regrowth/mapbiomas/regrowth_rasters/", location, '_', c(1988:2019), '.tif')
  mask_raster_list <- lapply(filepaths, raster)
  mask_raster_stack <- stack(mask_raster_list)
  regrowth_mask <- merge(mask_raster_stack)
  regrowth_mask <- writeRaster(regrowth_mask, file.path(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))) # mask made and saved.

