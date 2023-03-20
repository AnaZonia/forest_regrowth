
################################## TESTING ZONE  ######################################################
tst <- readRDS('./test_files/agb_forest_age_santoro_paragominas.rds')

head(agb_forest_age)

biomass[biomass$agbd == 49 & biomass$x == 213210 & biomass$y == 9663510,]

regrowth[rownames(regrowth) == 505845282, ] 
biomass[biomass$xy == 228090109643800, ] 



# lots of low biomasses at all ages... why? Biomass is showing as lower than expected.
# could be misalignment - showing nearby patches of low biomass.
# could be misclassification
#remove all cases of 305
#remove all cases with 300-500-300
#remove all cases of less than 4 consecutive years of 100-class before regrowth (seem to be accounting for high biomass, low age)
# GEDI is showing some suspiciously low AGBD readings.


tst = subset(agb_forest_age, agbd < 50)
tst = subset(tst, agbd > 13)
tst = subset(tst, forest_age == 29)


hist(biomass$agbd, breaks=20)


biomass |>
  subset(abs(N00W060_agb - 67) < 1e-3)


biomass |>
  subset(abs(Lon - -100.7) < 1e-5 & abs(Lat - 59.6) < 1e-5)

Q <- quantile(biomass$agbd, probs=c(.1, .9), na.rm = FALSE)
iqr <- IQR(biomass$agbd)
up <-  Q[2]+1.5*iqr # Upper Range  
eliminated<- subset(biomass, biomass$agbd < (Q[2]+1.5*iqr))

##################

extent_list = lapply(tmp_dfs, extent)
# make a matrix out of it, each column represents a raster, rows the values
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")

best_extent = extent(min(matrix_extent[1,]), max(matrix_extent[3,]),
min(matrix_extent[2,]), max(matrix_extent[4,]))

# the range of your extent in degrees
ranges<-apply(as.matrix(best_extent), 1, diff)
# the resolution of your raster (pick one) or add a desired resolution
reso<-res(tmp_dfs[[1]])
# deviding the range by your desired resolution gives you the number of rows and columns
nrow_ncol<-ranges/reso

# create your raster with the following
s<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=tmp_dfs[[1]]@crs)

##########

results <- list()

for(i in 1:length(tmp_dfs)) {
  print(i)
  e <- extent(s)
  r <-tmp_dfs[[i]] # raster(files[i])
  rc <- crop(tmp_dfs[[i]], e)
  if(sum(as.matrix(extent(rc))!= as.matrix(e)) == 0){ # edited
    rc <- mask(rc, a) # You can't mask with extent, only with a Raster layer, RStack or RBrick
  }else{
    rc <- extend(rc,s)
    rc<- mask(rc, s)
  }

  # commented for reproducible example      
  results[[i]] <- rc # rw <- writeRaster(rc, outfiles[i], overwrite=TRUE)
  # print(outfiles[i])
}



for (i in 1988:2019){
  tmp = raster::stack(list.files(path = './mapbiomas/regrowth_amazon', pattern='2016', full.names=TRUE))
  }




cores <- 50
cl <- makeCluster(cores) #output should make it spit errors
registerDoParallel(cl)

############ Splitting a raster for easier handling

# The function spatially aggregates the original raster
# it turns each aggregated cell into a polygon
# then the extent of each polygon is used to crop
# the original raster.
# The function returns a list with all the pieces
# in case you want to keep them in the memory. 
# it saves and plots each piece
# The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)
SplitRas <- function(raster,ppside){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  return(r_list)
}

split_list = SplitRas(r2, 4)

foreach(i=1:length(split_list)) %dopar% {
  writeRaster(mask(split_list[i], Brazil))
}



############ Creating one unified dataframe from multiple raw Mapbiomas .tif files


# It intakes:
# file = the path to the directory (string)
# crop = whether we are subsetting the data into our region of interest or maintaining it whole (Boolean)
# It outputs:
# Converts into dataframes with the year as column
making_df = function(file, crop){

  files = list.files(path = file, pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 

  tmp_rasters = lapply(files, raster)

  # if we are subsetting the data into our region of interest
  coord_oi = c(xmin, xmax, ymin, ymax) #this specifies the coordinates of Paragominas.
  if(crop == T){
    e = as(extent(coord_oi), 'SpatialPolygons')
    crs(e) = "+proj=longlat +datum=WGS84 +no_defs"
    tmp_rasters = lapply(tmp_rasters, crop, e) # subsects all rasters to area of interest
  }

  tmp_dfs = lapply(tmp_rasters, as.data.frame, xy=T)

  merged_df = df_merge(tmp_dfs)

  colnames(merged_df) = str_sub(colnames(merged_df), start= -4)   # makes column names only "yyyy"

  merged_df = merged_df[order(merged_df$x),]   #order by longitude, so the pixels are separated by UTM zones.

  merged_df = cbind(merged_df, LongLatToUTM(merged_df$x, merged_df$y))   # converts lat, long coordinates to UTM

  return(merged_df)
}


############ DIAGNOSE CORRUPTED FILES

for (i in 1:12){
  location = locations[i]
  # the Amazon is divided in 12 parts, each with their own identifier
  # each location is an identifier.
  #for (location in locations){
    files_tmp <- list.files(path = './mapbiomas/regrowth_amazon', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
    files_tmp <- sort(files_tmp)

    #dir.create(paste0('./regrowth_dataframes/', location))

    # obtain the raster file for all years within that location

    raster2 <- function(file, ...) {
    tryCatch(raster(file, ...), error = function(c) {
      c$message <- paste0(c$message, " (in ", file, ")")
      print(c)
    })
    }

    for (j in 1:length(files_tmp)){
      tmp <- raster2(files_tmp[j])
    }
}

############ DOWNLOAD AND MANIPULATE LANDSAT FILES

library(reticulate)
py_config()

Sys.setenv(RETICULATE_PYTHON = '/usr/local/lib/python2.7/')
#use_python('/usr/local/lib/python2.7/')

repl_python()

# ~~~~ Begin Python environment


import ee

#automate
folder <- 'LE07_L2SP_012054_20120217_20200909_02_T1'
bands <- c('_SR_B1') #etc
paste(folder, bands)

x <- list(r1, r2)
names(x) <- c("x", "y")
x$filename <- 'test.tif'
x$overwrite <- TRUE
m <- do.call(merge, x)


#how to deal if the UTM coordinates span zones?

labels <- read.csv('./data/training.csv' )
training <- getKMLcoordinates('./data/training.kml', ignoreAltitude=TRUE)
names(training) <- labels$Name
training <- training[5:length(training)]

training <- lapply(training, as.data.frame)
for (i in 1:length(training)){
  training[[i]]$label <- names(training[i])
}
training <- do.call(rbind, training)
colnames(training) <- c('x','y', 'label')

utm <- LongLatToUTM(training$x,training$y,17)
training <- cbind(utm, training)

training <- training %>%
  rename(x = 1,
         y = 2,
         long = 3,
         lat = 4)
training$label = as.factor(tolower(training$label))

training <- training[c(5:nrow(training)), -c(5,7,8,9)]
saveRDS(training, 'training.rds')


setwd("/home/aavila/Documents/forest_regrowth")



biomass <- readRDS('df_unified.rds')
regrowth <- readRDS('regrowth_cleaned.rds')


biomass <- cbind(biomass, LongLatToUTM(biomass$lon, biomass$lat))

biomass$xy <- paste0(biomass$zone, biomass$x, biomass$y)

agb_forest_age <- cbind(regrowth, agbd = biomass[match(regrowth$xy,biomass$xy),c("agbd")])

agb_forest_age <- agb_forest_age[complete.cases(agb_forest_age[, ncol(agb_forest_age)]), ]

plot(agb_forest_age$forest_age, agb_forest_age$agbd)


sds <- aggregate(agbd ~ forest_age, agb_forest_age, sd)
means <- aggregate(agbd ~ forest_age, agb_forest_age, mean)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'mean', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = mean)) + 
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd)) +
  geom_point() + theme(text = element_text(size = 20))  

tst = subset(agb_forest_age, forest_age == 28)
tst = subset(tst, agbd < 25)


regrowth[rownames(regrowth) == 473167959, ] 
biomass[biomass$xy == 231981809940500, ] 


##################################################################

setwd("/home/aavila/Documents/forest_regrowth")

fire = readRDS('./mapbiomas/dataframes/0000000000-0000095232_fire.rds')
# > range(fire$x)
# [1] -48.31863 -43.99998
# > range(fire$y)
# [1] -3.2823093 -0.5377764
lulc = readRDS('./test_files/lulc.rds')
# > range(lulc$lat)
# [1] -3.2823093 -0.5377764
# > range(lulc$lon)
# [1] -48.32644 -43.99998
regrowth_mask <- raster('./mapbiomas/regrowth_masks/0000000000-0000095232_mask.tif')
# class      : Extent 
# xmin       : -48.32658 
# xmax       : -43.99984 
# ymin       : -3.282444 
# ymax       : 5.272392 
GEDI = readRDS('./GEDI_dataframes/0000000000-0000095232_GEDI.rds')
# > range(GEDI$lat_lowestmode)
# [1] -3.2823053 -0.5379698
# > range(GEDI$lon_lowestmode)
# [1] -48.32643 -43.99999


check = readRDS('./test_files/regrowth_cleaned.rds')
head(check)
# 22 833940 9940500 

check2 = raster('./test_files/merged_years.tif')
# merged for mask of different year

sds <- aggregate(agbd ~ forest_age, agb_forest_age, sd)
means <- aggregate(agbd ~ forest_age, agb_forest_age, mean)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'mean', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = mean)) + 
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd)) +
  geom_point() + theme(text = element_text(size = 20))  


install.packages('pbapply')

######################

path <- './mapbiomas/regrowth_rasters'
files <- list.files(path)
locations <- str_sub(files, end = -10)
locations <- unique(locations)
locations[3]


list.files(path = './mapbiomas/regrowth_rasters', pattern=locations[6], full.names=TRUE)   # obtain paths for all files for that location

raster::extent(raster('0000000000-0000095232_lulc.tif'))

# > range(regrowth$lat)
# [1] -11.141041  -3.282579
# > range(regrowth$lon)
# [1] -73.86052 -65.43638


# [1] "/home/aavila/Documents/forest_regrowth/mapbiomas/regrowth_rasters/0000000000-0000031744_2006.tif"
# has different extent - something went wrong there.
# [31] "./mapbiomas/regrowth_rasters/0000000000-0000000000_2019.tif" is missing

######################

path <- './mapbiomas/regrowth_rasters'
files <- list.files(path)
locations <- str_sub(files, end = -10)
locations <- unique(locations)
locations[3]


list.files(path = './mapbiomas/regrowth_rasters', pattern=locations[6], full.names=TRUE)   # obtain paths for all files for that location

raster::extent(raster('./mapbiomas/regrowth_rasters/0000031744-0000000000_1988.tif'))

# > range(regrowth$lat)
# [1] -11.141041  -3.282579
# > range(regrowth$lon)
# [1] -73.86052 -65.43638


# [1] "/home/aavila/Documents/forest_regrowth/mapbiomas/regrowth_rasters/0000000000-0000031744_2006.tif"
# has different extent - something went wrong there.
# [31] "./mapbiomas/regrowth_rasters/0000000000-0000000000_2019.tif" is missing

# central_df <- cbind(agb_forest_age, last_burn = fire[match(agb_forest_age$xy,fire$xy),c("last_burn")])
# central_df <- cbind(central_df, num_fires = fire[match(central_df$xy,fire$xy),c("num_fires")])
# central_df <- cbind(central_df, pasture = lulc[match(central_df$xy,lulc$xy),c("pasture")])
# central_df <- cbind(central_df, soy = lulc[match(central_df$xy,lulc$xy),c("soy")])
# central_df <- cbind(central_df, other_perennial = lulc[match(central_df$xy,lulc$xy),c("other_perennial")])
# central_df <- cbind(central_df, other_annual = lulc[match(central_df$xy,lulc$xy),c("other_annual")])
# central_df <- cbind(central_df, tavg = temp[match(central_df$xy,temp$xy),c("mean")])
# central_df <- cbind(central_df, prec = prec[match(central_df$xy,prec$xy),c("mean")])
# central_df$last_burn[is.na(central_df$last_burn)] <- 1
# central_df <- lapply(central_df, as.numeric)


writeRaster(df_prec, "df_prec_BRA.tif")
writeRaster(df_tmin, "df_tmin_BRA.tif")
writeRaster(df_tmax, "df_tmax_BRA.tif")

tmin <- as.list(intersect(list.files(path = outdir, pattern = "*.tif", full.names = TRUE), list.files(path = outdir, pattern = "tmin", full.names = TRUE)))
tmax <- as.list(intersect(list.files(path = outdir, pattern = "*.tif", full.names = TRUE), list.files(path = outdir, pattern = "tmax", full.names = TRUE)))
prec <- as.list(intersect(list.files(path = outdir, pattern = "*.tif", full.names = TRUE), list.files(path = outdir, pattern = "prec", full.names = TRUE)))
# read the files
raster_tmin <- lapply(tmin, raster)
raster_tmax <- lapply(tmax, raster)
raster_prec <- lapply(prec, raster)

# if we are subsetting the data into our region of interest
BRA <- subset(wrld_simpl, NAME=="Brazil") # get Brazil mask
# coord_oi = c(xmin, xmax, ymin, ymax)
# e = as(extent(coord_oi), 'SpatialPolygons')
# crs(e) = "+proj=longlat +datum=WGS84 +no_defs"
df_tmin = pbapply::pblapply(raster_tmin, terra::crop, BRA) # subsects all rasters to area of interest
df_tmin = pbapply::pblapply(raster_tmax, terra::crop, BRA) # subsects all rasters to area of interest
df_tmin = pbapply::pblapply(raster_prec, terra::crop, BRA) # subsects all rasters to area of interest

df_tmin = pbapply::pblapply(raster_tmin, terra::mask, BRA) # subsects all rasters to area of interest
df_tmin = pbapply::pblapply(raster_tmax, terra::mask, BRA) # subsects all rasters to area of interest
df_tmin = pbapply::pblapply(raster_prec, terra::mask, BRA) # subsects all rasters to area of interest

df_clim = lapply(df_clim, as.data.frame, xy=T)



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
    agbd = nmv[grepl("/agbd$", nmv)]
    response_limit_flag = nmv[grepl("/response_limit_flag$", nmv)]

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
l4 <- c(paste0("./GEDI_raw/", list.files("./GEDI_raw/", pattern = '.h5')))

level4a_h5 <- hdf5r::H5File$new(l4[[1]], mode = "r")#mode="r" indica aprire in sola lettura
groups_id <- grep("BEAM\\d{4}$", gsub("/", "", hdf5r::list.groups(level4a_h5,recursive = F)), value = T)

beam = rep(i, length(level4a_i[["shot_number"]][])),
#       if (any(hdf5r::list.datasets(level4a_i) == "shot_number")) {

hdf5r::list.datasets(level4a_i) == "shot_number"

rhs <- data.table::data.table(level4a_h5[["shot_number"]][])



rasterFromXYZ_irr <- function(df, ref_raster, field){   # dataframe to be converted, reference raster
  proj <- "+proj=longlat +elips=WGS84"
  presence <- terra::vect(df[,c("lon", "lat", field)], crs = proj)
  # sum within cells, using field "value"
  convertedRaster <- terra::rasterize(presence, ref_raster, field = field, fun = 'max')
  return(convertedRaster)
}

df = GEDI
field = 'agbd'
presence <- terra::vect(df[,c("lon", "lat", field)], crs = proj)
# sum within cells, using field "value"
convertedRaster <- terra::rasterize(presence, ref_raster, field = field, fun = 'max')
#Error in p@ptr$read(x, layer, query, extent, filter, proxy, what) : 
#  Expecting a single string value: [type=character; extent=3].
GEDI_raster <- rasterFromXYZ_irr(GEDI, age_raster, "agbd")
soil_raster <- rasterFromXYZ_irr(soil, age_raster, 'type')


regrowth_list <- c()
for (i in 1:length(files_tmp)){
  regrowth_list <- c(regrowth_list, rast(files_tmp[i]))
  print(i)
}