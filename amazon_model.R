####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Sept 2022
####################################################################

library(raster) # handling spatial data
library(geodata) # to extract worldclim with getData
library(sp) # to extract worldclim with getData
library(terra)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("grimbough/rhdf5")
library(rhdf5)
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R)
library(stringr)
library(tidyverse)
library(plyr)
library(foreach)
library(doParallel)
## Brazil shapefile mask
library(maptools)  ## For wrld_simpl
data(wrld_simpl)
SPDF <- subset(wrld_simpl, NAME=="Brazil")

setwd("/home/aavila/Documents/forest_regrowth")

#specifying coordinates of interest - where would you like to crop the dataframe to?
xmin = -48.89171
xmax = -46.41909
ymin = -3.837333
ymax = -2.41036

####################################################################
########## FUNCTIONS ##########
####################################################################

# find utm zone from longitude - allows to analyze data spanning multiple UTM zones
long2UTMzone = function(long) {
  ## Function to get the UTM zone for a given longitude
  (floor((long + 180)/6) %% 60) + 1
}

# Convert coodinates from lat/long to UTM.
# this allows merging dataframes with different coordinates (converting all to 30m resolution)
LongLatToUTM = function(x,y){
  xy = data.frame(x = x, y = y)
  xy$zone = long2UTMzone(x)

  # split the dataframe by UTM zones
  list_zones = split(xy, xy$zone)

  res_list = list()

  #must convert coodinates separately for different zones
  for (i in 1:length(list_zones)){
    z = list_zones[[i]][1,ncol(list_zones[[i]])] #obtain zone value
    coordinates(list_zones[[i]]) = c("x", "y")
    proj4string(list_zones[[i]]) = CRS("+proj=longlat +datum=WGS84") #convert to spatial object
    # EPSG code calculated for the southern hemisphere as 32700+UTM zone
    # add converted coordinates back into the list
    # obtain list of SpatialObjects, one per UTM zone, with the UTM coordinates.
    res_list = append(res_list, spTransform(list_zones[[i]], CRS(paste("+proj=utm +zone=", z, " +init=epsg:327", z, sep=''))) )
  }

  #convert SpatialObjects into data frames
  res_list = lapply(res_list, as.data.frame)

  #if our data spans more than one UTM zone, res_list will have more than one element.
  #in this case, we unite them all into a single dataframe with zone, x and y columns.
  if (length(res_list) != 1){
    for (i in 2:length(res_list)){
    res_list[[1]] = rbind(res_list[[1]], res_list[[i]])
    }
  }
  
  #convert all coordinates to the nearest multiple of 30 (nearest pixel present in Landsat resolution)
  res_list[[1]]$x = round( res_list[[1]]$x/30 ) * 30
  res_list[[1]]$y = round( res_list[[1]]$y/30 ) * 30

  result = res_list[[1]][ order(as.numeric(row.names(res_list[[1]]))), ]

  #returns dataframe with zone, x and y columns.
  return(result)
}

# unifies all dataframes in a list into a single dataframe with one column per year.
df_merge = function(df){
  for (i in 2:length(df)){
    df[[1]] = cbind(df[[1]], df[[i]][3])
  }
  return(df[[1]])
}

# creating one processed dataframe from multiple raw Mapbiomas .tif files.
# file is the path to the directory, and crop a Boolean (whether we are subsetting the data into our region of interest)
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

# reduce date from "%Y-%m-%d %H:%M:%S" format into just the year
extract_year = function(df, in_format, out_format){
  df$date = as.POSIXct(df$date, format = in_format)
  df$date = format(df$date, format=out_format) 
  return(df)
}

# find last instance of a value in a dataframe, and returns the year (column name) of that occurrence
find_last_instance = function(df, fun){
  instances = sapply(apply(df[3:(ncol(df)-3)], 1, fun),names) # returns rowname/index as well as years flagged as FIRE events
  last = lapply(lapply(instances, unlist), max) # select for only the most recent observed regrowth
  last_df = as.data.frame(unlist(last))
  last_df = as.data.frame(lapply(last_df,as.numeric))
  last_df[is.na(last_df)] = 0
  return(last_df)
}

####################################################################
########## SWITCHES ##########
####################################################################
import_mapbiomas = T
import_clim == T
import_dubayah = F
import_santoro = F
#import_potapov = F

####################################################################
########## EXTRACTING DATA ##########
####################################################################

##########  REGROWTH/DEFORESTATION, FIRE AND LAND USE ##########

# Mapbiomas data was manually downloaded.
if (import_mapbiomas == T){
  regrowth = readRDS("./scripts/regrowth.rds")
  fire = readRDS("./scripts/fire.rds")
  lulc = readRDS("./scripts/lulc.rds")
}else{
  regrowth = making_df('./mapbiomas/regrowth', crop=F)
  # subset only for the pixels that show some regrowth history
  regrowth = regrowth[rowSums(sapply(regrowth, '%in%', seq(500, 599, by=1))) > 0,]
  saveRDS(regrowth, "regrowth.rds")

  fire = making_df('./mapbiomas/fire', crop=F)
  fire = fire[rownames(regrowth),] # subset only the coordinates that have some history of regrowth.
  saveRDS(fire, "fire.rds")

  lulc = making_df('./mapbiomas/lulc', crop=T)
  lulc = lulc[rownames(regrowth),] # subset only the coordinates that have some history of regrowth.
  saveRDS(lulc, "lulc.rds")
}

########### NOTE: Had issues using making_df() for the newly imported files, so I rewrote the code:
# A few files were corrupted. Need to redownload the data.

files = list.files(path = './mapbiomas/regrowth_amazon', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 

regrowth_list = c()
for (i in 1:length(files)){
  regrowth_list[[i]] = try(raster(files[i]))
  print(i)
}

tmp_dfs <- discard(regrowth_list, inherits, 'try-error')

#selecting for southern amazon

#extent = extent(c(xmin = -65, xmax = -45, ymin = -15, ymax = 0))
#tst = crop(tmp_dfs[[1]], extent)

tmp_dfs = lapply(tmp_dfs, crop, xy=T)


tmp_dfs = lapply(tmp_dfs, as.data.frame, xy=T)



merged_df = df_merge(tmp_dfs)

colnames(merged_df) = str_sub(colnames(merged_df), start= -4)   # makes column names only "yyyy"

merged_df = merged_df[order(merged_df$x),]   #order by longitude, so the pixels are separated by UTM zones.

merged_df = cbind(merged_df, LongLatToUTM(merged_df$x, merged_df$y))   # converts lat, long coordinates to UTM





##########  TEMPERATURE AND RAINFALL ##########

if (import_clim == T){

  climate_data_import = function(df){
    colnames(df) = str_sub(colnames(df), start= -7)   # makes column names only "yyyy.mm"
    result = t(apply(df[3:ncol(df)], 1, tapply, gl(34, 12), mean))
    colnames(result) = c(1985:2018)
    return(as.data.frame(cbind(df[,1:2], result)))
  }

  df_prec = readRDS('df_prec.rds')
  df_tmin = readRDS('df_tmin.rds')
  df_tmax = readRDS('df_tmax.rds')

  prec = climate_data_import(df_prec)
  tmax = climate_data_import(df_tmax)
  tmin = climate_data_import(df_tmin)
  temp = (tmax + tmin) / 2

  climate_data_cleanup = function(df){
    df$mean = colMeans(df[3:ncol(df)])
    df = cbind(df, LongLatToUTM(df$x, df$y))
    df = df[,3:ncol(df)]
    df$xy = paste0(df$zone, df$x, df$y)
    return(df)
  }

  prec = climate_data_cleanup(prec)
  temp = climate_data_cleanup(temp)

  #resolution is about 4.5 km.
  #Need to add temperature data for ALL cells within 4.5km of center point.





}else{

  outdir <- "./worldclim_brazil" # Specify your preferred working directory
  vars = c("tmin_", "tmax_", "prec_")
  ranges = c("1980-1989", "1990-1999", "2000-2009", "2010-2018")

  url = paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", do.call(paste0, expand.grid(vars, ranges)), ".zip")
  zip <- file.path(outdir, basename(url))

  for (i in 1:length(url)){
    download.file(url[i], zip[i])
  }

  # get all the zip files
  zipF <- list.files(path = outdir, pattern = "*.zip", full.names = TRUE)
  # unzip all your files
  ldply(.data = zipF, .fun = unzip, exdir = outdir)
  # get the csv files
  #remove all datapoints before landsat (1985)
  file.remove(list.files(path = outdir, pattern = "1980|1981|1982|1983|1984", full.names = TRUE))

  clim_files <- as.list(list.files(path = outdir, pattern = "*.tif", full.names = TRUE))
  # read the csv files
  raster_clim <- lapply(clim_files, raster)
  
  # if we are subsetting the data into our region of interest
  coord_oi = c(xmin, xmax, ymin, ymax) #this specifies the coordinates of Paragominas.
  e = as(extent(coord_oi), 'SpatialPolygons')
  crs(e) = "+proj=longlat +datum=WGS84 +no_defs"
  df_clim = lapply(raster_clim, crop, e) # subsects all rasters to area of interest

  df_clim = lapply(df_clim, as.data.frame, xy=T)

  #we have 409 entries for prec, tmax and tmin.
  df_prec = df_merge(df_clim[1:408])
  df_tmax = df_merge(df_clim[410:817])
  df_tmin = df_merge(df_clim[819:c(length(df_clim))])

  saveRDS(df_prec, "df_prec.rds")
  saveRDS(df_tmax, "df_tmax.rds")
  saveRDS(df_tmin, "df_tmin.rds")
}

##########  SOIL TYPE ##########



##########  BIOMASS ##########

# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
# 1km resolution, 2019-2021
if (import_dubayah == T){
  biomass = readRDS("./biomass_dubayah.rds")
}else{
  
  lat_lon_oi = c(ymax, ymin, xmin, xmax)  #lat and lon of interest

  GEDI_download = l4_download(
    lat_lon_oi,
    outdir = "./GEDI",
    from = "2020-01-01",
    to = "2020-07-31",
    just_path = F
  )

  filepaths = paste0('./GEDI/', list.files("./GEDI"))
  GEDI_list = lapply(filepaths, l4_get, just_colnames = F)
  GEDI_list = lapply(GEDI_list, subset, select=c("date", "lat_lowestmode", "lon_lowestmode", "agbd_se", "agbd"))
  GEDI_list = lapply(GEDI_list, extract_year, in_format = "%Y-%m-%d %H:%M:%S", out_format = "%Y")

  select_range = function(df){
    df = df[ymin < lat_lowestmode & lat_lowestmode < ymax & xmin < lon_lowestmode & lon_lowestmode < xmax,]
    return(df)
  }

  GEDI_list = lapply(GEDI_list, select_range)

  for (i in 2:length(GEDI_list)){
  GEDI_list[[1]] = rbind(GEDI_list[[1]], GEDI_list[[i]])
  }

  biomass = GEDI_list[[1]]

  biomass = cbind(biomass, LongLatToUTM(biomass$lon_lowestmode, biomass$lat_lowestmode))

  saveRDS(biomass, "biomass_dubayah.rds")
}

# Santoro et al 2018 data -> GlobBiomass ESA (Mg/ha)
# 100m resolution, 2010
if (import_santoro == T){
  biomass = readRDS("biomass_santoro.rds")
}else{
  biomass1 = raster("N00W060_agb.tif")
  biomass2 = raster("N00W100_agb.tif")
  biomass3 = raster("N40W060_agb.tif")

  biomass = merge(biomass1, biomass2, biomass3)

  ## crop and mask
  r2 <- crop(biomass, extent(SPDF))
  r3 <- mask(r2, SPDF)

  biomass = as.data.frame(biomass, xy = TRUE)

  #biomass = biomass[ymin < biomass$y & biomass$y < ymax & xmin < biomass$x & biomass$x < xmax,]

  biomass = cbind(biomass, LongLatToUTM(biomass$x, biomass$y))

  colnames(biomass) = c('lon', 'lat', 'agbd', 'zone', 'x', 'y')
  saveRDS(biomass, "biomass_santoro_Brazil.rds")
}

# Potapov et al 2020 -> GLAD Forest Canopy Height (m)
# 30m resolution, 2019
if (import_potapov == T){
  biomass = readRDS("biomass_potapov.rds")
}else{
  biomass = raster("Forest_height_2019_SAM.tif")

  ## crop and mask
  r2 <- crop(biomass, extent(SPDF))

  #writeRaster(r3, "Forest_height_2019_Brazil.tif")

cores <- 50
cl <- makeCluster(cores) #output should make it spit errors
registerDoParallel(cl)

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
  writeRaster(mask(split_list[i], SPDF))
}

####################################################################
########## ASSEMBLING THE DATASET ##########
####################################################################

########## REGROWTH ##########
# VARIABLES OBTAINED
# forest_age (time of last observed regrowth)

# INDEX #
# 100 anthropic
# 200 mature
# 300 secondary
# 400 deforestation
# 500 regrowth
# 600 secondary suppression

if (import_santoro == T){
  regrowth = cbind(regrowth[,1:25], regrowth[,(ncol(regrowth)-2):ncol(regrowth)])
  fire = cbind(fire[,1:28], fire[,(ncol(fire)-2):ncol(fire)])
  lulc = cbind(lulc[,1:28], lulc[,(ncol(lulc)-2):ncol(lulc)])
}

# select only years that have regrowth that hasn't been suppressed.

regrowth_last_instance = find_last_instance(regrowth, function(x) which(x == 503))
colnames(regrowth_last_instance) = "last_regrowth"

suppression_last_instance = find_last_instance(regrowth, function(x) which(600 <= x & x < 700))
colnames(suppression_last_instance) = "last_suppression"


regrowth_unsuppressed = cbind(regrowth[(ncol(regrowth)-2):ncol(regrowth)],regrowth_last_instance)
regrowth_unsuppressed = subset(regrowth_unsuppressed, regrowth_unsuppressed$last_regrowth-suppression_last_instance > 0)
regrowth_unsuppressed = cbind(regrowth_unsuppressed[,(ncol(regrowth_unsuppressed)-2):ncol(regrowth_unsuppressed)], 'forest_age' = max(regrowth_unsuppressed$last_regrowth)-regrowth_unsuppressed$last_regrowth)

#############################


# this removes any instance of -600 coming after -500
regrowth2 = regrowth[rownames(regrowth) %in% rownames(regrowth_unsuppressed), ] 
regrowth2[regrowth2 < 400 & regrowth2 >= 300] = 0 #secondary
regrowth2[100 <= regrowth2 & regrowth2 < 200] = 1 #anthropic
regrowth2[regrowth2 == 515] = 1 #anthropic - fixing typos
regrowth2[regrowth2 == 215] = 1 #anthropic - fixing typos

regrowth2[regrowth2 > 700 & regrowth2 < 800] <- NA # remove values that account for misclassification
regrowth2[regrowth2 > 400 & regrowth2 < 500] <- NA # remove values that account for urban areas/misclassification
regrowth2 = regrowth2[complete.cases(regrowth2),]

# now, we remove pixels that show unflagged moments of regrowth or repression.
tmp = cbind(NA, regrowth2[,3:(ncol(regrowth)-3)])
tmp2 = cbind(regrowth2[,3:(ncol(regrowth)-3)],NA)
tmp3 = tmp-tmp2

# -1 is evidence of unflagged repression
# 1 is unflagged regrowth
tmp3[tmp3 == 1] <- NA
tmp3[tmp3 == -1] <- NA
tmp3 = tmp3[,2:(ncol(tmp)-1)]
tmp3 = tmp3[complete.cases(tmp3),]

#selects for cleaned rows
regrowth = regrowth[rownames(regrowth) %in% rownames(tmp3), ] 

regrowth_last_instance = find_last_instance(regrowth, function(x) which(x == 503))
colnames(regrowth_last_instance) = "last_regrowth"

regrowth_cleaned = cbind(regrowth[(ncol(regrowth)-2):ncol(regrowth)],regrowth_last_instance)
regrowth_cleaned = cbind(regrowth_cleaned[,1:3], 'forest_age' = max(regrowth_cleaned$last_regrowth)-regrowth_cleaned$last_regrowth)

regrowth_cleaned$xy = paste0(regrowth_cleaned$zone, regrowth_cleaned$x, regrowth_cleaned$y)
biomass$xy = paste0(biomass$zone, biomass$x, biomass$y)

agb_forest_age = cbind(regrowth_cleaned, agbd = biomass[match(regrowth_cleaned$xy,biomass$xy),c("agbd")])
agb_forest_age = agb_forest_age[complete.cases(agb_forest_age[, ncol(agb_forest_age)]), ]

plot(agb_forest_age$forest_age, agb_forest_age$agbd)

# how to add distance to nearest mature forest?

sds = aggregate(agbd ~ forest_age, agb_forest_age, sd)
means = aggregate(agbd ~ forest_age, agb_forest_age, mean)
sum_stats = cbind(means, sds[,2])
colnames(sum_stats) = c('age', 'mean', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = mean)) + 
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd)) +
  geom_point()




########## LAND USE ##########
# VARIABLES OBTAINED
# number of years under each land use tnrype
# time since last observation of each land use type

# INDEX ## 3 = forest
# 15 = pasture
# 39 = soy
# 46 = coffee
# 20 = sugar cane
# 41 = other annual crop
# 48 = other perennial crop

# total years under each land use type
lulc$pasture = rowSums(lulc == 15)
lulc$soy = rowSums(lulc == 39)
lulc$coffee = rowSums(lulc == 46)
lulc$sugar = rowSums(lulc == 20)
lulc$other_perennial = rowSums(lulc == 48)
lulc$other_annual = rowSums(lulc == 41)


# time since last observation of each land use type
ts_pasture = find_last_instance(lulc, function(x) which(x == 15))
lulc$ts_pasture = max(ts_pasture)-ts_pasture

ts_soy = find_last_instance(lulc, function(x) which(x == 39))
lulc$ts_soy = max(ts_soy)-ts_soy

#ts_coffee = find_last_instance(lulc, function(x) which(x == 46))
#ts_sugar = find_last_instance(lulc, function(x) which(x == 20))

ts_other_perennial = rowSums(lulc == 48)
lulc$ts_other_perennial = max(ts_other_perennial)-ts_other_perennial

ts_other_annual = rowSums(lulc == 41)
lulc$ts_other_annual = max(ts_other_annual)-ts_other_annual

lulc$xy = paste0(lulc$zone, lulc$x, lulc$y)

########## FIRE ##########

#count number of total fires
fire$num_fires = rowSums(fire[3:(ncol(fire)-3)])
fire_last_instance = find_last_instance(fire, function(x) which(x == 1))
colnames(fire_last_instance) = "last_burn"
fire = cbind(fire,fire_last_instance)
fire$last_burn = max(fire$last_burn) - fire$last_burn
fire$last_burn[fire$last_burn == max(fire$last_burn)] = NA
fire$xy = paste0(fire$zone, fire$x, fire$y)

central_df = cbind(agb_forest_age, last_burn = fire[match(agb_forest_age$xy,fire$xy),c("last_burn")])
central_df = cbind(central_df, num_fires = fire[match(central_df$xy,fire$xy),c("num_fires")])
central_df = cbind(central_df, pasture = lulc[match(central_df$xy,lulc$xy),c("pasture")])
central_df = cbind(central_df, soy = lulc[match(central_df$xy,lulc$xy),c("soy")])
central_df = cbind(central_df, other_perennial = lulc[match(central_df$xy,lulc$xy),c("other_perennial")])
central_df = cbind(central_df, other_annual = lulc[match(central_df$xy,lulc$xy),c("other_annual")])
central_df = cbind(central_df, tavg = temp[match(central_df$xy,temp$xy),c("mean")])
central_df = cbind(central_df, prec = prec[match(central_df$xy,prec$xy),c("mean")])
central_df$last_burn[is.na(central_df$last_burn)] = 1

central_df = lapply(central_df, as.numeric)
################# passing into the model ##################

# pars[2]*tavg+pars[3]*prec+
# central_df$tavg, central_df$prec, 
#  0.001, 0.0001, 




G = function(pars, other_perennial, other_annual, pasture, soy, forest_age, num_fires, last_burn) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  k = pars[2]*other_perennial+pars[3]*pasture+pars[4]*soy+pars[5]*other_annual + pars[6]*other_annual + pars[7]*num_fires^(pars[8]*last_burn)
  # Prediction of the model
  Gmax * (1 - exp(-k*(forest_age)))
}

####################
#finding the right parameters
par0 = c(50, 0.0001,  0.0001,  0.001, 0.001,  0.001, 0.1, 0.1)

Gpred = G(par0, central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

# dat$tavg, dat$prec, 

NLL = function(pars, dat) {
  # Values prediced by the model
  if(pars[9] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars, dat$other_perennial, dat$other_annual, dat$pasture, dat$soy, dat$forest_age, dat$num_fires, dat$last_burn)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = dat$agbd, mean = Gpred, sd = pars[9], log = TRUE))
  #print(pars)
  return(fun)
}

par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.1, 0.1)
NLL(par0, dat=central_df)

o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "CG")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG")
for (i in meth0){
  o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

#  central_df$tavg, central_df$prec, 

pred = G(o$par[1:9], central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

plot(central_df$agbd, pred, abline(0,1))



########################################################################################


biomass[biomass$agbd == 49 & biomass$x == 213210 & biomass$y == 9663510,]

regrowth[rownames(regrowth) == 11768979, ] 


######################### WHAT'S WRONG
# 34182406   23 233190 9622260          1 232331909622260 2205.453
# 15229309   23 306570 9684000         30 233065709684000 0.0192384105
regrowth[rownames(regrowth) == 40507664, ] 
# -48.85896 -3.60705
# 449.0891 @ 3 yrs old
# -3.607193 -48.85909

tst = subset(amazon, agbd > 350)
# 18412814   23 300000 9673650          2  233e+059673650  435.5806 435.5805969
tst = subset(tst, forest_age == 1)

# 30554059   22 773460 9634020         22 227734609634020    0

biomass |>
  subset(abs(N00W060_agb - 67) < 1e-3)


biomass |>
  subset(abs(Lon - -100.7) < 1e-5 & abs(Lat - 59.6) < 1e-5)

Q <- quantile(biomass$agbd, probs=c(.1, .9), na.rm = FALSE)
iqr <- IQR(biomass$agbd)
up <-  Q[2]+1.5*iqr # Upper Range  
eliminated<- subset(biomass, biomass$agbd < (Q[2]+1.5*iqr))

