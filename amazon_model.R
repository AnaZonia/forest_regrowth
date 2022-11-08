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

  files = list.files(path = './mapbiomas/regrowth', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 

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


####################################################################
########## SWITCHES ##########
####################################################################
import_mapbiomas = T
import_GEDI = T

####################################################################
########## EXTRACTING DATA ##########
####################################################################

# REGROWTH/DEFORESTATION, FIRE AND LAND USE.
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

# biomass data
if (import_GEDI == T){
  biomass = readRDS("./biomass.rds")
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

  saveRDS(biomass, "biomass.rds")
}



# TEMPERATURE AND RAINFALL

# Downloading and unzipping

#if (import_clim == T){
#}else{
  outdir <- "./worldclim_brazil" # Specify your preferred working directory
  # vars = c("tmin_", "tmax_", "prec_")
  # ranges = c("1980-1989", "1990-1999", "2000-2009", "2010-2018")

  # url = paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", do.call(paste0, expand.grid(vars, ranges)), ".zip")
  # zip <- file.path(outdir, basename(url))

  # for (i in 1:length(url)){
  #   download.file(url[i], zip[i])
  # }

  # get all the zip files
  #zipF <- list.files(path = outdir, pattern = "*.zip", full.names = TRUE)
  # unzip all your files
  #ldply(.data = zipF, .fun = unzip, exdir = outdir)
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
  df_prec = df_clim[1:409]
  df_tmax = df_clim[410:818]
  df_tmin = df_clim[819:length(df_clim)]

  prec = df_merge(df_prec)
  tmax = df_merge(df_tmax)
  tmin = df_merge(df_tmin)

  colnames(prec) = str_sub(colnames(prec), start= -7)   # makes column names only "yyyy"
  colnames(tmax) = str_sub(colnames(tmax), start= -7)   # makes column names only "yyyy"
  colnames(tmin) = str_sub(colnames(tmin), start= -7)   # makes column names only "yyyy"

  # reduce date from "%Y-%m-%d %H:%M:%S" format into just the year
  extract_year = function(df, in_format, out_format){
  colnames(prec[3:ncol(prec)]) = as.POSIXct(df$date, format = in_format)
  df$date = format(df$date, format=out_format) 
  return(df)



#}




#tst = as.POSIXct(colnames(, format = "%Y.%m")

#as.Date(paste(prec[3:ncol(prec)]),1,sep="."),"%Y.%m.%d")

#class(colnames(prec[3:ncol(prec)]))



tavg = readRDS("./mapbiomas/tavg")
prec = readRDS("./mapbiomas/prec")
tavg$mean = rowMeans(tavg[,c(3:14)])
prec$mean = rowMeans(prec[,c(3:14)])
prec = cbind(prec, LongLatToUTM(prec$x, prec$y))
tavg = cbind(tavg, LongLatToUTM(tavg$x, tavg$y))
prec = prec[,-c(1:2)]
tavg = tavg[,-c(1:2)]

# SOIL TYPE


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

regrowth2 = regrowth[rownames(regrowth) %in% rownames(amazon), ] 
regrowth2[regrowth2 < 400 & regrowth2 >= 300] = 1 #secondary
regrowth2[100 < regrowth2 & regrowth2 < 200] = 0 #anthropic


tmp = cbind(regrowth2[,3:(ncol(regrowth)-3)],NA)
tmp2 = cbind(NA, regrowth2[,3:(ncol(regrowth)-3)])
res = apply(tmp-tmp2, 1, function(x) {any(x == -1)})
res[is.na(res)] = FALSE

tst = regrowth2[res,]


# within these, find 2019-year of last pixel flagged as "regrowth" -> this will give forest age.
regrowth_instances = sapply(apply(tst[3:(ncol(tst)-4)], 1, function(x) which(600 > x & x >= 500)),names) # returns rowname/index as well as years flagged as regrowth events
suppression_instances = sapply(apply(regrowth[3:(ncol(regrowth)-4)], 1, function(x) which(600 <= x & x < 700)),names) # returns rowname/index as well as years flagged as regrowth events

last_regrowth = lapply(lapply(regrowth_instances, unlist), max) # select for only the most recent observed regrowth
last_regrowth_df = as.data.frame(unlist(last_regrowth))
last_regrowth_df = as.data.frame(lapply(last_regrowth_df,as.numeric))
last_regrowth_df = cbind(regrowth[(ncol(tst)-2):ncol(tst)],last_regrowth_df)

suppression_instances = lapply(lapply(suppression_instances, unlist), max) # select for only the most recent observed regrowth
suppression_instances_df = as.data.frame(unlist(suppression_instances))
suppression_instances_df[is.na(suppression_instances_df)] = 0
suppression_instances_df = as.data.frame(lapply(suppression_instances_df,as.numeric))

amazon = subset(last_regrowth_df, last_regrowth_df[,4]-suppression_instances_df > 0)

amazon = cbind(amazon[,1:3], 'forest_age' = 2019-amazon[,4]) # since AGB data is from 2019

#############################




########## LAND USE ##########
# VARIABLES OBTAINED
# number of years under each land use type
# time since last observation of each land use type

# INDEX ## 3 = forest
# 15 = pasture
# 39 = soy
# 46 = coffee
# 20 = sugar cane
# 41 = other annual crop
# 48 = other perennial crop

amazon$pasture = rowSums(lulc == 15)
amazon$soy = rowSums(lulc == 39)
amazon$coffee = rowSums(lulc == 46)
amazon$sugar = rowSums(lulc == 20)
amazon$other_perennial = rowSums(lulc == 48)
amazon$other_annual = rowSums(lulc == 41)

########## FIRE ##########

#count number of total fires
fire$num_fires = rowSums(fire[3:(ncol(fire)-3)])


########## ENVIRONMENTAL VARIABLES ##########

amazon$tavg = mean(tavg$mean) #******* all the same since we're looking at only a city for now, to simplify
amazon$prec = mean(prec$mean) #******* find temp and prec data for Brazil with time resolution

amazon = amazon[,3:ncol(amazon)]
amazon = amazon[,-c(1:2)]
prec = prec[,-c(ncol(prec))]

amazon$xy = paste0(amazon$zone, amazon$x, amazon$y)
biomass$xy = paste0(biomass$zone, biomass$x, biomass$y)
regrowth$xy = paste0(regrowth$zone, regrowth$x, regrowth$y)



amazon = cbind(amazon, agbd = biomass[match(amazon$xy,biomass$xy),c("agbd")])
agb_amazon = amazon[complete.cases(amazon[, ncol(amazon)]), ]

plot(agb_amazon$forest_age, agb_amazon$agbd)

# repression without flagging is still an issue. fix it with jacqueline's code




################# passing into the model ##################


G = function(pars, tavg, rain, other_perennial, other_annual, pasture,soy, forest_age) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  k = pars[2]*tavg+pars[3]*rain+pars[4]*other_perennial+pars[5]*pasture+pars[6]*soy+pars[7]*other_annual
  # Prediction of the model
  Gmax * (1 - exp(-k*(forest_age)))
}

####################
#finding the right parameters
par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001, 0.0001, 0.1, 0.1)
val = G(par0, mean(complete$tavg), mean(complete$prec), mean(complete$other_perennial), mean(complete$other_annual), mean(complete$pasture), mean(complete$soy), mean(complete$forest_age))

NLL = function(pars, dat) {
  # Values prediced by the model
  if(pars[9] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars, dat$tavg, dat$prec, dat$other_perennial, dat$other_annual, dat$pasture, dat$soy, dat$forest_age)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = dat$agbd, mean = Gpred, sd = pars[9], log = TRUE))
  #print(pars)
  return(fun)
}

par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.1, 0.1)
NLL(par0, dat=complete)

o = optim(par = par0, fn = NLL, dat = complete, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "SANN")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG", "SANN")
for (i in meth0){
  o = optim(par = par0, fn = NLL, dat = complete, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

pred = G(o$par[1:9], complete$tavg, complete$prec, complete$other_perennial, complete$other_annual, complete$pasture, complete$soy, complete$forest_age)

plot(complete$agbd, pred, abline(0,1))
