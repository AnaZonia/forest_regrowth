####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Ana Avila - Jan 2023
# Functions
####################################################################

# Finds utm zone from longitude - allows to analyze data spanning multiple UTM zones
# This function is necessary for LongLatToUTM (below)
# It intakes:
# Longitude (numeric)
# It outputs:
# UTM zone (numeric)
long2UTMzone <- function(long) {
  ## Function to get the UTM zone for a given longitude
  (floor((long + 180)/6) %% 60) + 1
}

# Converts coodinates from lat/long to UTM.
# this allows merging dataframes with different coordinates.
# It converts all points to to 30m resolution by finding the nearest coordinate multiple of 30.
# It intakes:
# x and y = longitude and latitude, respectively (numeric or vector)
# It outputs:
# Dataframe with zone, x and y columns in UTM format.
LongLatToUTM <- function(x,y){
  xy <- data.frame(x = x, y = y)
  xy$zone <- long2UTMzone(x)

  # split the dataframe by UTM zones
  list_zones <- split(xy, xy$zone)

  res_list <- list()

  #must convert coodinates separately for different zones
  for (i in 1:length(list_zones)){
    z <- list_zones[[i]][1,ncol(list_zones[[i]])] #obtain zone value
    coordinates(list_zones[[i]]) <- c("x", "y")
    proj4string(list_zones[[i]]) <- CRS("+proj=longlat +datum=WGS84") #convert to spatial object
    # EPSG code calculated for the southern hemisphere as 32700+UTM zone
    # add converted coordinates back into the list
    # obtain list of SpatialObjects, one per UTM zone, with the UTM coordinates.
    res_list <- append(res_list, spTransform(list_zones[[i]], CRS(paste("+proj=utm +zone=", z, " +init=epsg:327", z, sep=''))) )
  }

  #convert SpatialObjects into data frames
  res_list <- lapply(res_list, as.data.frame)

  #if our data spans more than one UTM zone, res_list will have more than one element.
  #in this case, we unite them all into a single dataframe with zone, x and y columns.
  if (length(res_list) != 1){
    for (i in 2:length(res_list)){
    res_list[[1]] <- rbind(res_list[[1]], res_list[[i]])
    }
  }
  
  #convert all coordinates to the nearest multiple of 30 (nearest pixel present in Landsat resolution)
  res_list[[1]]$x <- round( res_list[[1]]$x/30 ) * 30
  res_list[[1]]$y <- round( res_list[[1]]$y/30 ) * 30

  result <- res_list[[1]][ order(as.numeric(row.names(res_list[[1]]))), ]

  #returns dataframe with zone, x and y columns.
  return(result)
}

# Unifies all dataframes in a list into a single dataframe with one column per year.
# Useful when you have a list with one dataframe per moment in time.
# It intakes:
# List of dataframes
# It outputs:
# One combined dataframe containing the THIRD column of every dataframe in the list
df_merge <- function(df){
  for (i in 2:length(df)){
    df[[1]] <- cbind(df[[1]], df[[i]][3])
  }
  return(df[[1]])
}

# reduce date from "%Y-%m-%d %H:%M:%S" format into just the year
extract_year <- function(df, in_format, out_format){
  df$date <- as.POSIXct(df$date, format = in_format)
  df$date <- format(df$date, format = out_format) 
  return(df)
}

# makes dataframe from large raster files, substituting as.data.frame()
# assumes cells in the raster start from 1 for correct assignment of coordinates.
df_from_raster <- function(raster){
  bm_test <- getValues(raster)
  bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
  bm_test <- na.omit(bm_test)
  bm_test[,c("x","y")] <- xyFromCell(raster, bm_test$cell)
  return(bm_test)
}

# find last instance of a value in a dataframe, and returns the year (column name) of that occurrence
find_last_instance <- function(df, fun){
  instances <- sapply(apply(df[3:(ncol(df)-3)], 1, fun),names) # returns rowname/index as well as years flagged as FIRE events
  last <- lapply(lapply(instances, unlist), max) # select for only the most recent observed regrowth
  last_df <- as.data.frame(unlist(last))
  last_df <- as.data.frame(lapply(last_df,as.numeric))
  last_df[is.na(last_df)] <- 0
  return(last_df)
}
