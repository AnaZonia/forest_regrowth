library(zoo)
library(raster)
library(rgdal)

setwd("/home/aavila/Documents/land_classification")

# select the raw landsat files for a single landsat footprint (over the needed time frame)
panama <- list.files(path = './landsat/bayano', pattern='\\.TIF$', full.names=TRUE)


# make a list of the raster files for a single footprint
panama_rasters = list()
for (i in 1:length(panama)){
  panama_rasters = append(panama_rasters,raster(panama[i]))
}

e = raster::extent(c(max(x_min), min(x_max), max(y_min), min(y_max)))

#range 0-255
tst = as.data.frame(panama_rasters[[1]], xy = T)
tst = subset(tst, tst$x == 734550 & tst$y == 1062270)
tst

# calculate NDMI values from original raster values
panama_NDMI = list()
for (i in seq(1, length(panama_rasters), 2)){
  print(i)
  panama_NDMI = append(panama_NDMI, (panama_rasters[[i]]-panama_rasters[[i+1]])/(panama_rasters[[i]]+panama_rasters[[i+1]]))
}


# because of cloud cover, the extent is not exactly the same for all images throughout time.
# to run through the algorithm, we need the exact same extent

x_min = c()
x_max = c()
y_min = c()
y_max = c()


for (i in 1:length(panama_NDMI)){
  x_min = append(x_min, xmin(panama_NDMI[[i]]))
  x_max = append(x_max, xmax(panama_NDMI[[i]]))
  y_min = append(y_min, ymin(panama_NDMI[[i]]))
  y_max = append(y_max, ymax(panama_NDMI[[i]]))
}

e = raster::extent(c(max(x_min), min(x_max), max(y_min), min(y_max)))

for (i in 1:length(panama_NDMI)){
  panama_NDMI[[i]] <- raster::crop(panama_NDMI[[i]], e)
}

NDMI_df = as.data.frame(panama_NDMI[[1]], xy=TRUE)

#make NDMI change df for the plot
for (i in 2:length(panama_NDMI)){
    tst = as.data.frame(panama_NDMI[[i]], xy=TRUE)
    NDMI_df = cbind(NDMI_df, tst[,3])
}

NDMI_df_complete = NDMI_df[complete.cases(NDMI_df),]

threshold = -0.2


df = NDMI_df_complete[rollapply(NDMI_df_complete[,3:ncol(NDMI_df_complete)] > threshold, width = 3, all),]

library(doParallel)

cl <-makeCluster(8, type ="PSOCK")
res.mat <- parLapply(x = 1:10,
                     fun = function(x) my.fun(x),
                     cl = cl
                     )





#######################

setwd("/home/aavila/Documents/GFW")

height = raster('Forest_height_2019_SAM.tif')







