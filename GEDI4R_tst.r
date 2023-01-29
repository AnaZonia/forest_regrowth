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
#coords <- c(-0.5377764, -3.2823093, -48.32644, -43.99998) #coordinates for a test region in the northeastern Amazon
coords <- c(-1, -8.3, -56, -49) #coordinates for a test region in the northeastern Amazon
#coords <- c(max(age$lat), min(age$lat), min(age$lon), max(age$lon)) #coordinates for a test region in the northeastern Amazon


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

saveRDS(clipped, '0000000000-0000095232_GEDI.rds')
#l4_convert(clipped,4326,filename=paste0(outdir,"amazon_gedi_l4.shp"), append=FALSE)



# in the mature patch:
# 15% are smaller than 5
# 17% are smaller than 10
# 30% are smaller than 50
# 3.8% are larger than 400
# 0.7% are larger than 600

head(clipped)

tst <- subset(clipped, agbd > 600)
nrow(tst)/nrow(clipped)
hist(tst$agbd, breaks = 500)

# in the younger patch:
# 42% are smaller than 5
# 47% are smaller than 10
# 66.5% are smaller than 50
# 0.09% are larger than 400
# 0.03% are larger than 600
