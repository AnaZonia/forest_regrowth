library(GEDI4R)
#setwd("/home/aavila/forest_regrowth/GEDI_italy_tst")
outdir = tempdir()
#using Italy as bounding box for search GEDI data
# ita <- sf::st_as_sf(raster::getData('GADM', country = 'ITA', level = 1))
# #extract extent
# e <- raster::extent(ita)
# ul_lat <- e@ymax
# lr_lat <- e@ymin
# ul_lon <- e@xmin
# lr_lon <- e@xmax

# file_download <- l4_download(
#   ul_lat,
#   lr_lat,
#   ul_lon,
#   lr_lon,
#   ncore = parallel::detectCores()-1,
#   outdir = outdir,
#   from = "2020-01-01",
#   to = "2020-01-31",
#   just_path = F,
#   subset=1:4
# )

l4_zip <- system.file("extdata",
                      c("GEDI04_A_2020036151358_O06515_02_T00198_02_002_01_V002.zip",
                        "GEDI04_A_2021150031254_O13948_03_T06447_02_002_01_V002.zip"
                      ),
                      package="GEDI4R")

#Unzipping GEDI level4A data
l4 <- lapply(l4_zip,unzip,exdir = outdir)
list.files(outdir)

#list all dataset in h5 file
dataname <- l4_getmulti(l4[[1]],just_colnames = T)
head(dataname,10)
#read all footprint and merge them together.
gediL4_path <- l4_getmulti(l4,merge=T)

hist(gediL4_path$agbd,breaks <- 5000)


#select other columns to add to the default output.
#if columns are already present in the default output they will be dropped
col <-
  c("land_cover_data/leaf_off_flag",
    "agbd_pi_lower",
    "agbd_pi_upper",
    "agbd"#this will be dropped as it is already included by default in the output.
    )
#get level 4 data with the user defined column binded and with the source path of each observation
#with source=T a column with the source path for each observation will be added
gediL4 <- l4_getmulti(l4,add_col = col,source=T)

knitr::kable(head(gediL4[,c("date","tree_cover","agbd","agbd_se")]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`


outdir = tempdir()
l4_zip <- system.file("extdata",
                     "GEDI04_A_2020036151358_O06515_02_T00198_02_002_01_V002.zip",
                     package="GEDI4R")
l4 <- unzip(l4_zip,exdir = outdir)
#get GEDI data
l4_data <- l4_getmulti(l4)
#clip using vector of coordinates
b_box <- c(-50,35,52,37)
clipped <- l4_clip(l4_data,clip=b_box)
#using Shapefile to clip
bound <- system.file("extdata","bound4326.shp",package="GEDI4R")
#with  extension
clipped <- l4_clip(l4_data,clip=bound,usegeometry = F)
#with  polygon boundaries
clipped2 <- l4_clip(l4_data,clip=bound,usegeometry = T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

gediL4 <- l4_getmulti(l4)
#footprints locations and AGBD distribution against elevation
jpeg('italy.png')
l4_plotprofile(gediL4,beam_id="all")
dev.off()
