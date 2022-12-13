####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Dec 2022
####################################################################


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("grimbough/rhdf5")
library(rhdf5) # for handling raw GEDI data
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R) # for extracting raw GEDI data

setwd("/home/aavila/Documents/forest_regrowth")


# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
ul_lat = -3.68658
lr_lat = -6.98754
ul_lon = -60.94691
lr_lon = -53.48623

outdir = tempdir()
  GEDI_download = l4_download(
  ul_lat,
  lr_lat,
  ul_lon,
  lr_lon,
  outdir = outdir,
    from = "2020-01-01",
    to = "2020-07-31",
    just_path = T
  )


  filepaths = c('./GEDI_mature/GEDI04_A_2020193104720_O08946_04_T01367_02_002_02_V002.h5', './GEDI_mature/GEDI04_A_2020193231025_O08954_01_T00916_02_002_02_V002.h5')

dataname <- l4_getmulti(filepaths[1],just_colnames = T)


  GEDI_list_mature = lapply(filepaths, l4_get, just_colnames = T)


  GEDI_list = lapply(GEDI_list, subset, select=c("date", "lat_lowestmode", "lon_lowestmode", "agbd_se", "agbd"))
  GEDI_list = lapply(GEDI_list, extract_year, in_format = "%Y-%m-%d %H:%M:%S", out_format = "%Y")


hist(GEDI_list[[5]]$agbd, breaks=2000, xlim = c(0, 80)) +
 theme(text = element_text(size = 20))  


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

  saveRDS(biomass, "biomass_0000000000.0000095232.rds")


  

# Santoro et al 2018 data -> GlobBiomass ESA (Mg/ha)
# 100m resolution, 2010
if (import_santoro == T){
  biomass = readRDS("biomass_santoro_Brazil.rds")
}else{
  biomass1 = raster("./santoro/N00W060_agb.tif")
  biomass2 = raster("./santoro/N00W100_agb.tif")
  biomass3 = raster("./santoro/N40W060_agb.tif")

  biomass = merge(biomass1, biomass2, biomass3)

  ## crop and mask
  r2 <- crop(biomass, extent(BRA))
  r3 <- mask(r2, BRA) #save this somewhere
  e <- extent(-48.32644, -43.99998, -3.2823093, -0.5377764)
  r4 <- crop(biomass1, e)


  bm_test <- getValues(r4)
  bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
  bm_test <- na.omit(bm_test)
  bm_test[,c("x","y")] <- xyFromCell(r4, bm_test$cell)

  #biomass = biomass[ymin < biomass$y & biomass$y < ymax & xmin < biomass$x & biomass$x < xmax,]

  biomass = cbind(biomass, LongLatToUTM(biomass$x, biomass$y))

  biomass <- as.data.frame(r4, xy = TRUE)
  colnames(biomass) = c('lon', 'lat', 'agbd', 'zone', 'x', 'y')
  saveRDS(biomass, "biomass_test.rds")
}


# Potapov et al 2020 -> GLAD Forest Canopy Height (m)
# 30m resolution, 2019
if (import_potapov == T){
  biomass = readRDS("Forest_height_2019_Brazil.rds")
}else{
  biomass = raster("Forest_height_2019_SAM.tif")

  ## crop and mask
  r2 <- crop(biomass, extent(BRA))
  r3 <- mask(r2, BRA)

  #writeRaster(r3, "Forest_height_2019_Brazil.tif")

}


agb_forest_age = readRDS('agb_forest_age.rds')

plot(agb_forest_age$forest_age, agb_forest_age$agbd)

agb_forest_age_2 = subset(agb_forest_age, agbd > 100)

plot(agb_forest_age_2$forest_age, agb_forest_age_2$agbd)



agb_forest_age = readRDS('agb_forest_age.rds')

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

############################################ POTAPOV #############################################


e <- extent(xmin, xmax, ymin, ymax)

biomass = raster("Forest_height_2019_Brazil.tif")
biomass_cropped <- crop(biomass,e)
#biomass_df <- as.data.frame(biomass_cropped, xy=T, na.rm = TRUE)
bm_test <- values(biomass_cropped)

bm_tst_complete <- na.omit(bm_test)


bm_test <- getValues(biomass_cropped)
bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
bm_test <- na.omit(bm_test)
bm_test[,c("x","y")] <- xyFromCell(biomass_cropped, bm_test$cell)

biomass = cbind(biomass_with_data, LongLatToUTM(biomass_with_data$x, biomass_with_data$y))

saveRDS(biomass, "biomass_potapov_tst.rds")
