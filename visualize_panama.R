###################

library(raster) # For convenience, shapefile function and show method for Spatial objects
library(rgeos)
library(magrittr)
library(sf)
library(rgdal)
library(tmap)
library(ggplot2)

setwd("/home/aavila/Documents/forest_regrowth")

panama <- sf::st_read( 
  dsn= paste0(getwd(),"/2021map/") , 
  layer="CoberturaBoscosaUsoSuelo_2021_25k"
)


jpeg(file="panama_2021_lulc.jpeg")

tm_shape(panama) +
  tm_fill("Categoria",style="cat",palette="Paired") +
  tm_layout(legend.outside = TRUE,
            legend.width = 3,
            legend.text.size = 25,
            legend.outside.size = 1,
            legend.outside.position = "bottom",
            main.title = "Panama 2021 LULC",
            main.title.position = "centre") +
  tmap_options(max.categories = 33)

dev.off()
