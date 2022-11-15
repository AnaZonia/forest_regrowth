###################

library(raster) # For convenience, shapefile function and show method for Spatial objects
library(rgeos)
library(magrittr)
library(sf)
library(rgdal)
library(tmap)
library(ggplot2)
library(mapview)

setwd("/home/aavila/Documents/forest_regrowth")

panama <- sf::st_read( 
  dsn= paste0(getwd(),"/2021map/") , 
  layer="CoberturaBoscosaUsoSuelo_2021_25k"
)

# tm_shape(panama) +
#   tm_fill("Categoria",style="cat",palette="Paired") +
#     tm_layout(legend.outside = TRUE,
#             legend.width = 3,
#             legend.text.size = 25,
#             legend.outside.size = 1,
#             legend.outside.position = "bottom",
#             main.title = "Panama 2021 LULC",
#             main.title.position = "centre") +
#   tmap_options(max.categories = 33)


cafe = panama[panama$Categoria == "Café",]
pasto = panama[panama$Categoria == "Pasto",]
anual = panama[panama$Categoria == "Otro cultivo anual",]
permanente = panama[panama$Categoria == "Otro cultivo permanente",]
cana = panama[panama$Categoria == "Caña de azúcar",]
banana = panama[panama$Categoria == "Plátano/banano" ,]
citric = panama[panama$Categoria == "Cítrico" ,]
maiz = panama[panama$Categoria == "Maíz" ,]

class(cafe)

pasto_half = st_crop(pasto, c(xmin=273983.2, xmax=595076.6, ymin=796580, ymax=1066220))

mapview(pasto_half, col.regions = sf.colors(10))