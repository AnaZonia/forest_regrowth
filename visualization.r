library(ggplot2)
setwd("/home/aavila/forest_regrowth/dataframes")

#agb_forest_age <- readRDS("agb_forest_age.rds")
agb_forest_age <- readRDS('./test_files/agb_forest_age_santoro_paragominas.rds')


GEDI <- cbind(GEDI, LongLatToUTM(GEDI$lon, GEDI$lat))
GEDI$xy <- paste0(GEDI$zone, GEDI$x, GEDI$y)

agb_forest_age <- cbind(age, agbd = GEDI[match(age$xy,GEDI$xy),c("agbd")])
agb_forest_age <- agb_forest_age[complete.cases(agb_forest_age[, ncol(agb_forest_age)]), ]

# max(agb_forest_age$agbd) ----> [1] 3170.47
# there are some absurd outliers in the data.
agb_forest_age2 <- subset(agb_forest_age, agbd < 500)

plot(agb_forest_age2$forest_age, agb_forest_age2$agbd)

# saveRDS(agb_forest_age, 'agb_forest_age.rds')
# agb_forest_age = readRDS('agb_forest_age.rds')

ggplot(data = agb_forest_age, aes(x = forest_age, y = agbd)) + 
  geom_point() +
  stat_summary(geom = "point", fun = "mean", col = "black",
    size = 6, shape = 23,fill = "red")

# as we can see, the data is still really oddly distributed.

#GEDI_mid_amazon <- readRDS('GEDI_midamazon_dfunified.rds')
# > range(GEDI_mid_amazon$agbd) ---------> [1] 2.870049e-10 5.423377e+03
ggplot(agb_forest_age2, aes(x=agbd)) + 
    geom_histogram(binwidth=1) + 
    theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) + 
    scale_x_continuous(limits = c(0, 600)) + 
    scale_y_continuous(limits = c(0, 10000))

# coordinates of a large mature patch:
# -4, -60 (upper left)
# -6.5, -56.5
GEDI_mat <- GEDI_mid_amazon[,]
-4, -60

GEDI <- readRDS("0000000000-0000095232_GEDI.rds")




lulc <- rast('./model_ready_rasters/0000000000-0000095232_lulc_history.tif') 

layer_indices <- nlyr(lulc_brick_masked) - regrowth_mask # year before abandonment
last_LU <- selectRange(lulc_brick_masked, layer_indices)
last_LU <- selectRange(lulc_brick_masked, layer_indices+1)
last_LU <- selectRange(lulc_brick_masked, layer_indices)
last_LU <- selectRange(lulc_brick_masked, layer_indices)
last_LU <- selectRange(lulc_brick_masked, layer_indices)


