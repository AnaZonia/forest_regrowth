
library(tidyverse)
library(ggplot2)
# library(reshape2)

tst <- read.csv("./0_data/mapbiomas_GEDI.csv")
# tst <- tst %>% filter(biome == 4)
nrow(tst)

# Set up the plotting area to have 1 row and 2 columns
par(mfrow = c(1, 2))

# Plot histogram for GEDI_biomass with restricted x-axis
hist(tst$GEDI_biomass, main = "Histogram of GEDI_biomass", xlab = "GEDI_biomass", xlim = c(0, 500), col = "lightblue", border = "blue")

# Plot histogram for biomass with restricted x-axis
hist(tst$biomass, main = "Histogram of biomass", xlab = "biomass", xlim = c(0, 500), col = "lightgreen", border = "darkgreen")
