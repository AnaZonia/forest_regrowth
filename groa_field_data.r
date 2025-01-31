# Load necessary libraries
library(dplyr)
library(raster)
library(tidyverse)
library(sp)

# Read the CSV files
field <- read.csv("~/Documents/data/biomass_litter_CWD.csv")
sites <- read.csv("~/Documents/data/sites.csv")

# Filter the field data for aboveground biomass
field <- subset(field, variables.name == "aboveground_biomass")

# Merge lat, lon, and site.country from sites to field by matching site.id
field <- field %>%
    left_join(sites %>% dplyr::select(site.id, lat_dec, long_dec, site.country), by = "site.id")

# Keep only rows where site.country is Brazil
field_brazil <- subset(field, site.country == "Brazil")

# Step 1: Keep only entries with mean_ha less than 1000
field_brazil <- subset(field_brazil, mean_ha < 1000)

# Step 2: Remove entries where stand.age is 1 and mean_ha is greater than 100
field_brazil <- subset(field_brazil, !(stand.age == 1 & mean_ha > 100))

# Convert the data to a SpatialPointsDataFrame
coordinates(field_brazil) <- ~ long_dec + lat_dec # Assuming lat_dec and long_dec are the columns for coordinates

# Set the CRS (assuming WGS84)
proj4string(field_brazil) <- CRS("+proj=longlat +datum=WGS84")

# Export as a shapefile
shapefile(field_brazil, "~/Documents/data/field_brazil_filtered.shp")
