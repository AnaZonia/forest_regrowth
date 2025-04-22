# Data obtained from Global Reforestation Opportunity Assessment (GROA) https://github.com/forc-db/GROA/tree/master



# Load necessary libraries
library(dplyr)
library(terra)
library(tidyverse)


# -----------------------------------------------------
field <- read.csv("0_data/groa_field/biomass_litter_CWD.csv")
sites <- read.csv("0_data/groa_field/sites.csv")

head(sites)

# get the ones that have more precise coordinates
# 1km is the worst precision they get
# buffer around that region


# Merge lat, lon, and site.country from sites to field by matching site.id
field <- field %>%
    filter(site.id %in% sites$site.id) %>%
        left_join(sites %>% dplyr::select(site.id, site.country, site.state, lat_dec, long_dec), by = "site.id")

# Filter the field data for aboveground biomass
field <- subset(field, variables.name == "aboveground_biomass" & site.country == "Brazil")

# Rename and select relevant columns
field <- field %>%
    dplyr::rename(field_biomass = mean_ha, field_age = stand.age) %>%
    dplyr::select(c(field_age, field_biomass, date, lat_dec, long_dec))

# Create a SpatVector (terra object) using lat_dec and long_dec for coordinates
field_spat <- vect(field, geom = c("long_dec", "lat_dec"))

# Set the CRS (assuming WGS84)
crs(field_spat) <- "+proj=longlat +datum=WGS84"

# Export as a shapefile
writeVector(field_spat, "0_data/groa_field/field_biomass.shp", overwrite = TRUE)

# ------------------------------------------------------
# 




# Function to aggregate data based on age intervals
aggregate_biomass <- function(data, age_col, biomass_col, interval = 1) {
    data %>%
        mutate(age_interval = floor({{ age_col }} + 0.5)) %>% # Group into integer intervals
        group_by(age_interval) %>%
        summarise(mean_biomass = mean({{ biomass_col }}, na.rm = TRUE)) %>%
        rename(age = age_interval)
}

# Apply the function to both dataframes
field_aggregated <- aggregate_biomass(field, stand.age, mean_ha)
field_biomass_aggregated <- aggregate_biomass(field_biomass, field_age, field_biom)

# Combine the aggregated data
combined_data <- left_join(field_aggregated, field_biomass_aggregated, by = "age")



# # Calculate mean biomass per field_age
# mean_field_biomass <- field %>%
#     group_by(stand.age) %>%
#     summarise(
#         mean_field_biomass = mean(mean_ha, na.rm = TRUE)
#     )
# mean_field_biomass

field <- read.csv("./0_data/groa_field/field_biomass.csv") %>%
    # remove columns system.index and .geo
    select(-c(system.index, .geo)) %>%
    # make all values < 0 in column data NA
    mutate(across(everything(), ~ ifelse(. < 0, NA, .)))

field



field_age_lulc <- read.csv("./0_data/groa_field/field_age_lulc.csv") %>%
    # remove columns system.index and .geo
    select(-c(system.index, .geo)) %>%
    # keep only those with date > 0
    filter(date > 0)

nrow(field_age_lulc)

# how many rows in field_age_lulc have all column values between the 1985 and 2020 column be exactly the same
field_age_lulc_same <- field_age_lulc %>%
    filter(apply(select(., 1:36), 1, function(x) length(unique(x)) > 1))

head(field_age_lulc_same)


unified_field <- read.csv("./0_data/groa_field/unified_field.csv") %>%
    # remove columns system.index and .geo
    select(-c(system.index, .geo)) %>%
    # make all values < 0 in column data NA
    mutate(across(everything(), ~ ifelse(. < 0, NA, .)))

nrow(unified_field)

unified_field_date <- unified_field %>%
    # keep only those with date > 0
    filter(date > 0)

colnames(unified_field_date)


# identify which ones are in the same site
