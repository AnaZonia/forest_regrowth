# Data obtained from Global Reforestation Opportunity Assessment (GROA) https://github.com/forc-db/GROA/tree/master



# Load necessary libraries
library(dplyr)
library(terra)
library(tidyverse)

setwd("~/Documents/forest_regrowth")

# Read the CSV files
field <- read.csv("0_data/groa_field/biomass_litter_CWD.csv")

sites <- read.csv("0_data/groa_field/sites.csv")

# Filter the field data for aboveground biomass
aboveground_biomass <- subset(field, variables.name == "aboveground_biomass")

# Merge lat, lon, and site.country from sites to field by matching site.id
field <- field %>%
    left_join(sites %>% dplyr::select(site.id, lat_dec, long_dec, site.country), by = "site.id")

# Keep only rows where site.country is Brazil
field_brazil <- subset(field, site.country == "Brazil")

# Step 1: Keep only entries with mean_ha less than 1000
field_brazil <- subset(field_brazil, mean_ha < 1000)

# Step 2: Remove entries where stand.age is 1 and mean_ha is greater than 100
field_brazil <- subset(field_brazil, !(stand.age == 1 & mean_ha > 100))

field_brazil <- field_brazil %>%
    dplyr::rename(field_biomass = mean_ha, field_age = stand.age) %>%
    dplyr::select(c(field_age, field_biomass, lat_dec, long_dec))


# Convert the data to a SpatialPointsDataFrame
coordinates(field_brazil) <- ~ long_dec + lat_dec # Assuming lat_dec and long_dec are the columns for coordinates

# Set the CRS (assuming WGS84)
proj4string(field_brazil) <- CRS("+proj=longlat +datum=WGS84")

# Export as a shapefile
shapefile(field_brazil, "~/Documents/data/field_biomass.shp", overwrite = TRUE)


# -----------------------------------------------------
# Load necessary libraries

# Load your data
field <- read.csv("~/Documents/data/biomass_litter_CWD.csv")
field <- subset(field, variables.name == "aboveground_biomass")
field$field_biom <- field$mean_ha * 0.5

# Load the second dataset
field_biomass <- read.csv("~/Documents/data/field_biomass_with_biome.csv")
field_biomass$field_biom <- field_biomass$field_biom * 0.5

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
combined_data



# # Calculate mean biomass per field_age
# mean_field_biomass <- field_biomass %>%
#     group_by(field_age) %>%
#     summarise(
#         mean_field_biomass = mean(field_biom, na.rm = TRUE)
#     )
# mean_field_biomass

# # Calculate mean biomass per field_age
# mean_field_biomass <- field %>%
#     group_by(stand.age) %>%
#     summarise(
#         mean_field_biomass = mean(mean_ha, na.rm = TRUE)
#     )
# mean_field_biomass


# Calculate mean biomass and mean heinrich_biomass per mapbiomas_age
mean_biomass_neotropics <- field %>%
    group_by(stand.age) %>%
    summarise(
        mean_biomass_neotropics = mean(mean_ha, na.rm = TRUE)
    )

# Calculate mean biomass per field_age
mean_biomass_amazon <- field_biomass %>%
    group_by(field_age) %>%
    summarise(
        mean_biomass_amazon = mean(field_biom, na.rm = TRUE)
    )

# Combine the datasets for plotting
combined_data <- bind_rows(
    mean_biomass_neotropics %>%
        rename(age = stand.age, biomass_type = "mean_biomass_neotropics") %>%
        mutate(value = mean_biomass_neotropics),
    mean_biomass_amazon %>%
        rename(age = field_age, biomass_type = "mean_biomass_amazon") %>%
        mutate(value = mean_biomass_amazon)
)

# Create the combined plot
ggplot(combined_data) +
    geom_line(aes(x = age, y = value, color = biomass_type), size = 1) +
    labs(
        title = "Mean Biomass in Neotropics and Amazon",
        x = "Age",
        y = "Biomass",
        color = "Biomass Type"
    ) +
    scale_color_manual(values = c(
        "mean_biomass_neotropics" = "blue",
        "mean_biomass_amazon" = "red"
    )) +
    theme_minimal()
