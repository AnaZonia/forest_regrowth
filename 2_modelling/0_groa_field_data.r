# Data obtained from Global Reforestation Opportunity Assessment (GROA) https://github.com/forc-db/GROA/tree/master



# Load necessary libraries
library(dplyr)
library(terra)
library(tidyverse)


# -----------------------------------------------------
field <- read.csv("0_data/groa_field/biomass_litter_CWD.csv")
sites <- read.csv("0_data/groa_field/sites.csv")

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
