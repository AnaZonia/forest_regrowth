# ------------------------------------------------------
# Load libraries
# ------------------------------------------------------
library(tidyverse) # For data wrangling
library(terra)

# ------------------------------------------------------
# Load GROA field and site data
# Source: https://github.com/forc-db/GROA
# ------------------------------------------------------

field <- read.csv("0_data/groa_field/biomass_litter_CWD.csv")
sites <- read.csv("0_data/groa_field/sites.csv")

# ------------------------------------------------------
# Merge site-level metadata (lat/lon/country/state) into field data
# ------------------------------------------------------

sites_unique <- sites %>% distinct(site.id, .keep_all = TRUE)

field <- field %>%
    left_join(
        sites_unique %>% select(site.id, site.country, site.state, lat_dec, long_dec),
        by = "site.id"
    )

# ------------------------------------------------------
# Filter for aboveground biomass measurements in Brazil
# Units are Mg/ha, matching ESA CCI standards
# ------------------------------------------------------

field <- subset(field, variables.name == "aboveground_biomass" & site.country == "Brazil")

# ------------------------------------------------------
# Aggregate biomass data by field age
# This function calculates the mean biomass per field age, grouping by integer intervals.
# The data is used in 2_lag_field_data.r for visualization
# ------------------------------------------------------

aggregated_field <- field %>%
    mutate(age_interval = floor(`stand.age` + 0.5)) %>% # Group into integer intervals
    group_by(age_interval) %>%
    summarise(mean_biomass = mean(`mean_ha`, na.rm = TRUE)) %>%
    rename(age = age_interval)

write.csv(aggregated_field, "0_data/groa_field/average_biomass_per_age.csv", row.names = FALSE)

# ------------------------------------------------------
# Identify and handle plots with repeated measurements
# ------------------------------------------------------

# Identify plot IDs with more than one observation
plot_nums <- field %>%
    group_by(plot.id) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    arrange(desc(n))

# Extract only the rows with repeated plot IDs
field_repeats <- field %>%
    group_by(plot.id) %>%
    filter(n() > 1) %>%
    ungroup()

# there are 16 plots with repeated measurements, totalling 78 non-independent measurements out of 435 total.

# Export repeated measurements for separate analysis
write.csv(field_repeats, "0_data/groa_field/field_repeats.csv", row.names = FALSE)

# From each plot, randomly retain only one measurement (to avoid pseudoreplication)
field <- field %>%
    group_by(plot.id) %>%
    slice_sample(n = 1) %>%
    ungroup()

# ------------------------------------------------------
# Clean up columns and rename for clarity
# ------------------------------------------------------

field <- field %>%
    select(stand.age, mean_ha, date, plot.id, lat_dec, long_dec) %>%
        rename(
            field_biomass = mean_ha,
            field_age = stand.age,
            plot_id = plot.id
        )

plot(field$field_biomass ~ field$field_age, data = field, main = "Field Biomass by Age", xlab = "Field Age (years)", ylab = "Biomass (Mg/ha)")

# ------------------------------------------------------
# Convert to a spatial object (terra::SpatVector)
# and export as a shapefile
# ------------------------------------------------------

# Create spatial object using longitude and latitude
field_spat <- vect(field, geom = c("long_dec", "lat_dec"))

# Assign WGS84 coordinate reference system
crs(field_spat) <- "+proj=longlat +datum=WGS84"

# Export to shapefile
writeVector(field_spat, "0_data/groa_field/field_biomass.shp", overwrite = TRUE)

