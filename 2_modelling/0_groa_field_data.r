# makes shapefile of GROA field data for visualization
# and for keeping only the sites in the Amazon biome (subsequently done in Google Earth Engine)


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

field <- field %>%
    select(stand.age, mean_ha, site.id, plot.id, lat_dec, long_dec) %>%
    rename(
        field_biomass = mean_ha,
        age = stand.age,
        site_id = site.id,
        plot_id = plot.id
    )

# ------------------------------------------------------
# Convert to a spatial object (terra::SpatVector)
# and export as a shapefile
# ------------------------------------------------------

# Create spatial object using longitude and latitude
field_spat <- vect(field, geom = c("long_dec", "lat_dec"))

# Assign WGS84 coordinate reference system
crs(field_spat) <- "+proj=longlat +datum=WGS84"

# Export to shapefile
writeVector(field_spat, "0_data/groa_field/field.shp", overwrite = TRUE)


