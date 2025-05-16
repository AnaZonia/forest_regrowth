# Data obtained from Global Reforestation Opportunity Assessment (GROA) https://github.com/forc-db/GROA/tree/master



# -----------------------------------------------------
field <- read.csv("0_data/groa_field/biomass_litter_CWD.csv")
sites <- read.csv("0_data/groa_field/sites.csv")

# get the ones that have more precise coordinates
# 1km is the worst precision they get
# buffer around that region


# Merge lat, lon, and site.country from sites to field by matching site.id
field <- field %>%
    filter(site.id %in% sites$site.id) %>%
        left_join(sites %>% dplyr::select(site.id, site.country, site.state, lat_dec, long_dec), by = "site.id")

# Filter the field data for aboveground biomass
field <- subset(field, variables.name == "aboveground_biomass" & site.country == "Brazil")

plot_nums <- field %>%
    group_by(plot.id) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    arrange(desc(n))
# table(plot_nums$n)

# Rename and select relevant columns
field <- field %>%
    dplyr::rename(field_biomass = mean_ha, field_age = stand.age) %>%
    dplyr::select(c(field_age, field_biomass, date, lat_dec, long_dec))

head(field)

# Create a SpatVector (terra object) using lat_dec and long_dec for coordinates
field_spat <- vect(field, geom = c("long_dec", "lat_dec"))

# Set the CRS (assuming WGS84)
crs(field_spat) <- "+proj=longlat +datum=WGS84"

# # Export as a shapefile
# writeVector(field_spat, "0_data/groa_field/field_biomass.shp", overwrite = TRUE)

# ------------------------------------------------------
# 

# look into plots with repeated measurements

# select only rows with plot.id that show up more than once
field_repeats <- field %>%
    group_by(plot.id) %>%
    filter(n() > 1) %>%
    ungroup()
head(field_repeats)


# --------------------------



# Calculate mean biomass per field_age
aggregate_biomass <- function(data, age_col, biomass_col, interval = 1) {
    data %>%
        mutate(age_interval = floor({{ age_col }} + 0.5)) %>% # Group into integer intervals
        group_by(age_interval) %>%
        summarise(mean_biomass = mean({{ biomass_col }}, na.rm = TRUE)) %>%
        rename(age = age_interval)
}

# Apply the function to both dataframes
field_aggregated <- aggregate_biomass(field, field_age, field_biomass)

#plot field raw data scatterplot and averages from field_aggregated overlapping
plot(field$field_age, field$field_biomass, xlab = "Field Age", ylab = "Field Biomass", main = "Field Age vs Field Biomass")
points(field_aggregated$age, field_aggregated$mean_biomass, col = "red", pch = 19)


# --------------------------
