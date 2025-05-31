# ------------------------------------------------- #
# Figure - Model Performance Section
# Plot the importance of each category of variables
# ------------------------------------------------- #



tst <- tst[tst$importance_pct > 0.5, ]
# tst <- tst[tst$variable != "age", ]
tst

# Define separate lists for each category
landscape_vars <- c("sur_cover", "dist", "age")
climate_vars <- c("mean_srad", "mean_def", "mean_vpd", "mean_aet", "mean_pr", "mean_temp", "mean_soil")
disturbance_vars <- c("num_fires")
vegetation_vars <- c("floodable_forests")
protected_vars <- c("protec", "indig")
soil_vars <- c("phh2o", "sand", "clay", "soc", "ocs", "ocd", "cfvo", "nitro", "cec")

# Create a column to store the category for each variable
tst <- tst %>%
    mutate(category = case_when(
        variable %in% landscape_vars ~ "Landscape",
        variable %in% climate_vars ~ "Climate",
        variable %in% disturbance_vars ~ "Disturbance",
        variable %in% vegetation_vars ~ "Vegetation",
        variable %in% protected_vars ~ "Protected Area",
        variable %in% soil_vars ~ "Soil",
        TRUE ~ "Other" # Fallback if variable doesn't match
    ))

# Define custom colors for each category
category_colors <- c(
    "Landscape" = "#F0E442",
    "Climate" = "#0072B2",
    "Disturbance" = "#CC79A7",
    "Vegetation" = "#009E73",
    "Protected Area" = "#E69F00",
    "Soil" = "#D55E00"
)

# Create a mapping of short variable names to their full names
variable_names <- c(
    age = "Age",
    sur_cover = "Surrounding Mature Forest Cover",
    mean_srad = "Mean Solar Radiation",
    mean_def = "Mean Climatic Water Deficit",
    num_fires = "Number of Fires",
    phh2o = "Soil pH",
    mean_vpd = "Mean Vapor Pressure Deficit",
    mean_aet = "Mean Actual Evapotranspiration",
    floodable_forests = "Floodable Forests",
    sand = "Sand Content",
    mean_soil = "Mean Soil Moisture",
    protec = "Protected Area",
    ocs = "Organic Carbon Stock",
    ocd = "Organic Carbon Density",
    cfvo = "Coarse Fragments Volume",
    nitro = "Soil Nitrogen",
    dist = "Distance",
    indig = "Indigenous Area",
    cec = "Cation Exchange Capacity",
    clay = "Clay Content",
    mean_pr = "Mean Precipitation",
    mean_temp = "Mean Temperature",
    soc = "Soil Organic Carbon"
)

# Add full names to the importance dataframe
tst$full_name <- variable_names[tst$variable]
