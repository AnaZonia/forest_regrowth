# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Source other scripts
source("2_R_scripts/1_parameters.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_modelling.r")

# Get configuration
n_samples <- 10000

# Set up parallel processing
set.seed(1)
registerDoParallel(cores = 10)

# Load data
data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Function to run a single experiment
run_experiment <- function(basic_pars_name, data_pars_name, biome) {

    # Get parameters
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

    # Run cross-validation
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    # Return summary
    return(data.frame(
        basic_pars_name = basic_pars_name,
        data_pars_name = data_pars_name,
        biome = biome,
        mean_r2 = mean(cv_results)
    ))
}

results <- data.frame()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------- Check Importance of parameters included ---------------------------------#
# biome = Amazon
# basic_pars = intercept
# just vary data_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for (name in names(data_pars_options(colnames(data)))) {
    result <- run_experiment("intercept", name, 1)
    print(result)
    results <- rbind(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------- Check Model Forms ---------------------------------#
# biome = Amazon
# data_pars = all_mean_clim
# just vary basic_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for (name in names(basic_pars_options)) {
    result <- run_experiment(name, "all_mean_climate", 1)
    print(result)
    results <- rbind(result)
}

# write.csv(results, "0_results/model_comparisons.csv")


# get data for plotting
source("2_R_scripts/1_modelling.r")
# Load data
data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data
pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = data_pars_options(colnames(data))[["land_use_landscape_only"]], norm_data)

final_model <- run_optim(norm_data, pars_init, conditions)


data[["pred_lagged"]] <- growth_curve(final_model$par, norm_data, exp(final_model$par[["lag"]]))
data[["pred_nolag"]] <- growth_curve(final_model$par, norm_data) # with no lag, to give the expected values at low ages
data[["pred_future"]] <- growth_curve(final_model$par, norm_data, exp(final_model$par[["lag"]])+35) # with no lag, to give the expected values at low ages

print(mean(data[["pred_lagged"]]))
print(mean(data[["pred_nolag"]]))
print(mean(data[["pred_future"]]))



write.csv(data, paste0(c("0_results/lagged_nolag_unified_data.csv"), collapse = "_"), row.names = FALSE)



# --------------------------------- Plotting ---------------------------------#

# get data for plotting
source("2_R_scripts/1_modelling.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_parameters.r")

biome <- 1
n_samples <- 10000
# Load data

data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
    rename(biomass = b1) %>%
    # remove mean_pdsi column
    select(-c(mean_pdsi))

options <- data_pars_options(colnames(data))
options


norm <- normalize_independently(data)
norm_data <- norm$train_data
norm_stats <- norm$train_stats

pars_init <- find_combination_pars(basic_pars = c("B0", "k0", "theta"), "data_pars" = options[["all_no_categorical"]], norm_data)
pars_init


final_model <- run_optim(norm_data, pars_init, conditions)
final_model$par



permutation_importance <- calculate_permutation_importance(final_model, norm_data, options[["all_no_categorical"]])

library(ggplot2)


tst <- readRDS("tst.rds")

# ------------------------------------------------- #
# Plot the variable importance scores
# ------------------------------------------------- #
tst <- tst[tst$importance_pct > 0.5, ]
# tst <- tst[tst$variable != "age", ]
tst

# Define separate lists for each category
landscape_vars <- c("sur_cover", "dist", "age")
climate_vars <- c("mean_srad", "mean_def", "mean_vpd", "mean_aet", "mean_pr", "mean_temp")
disturbance_vars <- c("num_fires")
vegetation_vars <- c("floodable_forests")
protected_vars <- c("protec", "indig")
soil_vars <- c("phh2o", "sand", "clay", "soc", "ocs", "ocd", "cfvo", "nitro", "cec", "mean_soil")

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

# Create the plot with color-coded categories
importance_plot <- ggplot(tst, aes(x = reorder(full_name, importance), y = importance, fill = category)) +
    geom_col(width = 0.9) +  # Reduce bar spacing
    scale_fill_manual(values = category_colors) + # Apply custom colors
    coord_flip() +
    labs(
        title = "Permutation Importance Scores",
        x = "Variables",
        y = "Importance Score",
        fill = "Category"
    ) +
    theme_minimal(base_size = 16) + # Increases overall text size
    theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # Big, bold title, centered
        axis.title.x = element_text(size = 18, face = "bold"), # Bigger X-axis title
        axis.title.y = element_text(size = 18, face = "bold"), # Bigger Y-axis title
        axis.text.x = element_text(size = 14), # Increase X-axis text size
        axis.text.y = element_text(size = 14), # Increase Y-axis text size
        legend.title = element_text(size = 16, face = "bold"), # Bigger legend title
        legend.text = element_text(size = 14), # Bigger legend text
        panel.grid.major.y = element_blank(),  # Remove horizontal gridlines
        panel.grid.minor.y = element_blank()
    )


# Print the plot
# print(importance_plot)


# Optionally, save the plot to a file
ggsave("variable_importance_plot.png", plot = importance_plot, width = 10, height = 7)




analyze_variable_importance(final_model, options[["all_no_categorical"]])




library(ggplot2)

# Create a simple dataframe
test_data <- data.frame(
  x = 1:10,
  y = rnorm(10)
)

# Create a simple plot
test_plot <- ggplot(test_data, aes(x = x, y = y)) +
  geom_point() +
  labs(
    title = "Test Plot",
    x = "X Axis",
    y = "Y Axis"
  ) +
  theme_minimal()

# Print the plot
print(test_plot)
coeffs

tst

saveRDS(tst, "tst.rds")
