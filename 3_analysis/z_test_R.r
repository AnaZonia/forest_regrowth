library(foreach)
library(doParallel)
library(tidyverse)
library(ggplot2)
library(cowplot) # For legend extraction
library(ggpubr) # For legend extraction

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)
# For plotting
options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 20))


# folder <- "./0_data/grid_10k_amazon_secondary_old"

# # Get all CSV files in the folder
# files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

# for (f in files) {
#     df <- read.csv(f, check.names = FALSE)

#     # If the column "secondary_area" exists, rename it to "area"
#     if ("distance_deep_forest" %in% names(df)) {
#         names(df)[names(df) == "distance_deep_forest"] <- "dist"
#     }

#     # Overwrite the original CSV (no row names)
#     write.csv(df, f, row.names = FALSE)
# }





data <- import_data("grid_10k_amazon_secondary_old", biome = 1, n_samples = 150000)

indices <- sample(c(1:5), nrow(data), replace = TRUE)

index <- 3
data <- data[indices == index, ]

norm_data <- normalize_independently(data)
train_stats <- norm_data$train_stats
norm_data <- norm_data$train_data

pars_init <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = setdiff(data_pars_options(colnames(data))[["all"]], "floodable_forests"),
    data = norm_data
)

model <- run_optim(norm_data, pars_init[[1]], conditions)

model

# index = 1, lag = 25.63, distance_deep_forest = -0.011, mean(pred) = 184.7
# index = 2, lag = 26.5, distance_deep_forest = -0.0087, mean(pred) = 184.6
# r$> head(pred)
# [1] 238.8791 230.1101 224.6241 223.6540 230.9031 237.7478

# --------------------------------------#

data_1k <- import_data(paste0("grid_1k_amazon_pastureland"), biome = 1, n_samples = "all")

data_1k <- data_1k$df

data_1k <- apply_min_max_scaling(data_1k, train_stats)

data_1k$age <- 30

pred <- growth_curve(model$par, data_1k, model$par["lag"])
mean(pred)

