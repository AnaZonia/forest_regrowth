# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#           Field Data Analysis and Model Validation
#
#                 Ana Avila - August 2025
#
#     Fit the model to the field data
#     Find theta (shape parameter) value from the field data
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(dplyr)
library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

# Source other scripts
source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_feature_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# Load data
biome <- 1
n_samples <- 10000

# ------------------------------------------------------

# identify which ones are in the same site
field <- read.csv("./0_data/groa_field/field_predictors.csv")
field <- subset(field, biome == 1)

field <- field %>%
    rename(
        biomass = field_biom,
        asymptote = nearest_mature
    ) %>%
    mutate(age = floor(age + 0.5))

# pick one random row per unique value of site_id
field_non_repeats <- field %>%
    group_by(site_id) %>%
    slice_sample(n = 1) %>%
    ungroup()
# there are 44 unique sites in the Amazon (4 of them with one NA value for cec or something of the sort)
# remove rows with NA values in any colum
field <- field[complete.cases(field), ]

# remove columns biome, lat, lon, site_id, plot_id
field_non_repeats <- field_non_repeats %>%
    select(-biome, -lat, -lon, -site_id, -plot_id)

# Convert categorical to factors
field_non_repeats <- field_non_repeats %>%
    mutate(across(any_of(categorical), as.factor))

field_non_repeats <- dummy_cols(field_non_repeats,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

# ------------------------------------------------------


data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# remove the columns in data of the form ecoreg_xxx that have the value xxx not present in field_non_repeats$ecoreg
remove_unused_onehot <- function(data, field_non_repeats, prefix) {
    cols <- grep(paste0("^", prefix, "_"), colnames(data), value = TRUE)
    valid <- paste0(prefix, "_", unique(field_non_repeats[[prefix]]))
    cols_to_remove <- setdiff(cols, valid)
    data %>% select(-all_of(cols_to_remove))
}

data <- data %>%
    remove_unused_onehot(field_non_repeats, "ecoreg") %>%
    remove_unused_onehot(field_non_repeats, "topography")

# Fit the model on the full data
norm_data <- normalize_independently(data)
train_stats <- norm_data$train_stats
norm_data <- norm_data$train_data

pars_init <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(data))[["all_mean_climate"]],
    norm_data
)

final_model <- run_optim(norm_data, pars_init, conditions)

pred <- growth_curve(final_model$par, data = norm_data, lag = final_model$par["lag"])

# Calculate R2
r2 <- calc_r2(norm_data, pred)
r2 # 0.42

# ------------------------------------------------------
# Get R2 when predicting field data
# ------------------------------------------------------


# keep only the variables that are in the model
train_stats <- train_stats %>%
    filter(variable %in% names(final_model$par))

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}

field_non_repeats_scaled <- apply_min_max_scaling(field_non_repeats, train_stats)

pred_field <- growth_curve(final_model$par, data = field_non_repeats_scaled)

# get R2
r2 <- calc_r2(field_non_repeats_scaled, pred_field)
r2





# ------------------------------------------------------
# Get theta from field data
# we remove other predictors (keep only the basic)
# to avoid local minima
# ------------------------------------------------------

field_scaled <- apply_min_max_scaling(field, train_stats)

pars_init <- pars_init[c("k0")]
pars_init$theta <- 1
pars_init$B0 <- 0

final_model <- run_optim(field_scaled, pars_init, conditions)
final_model[['par']]
    
pred_field <- growth_curve(final_model$par, data = field_scaled)

# theta = 2


# ------------------------------------------------------
# Identify and handle plots with repeated measurements
# ------------------------------------------------------

# Identify plot IDs with more than one observation
plot_nums <- field %>%
    group_by(plot_id) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    arrange(desc(n))

# take only rows with more than 5 observations for the same plot_id
field_repeats <- field %>%
    group_by(plot_id) %>%
    filter(n() > 2) %>%
    ungroup()

plot(field_repeats$age, field_repeats$biomass,
     xlab = "Age", ylab = "Biomass",
     main = "Repeated Measurements for Plots with More than 5 Observations"
)


pred_repeats
table(field_repeats$site_id)
pred_repeats <- growth_curve(final_model$par, data = field_repeats)
plot(field_repeats$age, field_repeats$biomass, xlab = "Age", ylab = "Biomass", main = "Predictions for site with repeated measurements")
points(field_repeats$age, pred_repeats, col = "red", pch = 19)






# ------------------------------------------------------
# Make extended data figures
# ------------------------------------------------------



png("./0_results/figures/extended/field_predictions_scatterplot.png", width = 800, height = 600)
plot(field_non_repeats$age, field_non_repeats$biomass, xlab = "Field Age", ylab = "Field Biomass", main = "Age vs Biomass")
points(field_non_repeats$age, pred_field, col = "red", pch = 19)
dev.off()


png("./0_results/figures/extended/predicted_vs_observed_field.png", width = 800, height = 600)
plot(pred_field, field_non_repeats$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)
dev.off()



png("./0_results/figures/extended/field_age_histogram.png", width = 1800, height = 1400, res = 300)

ggplot(field, aes(x = age)) +
    geom_histogram(
        binwidth = 5,
        fill = "grey30",
        color = "white",
        boundary = 0
    ) +
    labs(
        x = "Forest age (years)",
        y = "Number of plots"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 20),
        axis.title = element_text(face = "bold", size = 22),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10)
    )

dev.off()
