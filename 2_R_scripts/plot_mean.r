# Load necessary libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# get data for plotting
source("2_R_scripts/1_modelling.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_parameters.r")

biome = 1
n_samples = 10000
# Load data
# data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)
set.seed(1)
registerDoParallel(cores = 4)

data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
        rename(biomass = b1) %>%
        # remove mean_pdsi column
        select(-c(mean_pdsi))

norm <- normalize_independently(data)
norm_data <- norm$train_data
norm_stats <- norm$train_stats

pars_init <- find_combination_pars(basic_pars = c("lag", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), norm_data)
pars_init


final_model <- run_optim(norm_data, pars_init, conditions)
final_model$par

baseline_pred <- growth_curve(final_model$par, norm_data)
calc_r2(norm_data, baseline_pred)



# for (i in 1:30) {
#     data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
#         rename(biomass = b1) %>%
#         # remove mean_pdsi column
#         select(-mean_pdsi)

#     # Fit the model on the full data
#     norm <- normalize_independently(data)
#     norm_data <- norm$train_data
#     norm_stats <- norm$train_stats
#     # print(norm_stats)

#     head(norm_data[, c("age", "num_fires", "sur_cover", "nearest_biomass", "biomass")])

#     # use optim GA

#     pars_init <- find_combination_pars(basic_pars = c("lag", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), norm_data)

#     print(pars_init)

#     final_model <- run_optim(norm_data, pars_init, conditions)

#     print(final_model$par)
#     pred <- growth_curve(final_model$par, norm_data, final_model$par[["lag"]])
#     print(calc_r2(norm_data, pred))

# }


# check R sauqred values for each

# check whether unnormalized age changes anything wiht B0

# constrain theta

# ages 1 - 35
data[["pred_lag"]] <- growth_curve(final_model$par, norm_data) # with no lag, to give the expected values at low ages

# intermediary ages
data[["pred"]] <- growth_curve(final_model$par, norm_data, final_model$par[["lag"]])



mean(data[["pred_lag"]])
mean(data[["pred"]])

# --------------------------------- Plotting ---------------------------------#

# Load the second dataset and modify as needed
field_biomass <- read.csv("0_data/groa_field/field_biomass_with_biome.csv")
field_biomass <- subset(field_biomass, biome == 1) # Filter for specific biome

# field_biomass$field_biom <- field_biomass$field_biom * 0.5

# Aggregate field biomass data by age
aggregate_biomass <- function(data, age_col, biomass_col, interval = 1) {
    data %>%
        mutate(age_interval = floor({{ age_col }} + 0.5)) %>% # Group into integer intervals
        group_by(age_interval) %>%
        summarise(mean_biomass = mean({{ biomass_col }} * (0.5), na.rm = TRUE)) %>%
        rename(age = age_interval)
}

field_aggregated <- aggregate_biomass(field_biomass, field_age, field_biom)

# ---------------------------- Plotting ----------------------------


# write.csv(data, paste0(c("0_results/lagged_nolag_unified_data.csv"), collapse = "_"), row.names = FALSE)

# Compute mean values and standard deviations per age
# Restructure the data to have separate columns for each biomass type
mean_biomass_data <- data %>%
    group_by(age) %>%
    summarise(
        mean_pred = median(pred, na.rm = TRUE),
        sd_pred = sd(pred, na.rm = TRUE),
        mean_biomass = median(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    ) %>%
    mutate(age = age + 34)

pred_data <- data %>%
    group_by(age) %>%
    summarise(
        mean_pred_lag = median(pred_lag, na.rm = TRUE),
        sd_pred_lag = sd(pred_lag, na.rm = TRUE),
    )

# Merge the data frames based on age using full_join
pred_data <- full_join(pred_data, mean_biomass_data, by = "age")

zero_intercept_pred <- mean_biomass_data$mean_pred - min(mean_biomass_data$mean_pred)
zero_intercept_pred <- c(zero_intercept_pred, zero_intercept_pred + max(zero_intercept_pred))

pred_data$zero_intercept_pred <- zero_intercept_pred[-1]
head(pred_data)
# Define colors

colors <- c(
    "mean_pred_lag" = "blue",
    "mean_pred" = "purple",
    "mean_biomass" = "red",
    "zero_intercept_pred" = "red"
)

# Create the plot
p <- ggplot(pred_data, aes(x = age)) +
    geom_line(aes(y = mean_pred_lag, color = "mean_pred_lag"), size = 1.5, linetype = "dashed") +
    geom_line(aes(y = mean_pred, color = "mean_pred"), size = 1.5, na.rm = TRUE) +
    geom_line(aes(y = mean_biomass, color = "mean_biomass"), size = 1.5, na.rm = TRUE) +
    geom_line(aes(y = zero_intercept_pred, color = "zero_intercept_pred"), size = 1.5) +
    geom_ribbon(aes(ymin = mean_pred_lag - sd_pred_lag, ymax = mean_pred_lag + sd_pred_lag, fill = "mean_pred_lag"),
        alpha = 0.2, color = NA
    ) +
    geom_ribbon(aes(ymin = mean_pred - sd_pred, ymax = mean_pred + sd_pred, fill = "mean_pred"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    geom_ribbon(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass, fill = "mean_biomass"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    scale_color_manual(values = colors, name = "Legend") + # Map colors correctly
    scale_fill_manual(values = colors, guide = "none") +
    labs(
        title = "Mean Biomass Predictions and Observations",
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)"
    ) +
    theme_minimal(base_size = 20) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = c(0.25, 0.75),
        legend.background = element_rect(fill = "white", color = "black", size = 1),
        aspect.ratio = 1 / 2
    )
p
