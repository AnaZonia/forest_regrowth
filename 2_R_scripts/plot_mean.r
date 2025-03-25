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

pars_init <- find_combination_pars(basic_pars = c("B0", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), norm_data)
pars_init


final_model <- run_optim(norm_data, pars_init, conditions)
final_model$par

baseline_pred <- growth_curve(final_model$par, norm_data)
calc_r2(norm_data, baseline_pred)



for (i in 1:30) {
    data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
        rename(biomass = b1) %>%
        # remove mean_pdsi column
        select(-mean_pdsi)

    # Fit the model on the full data
    norm <- normalize_independently(data)
    norm_data <- norm$train_data
    norm_stats <- norm$train_stats
    # print(norm_stats)

    head(norm_data[, c("age", "num_fires", "sur_cover", "nearest_biomass", "biomass")])

    # use optim GA

    pars_init <- find_combination_pars(basic_pars = c("lag", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), norm_data)

    print(pars_init)

    final_model <- run_optim(norm_data, pars_init, conditions)

    print(final_model$par)
    pred <- growth_curve(final_model$par, norm_data, final_model$par[["lag"]])
    print(calc_r2(norm_data, pred))

}


# check R sauqred values for each

# check whether unnormalized age changes anything wiht B0

# constrain theta




actual_lag_years <- 0.77*64+1
# actual_lag_years

# ages 1 - 35
data[["pred_lag"]] <- growth_curve(final_model$par, norm_data) # with no lag, to give the expected values at low ages

# intermediary ages
data[["pred"]] <- growth_curve(final_model$par, norm_data, final_model$par[["lag"]])


growth_curve_future <- function(pars, data, lag = 0) {

    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    k <- rep(pars[["k0"]], nrow(data))

    age <- data[["age"]]

    if ("lag" %in% names(pars)) {
        pars[["B0"]] <- 0
        age <- age + lag
        age <- age + max(age)
    }

    if (length(non_clim_pars) > 0) {
        k <- (k + rowSums(sapply(non_clim_pars, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE))) * age
    }

    # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
    for (clim_par in intersect(climatic_pars, names(pars))) {
        for (yrs in 1:max(data[["age"]])) {
            indices <- which(data[["age"]] == yrs)
            # Generate a sequence of years for the current age group
            # Starting from 2019 and going back 'yrs' number of years
            last_year <- max(2019 - yrs - round(lag) + 1, 1985)
            year_seq <- seq(2019, last_year, by = -1)
            clim_columns <- paste0(clim_par, "_", year_seq)
            # as.matrix(t()) is used to ensure that rowSums would work in cases with a single row
            k[indices] <- k[indices] + rowSums(as.matrix(t(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]][indices]))))
        }
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    return(pars[["B0"]] + (data[["nearest_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
}

#future forecast
data[["pred_future"]] <- growth_curve_future(final_model$par, norm_data, exp(final_model$par[["lag"]])) # with no lag, to give the expected values at low ages

print(mean(data[["pred_lag"]]))
print(mean(data[["pred"]]))


print(mean(data[["pred_future"]]))







# write.csv(data, paste0(c("0_results/lagged_nolag_unified_data.csv"), collapse = "_"), row.names = FALSE)

# Compute mean values and standard deviations per age
# Restructure the data to have separate columns for each biomass type
mean_biomass_data <- data %>%
    group_by(age) %>%
    summarise(
        mean_pred_lag = median(pred_lag, na.rm = TRUE),
        sd_pred_lag = sd(pred_lag, na.rm = TRUE),
        mean_pred = median(pred, na.rm = TRUE),
        sd_pred = sd(pred, na.rm = TRUE),
        mean_biomass = median(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE),
        mean_future = median(pred_future, na.rm = TRUE),
        sd_future = sd(pred_future, na.rm = TRUE)
    ) %>%
    mutate(age_pred_lag = age) %>%
    mutate(age_pred = age + actual_lag_years) %>%
    mutate(age_future = age + actual_lag_years + 35)
print(mean_biomass_data, n = 35)

# Reshape data to long format
mean_biomass_data_long <- mean_biomass_data %>%
    pivot_longer(
        cols = c(mean_pred_lag, mean_pred, mean_biomass, mean_future),
        names_to = "biomass_type",
        values_to = "mean_value"
    ) %>%
    mutate(
        sd_value = case_when(
            biomass_type == "mean_pred_lag" ~ sd_pred_lag,
            biomass_type == "mean_pred" ~ sd_pred,
            biomass_type == "mean_biomass" ~ sd_biomass,
            biomass_type == "mean_future" ~ sd_future,
            TRUE ~ NA_real_
        ),
        plot_age = case_when(
            biomass_type == "mean_pred_lag" ~ age_pred_lag,
            biomass_type == "mean_pred" ~ age_pred,
            biomass_type == "mean_biomass" ~ age_pred,
            biomass_type == "mean_future" ~ age_future
        )
    )
print(mean_biomass_data_long, n = 35)

mean_biomass_data_long$biomass_type <- factor(mean_biomass_data_long$biomass_type,
    levels = c("mean_pred_lag", "mean_pred", "mean_biomass", "mean_future")
)
colors <- c(
    "mean_pred_lag" = "green",
    "mean_pred" = "red",
    "mean_biomass" = "blue",
    "mean_future" = "purple"
)


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



p <- ggplot(mean_biomass_data_long, aes(x = plot_age, y = mean_value, color = biomass_type)) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = biomass_type),
        alpha = 0.2, color = NA
    ) +
    scale_color_manual(
        values = colors,
        breaks = names(colors), # Ensures the legend follows the correct order
        labels = c("Predicted", "Predicted (No Lag)", "Observed Biomass", "Future")
    ) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(
        title = "Mean Biomass Predictions and Observations",
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)",
        color = "Legend"
    ) +
    geom_point(data = field_aggregated, aes(x = age, y = mean_biomass), color = "purple", size = 2, alpha = 0.7) +
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

