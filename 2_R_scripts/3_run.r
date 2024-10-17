library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)
library(gridExtra)
library(scales)

# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_data_utils.r")
source("./2_R_scripts/2_model_utils.r")

set.seed(123)
re_base <- dnorm(1000)

# Import data
data <- import_data("./0_data/data_mapbiomas_4.csv", convert_to_dummy = TRUE)

# Define parameter sets
par_columns_simple <- c("sur_cover")
par_columns_all <- c(
    "sur_cover",
    "cwd",
    "num_fires_before_regrowth",
    colnames(data)[8:19]
)

# Function to run model and calculate R-squared
run_model <- function(model_type, par_columns) {
    if (model_type == "lag") {
        pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
        pars[c("theta", "k0", "m_base", "sd_base", "sd")] <- c(1, 0, 0, 1, 1)

        conditions <- list(
            function(pars) pars["theta"] > 10,
            function(pars) pars["theta"] < 0,
            function(pars) pars["sd"] < 0,
            function(pars) pars["sd_base"] < 0,
            function(pars) pars["m_base"] < 0
        )

        model <- optim(pars, function(pars) likelihood_lag(pars, data, conditions))
        pred <- unname(growth_curve_lag(model$par, data, exp((re_base + model$par["m_base"]) * model$par["sd_base"])))
    } else if (model_type == "B0") {
        pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
        pars[c("B0", "k", "sd")] <- c(0, 0.5, 1)

        conditions <- list(
            function(pars) pars["B0"] < 0,
            function(pars) pars["k"] < 0,
            function(pars) pars["sd"] < 0
        )

        model <- optim(pars, function(pars) likelihood_B0_theta(pars, data, conditions, growth_curve_B0))

        pred <- growth_curve_B0(model$par, data)
    } else if (model_type == "B0_theta") {
        pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
        pars[c("theta", "B0", "sd")] <- c(1, mean(data$agbd), 1)

        conditions <- list(
            function(pars) pars["theta"] > 10,
            function(pars) pars["theta"] < 0,
            function(pars) pars["sd"] < 0
        )

        model <- optim(pars, function(pars) likelihood_B0_theta(pars, data, conditions, growth_curve_B0_theta))
        pred <- growth_curve_B0_theta(model$par, data)
    }

    r_squared <- calc_rsq(data, pred)
    return(list(r_squared = r_squared, pred = pred))
}

# Run models and store results
models <- c("lag", "B0", "B0_theta")
par_sets <- list(par_columns_simple, par_columns_all)
results <- list()

for (model in models) {
    for (par_set in par_sets) {
        key <- paste(model, ifelse(length(par_set) == 1, "simple", "all"), sep = "_")
        results[[key]] <- run_model(model, par_set)
    }
}

# Print R-squared values
for (key in names(results)) {
    cat(sprintf("%s R-squared: %.4f\n", key, results[[key]]$r_squared))
}

# Function to plot histograms using ggplot2
plot_histograms <- function(pred, observed, title) {
    pred_df <- data.frame(value = pred, type = "Predicted")
    obs_df <- data.frame(value = observed, type = "Observed")
    combined_df <- rbind(pred_df, obs_df)

    ggplot(combined_df, aes(x = value, fill = type)) +
        geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
        scale_x_continuous(limits = c(0, 400), breaks = seq(0, 400, 100)) +
        scale_y_continuous(labels = comma) +
        labs(
            title = title,
            x = "Biomass (Mg/ha)",
            y = "Frequency",
            fill = "Data Type"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.title = element_text(face = "bold")
        ) +
        scale_fill_manual(values = c("Predicted" = "#1E90FF", "Observed" = "#FF6347"))
}

# Create a list to store all plots
plot_list <- list()

# Generate plots
for (model in models) {
    for (par_type in c("simple", "all")) {
        key <- paste(model, par_type, sep = "_")
        title <- sprintf("%s Model\n(%s parameters)", model, ifelse(par_type == "simple", "sur_cover only", "all"))
        plot_list[[key]] <- plot_histograms(results[[key]]$pred, data$agbd, title)
    }
}

# Arrange plots in a grid
grid_plot <- grid.arrange(grobs = plot_list, ncol = 2)
