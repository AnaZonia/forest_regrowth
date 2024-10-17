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
NLS <- TRUE # fit with Nonlinear Least Squares - IF FALSE, fit with Maximum Likelihood
basic_pars <- c("theta", "k", "m_base", "sd_base", "sd", "B0")

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

# Run models and store results
models <- c("lag", "B0_theta")
par_sets <- list(par_columns_simple, par_columns_all)
parameter_sets <- list(simple = par_columns_simple, all = par_columns_all)


results <- list()
for (model in models) {
    for (param_set_name in names(parameter_sets)) {

        if (model == "lag") {
            pars <- setNames(rep(0.0001, length(parameter_sets[[param_set_name]])), parameter_sets[[param_set_name]])
            pars[c("theta", "k", "m_base", "sd_base", "sd")] <- c(1, 0, 0, 1, 1)

            conditions <- list(
                function(pars) pars["theta"] > 10,
                function(pars) pars["theta"] < 0,
                function(pars) pars["sd"] < 0,
                function(pars) pars["sd_base"] < 0,
                function(pars) pars["m_base"] < 0
            )

            result <- cross_valid(data, pars, conditions, likelihood_lag, growth_curve_lag_k)

        } else if (model == "B0_theta") {

            pars <- setNames(rep(0.0001, length(parameter_sets[[param_set_name]])), parameter_sets[[param_set_name]])
            pars[c("theta", "B0", "sd")] <- c(1, mean(data$agbd), 1)

            conditions <- list(
                function(pars) pars["theta"] > 10,
                function(pars) pars["theta"] < 0,
                function(pars) pars["sd"] < 0
            )

            result <- cross_valid(data, pars, conditions, likelihood_B0_theta, growth_curve_B0_theta)
        }

        results[[paste(model, param_set_name, sep = "_")]] <- result
    }
}



# compare_results <- function(results) {
#     r_squared_values <- sapply(results, function(x) x$rsq)
#     r_squared_df <- data.frame(
#         Model = names(r_squared_values),
#         R_squared = unname(r_squared_values)
#     )
#     print(r_squared_df)
# }

# # Function to plot histograms using ggplot2
# plot_histograms <- function(pred, observed, title) {
#     pred_df <- data.frame(value = pred, type = "Predicted")
#     obs_df <- data.frame(value = observed, type = "Observed")
#     combined_df <- rbind(pred_df, obs_df)

#     ggplot(combined_df, aes(x = value, fill = type)) +
#         geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
#         scale_x_continuous(limits = c(0, 400), breaks = seq(0, 400, 100)) +
#         scale_y_continuous(labels = comma) +
#         labs(
#             title = title,
#             x = "Biomass (Mg/ha)",
#             y = "Frequency",
#             fill = "Data Type"
#         ) +
#         theme_minimal() +
#         theme(
#             legend.position = "bottom",
#             plot.title = element_text(hjust = 0.5, face = "bold"),
#             axis.title = element_text(face = "bold"),
#             legend.title = element_text(face = "bold")
#         ) +
#         scale_fill_manual(values = c("Predicted" = "#1E90FF", "Observed" = "#FF6347"))
# }

# # Create a list to store all plots
# plot_list <- list()

# # Generate plots
# for (model in models) {
#     for (par_type in c("simple", "all")) {
#         key <- paste(model, par_type, sep = "_")
#         title <- sprintf("%s Model\n(%s parameters)", model, ifelse(par_type == "simple", "sur_cover only", "all"))
#         plot_list[[key]] <- plot_histograms(results[[key]]$pred, data$agbd, title)
#     }
# }

# # Arrange plots in a grid
# grid_plot <- grid.arrange(grobs = plot_list, ncol = 2)
