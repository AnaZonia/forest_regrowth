library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)
library(gridExtra)
library(scales)

# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_data_utils.r")
source("./2_R_scripts/2_model_utils.r")

# Set seed for reproducibility
set.seed(123)

# Define constants
basic_pars <- c("theta", "k", "m_base", "sd_base", "sd", "B0")
re_base <- dnorm(1000)

# Import data
data <- import_data("./0_data/data_mapbiomas_4.csv", convert_to_dummy = TRUE)
unseen_data <- import_data("./0_data/unseen_data_mapbiomas_4.csv", convert_to_dummy = TRUE)

# Define parameter sets
# par_columns_simple <- c("sur_cover")
# par_columns_all <- setdiff(colnames(data), c("agbd", "age", "nearest_mature_biomass"))
# par_columns_all <- colnames(data)[3:23]

# Import data and convert categorical variables
categorical <- c("ecoreg", "topography")
data <- read_csv("./0_data/data_mapbiomas_4.csv", show_col_types = FALSE)
data[categorical] <- lapply(data[categorical], factor)

# Function to normalize from 0 to 1
normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
# Identify columns to normalize (exclude 'agbd' and categorical variables)
cols_to_normalize <- setdiff(names(data), c("agbd", "nearest_mature_biomass", categorical))
# Normalize the selected columns
data <- data %>%
    mutate(across(all_of(cols_to_normalize), normalize))

# Create dummy variables for present categorical columns
data <- dummy_cols(data,
    select_columns = categorical,
    remove_selected_columns = TRUE
)

# par_columns_all <- names(data)[names(data) != "agbd"]
par_columns_all <- c("age", "sur_cover", "nitro", "num_fires_before_regrowth", "lulc_sum_15")
pars <- setNames(rep(0, length(par_columns_all)), par_columns_all)
pars[c("theta", "B0", "k")] <- c(1, mean(data$agbd), - log(1 - (mean(data$agbd) / mean(data$nearest_mature_biomass))))
pars
conditions <- list(
    function(pars) pars["k"] < 0,
    function(pars) pars["theta"] > 10,
    function(pars) pars["theta"] < 0,
    function(pars) pars["age"] < 0,
    function(pars) pars["age"] > 5
)

model <- optim(pars, function(pars) likelihood_B0_theta(pars, data, conditions, growth_curve_B0_theta, FALSE))
model

pred <- growth_curve_B0_theta(model$par, data)
rsq <- calc_rsq(data, pred)
rsq

all_vars <- names(data)[names(data) != ["agbd", "nearest_mature_biomass"]]
formula <- as.formula(paste("agbd ~", paste(all_vars, collapse = " + ")))
summary(lm(formula, data = data))
summary(lm(agbd ~ age + sur_cover, data = data))












parameter_sets <- list(simple = par_columns_simple, all = par_columns_all)

# Models to evaluate
models <- c("lag", "B0_theta")

run_model <- function(model, pars, data, conditions, NLL) {
    # Define whether function is to be fit through
    # Maximum Likelihood Estimation (TRUE) or
    # Nonlinear Least Squares (FALSE)
    if (NLL) {
        pars["sd"] <- 1
        conditions <- c(conditions, function(pars) pars["sd"] < 0)
    }

    if (model == "lag") {
        pars[c("m_base", "sd_base")] <- c(0, 1)

        conditions <- c(conditions, function(pars) pars["sd_base"] < 0)

        return(cross_valid(data, pars, conditions, likelihood_lag, growth_curve_lag_k, NLL))

    } else if (model == "B0_theta") {
        pars["B0"] <- mean(data$agbd)

        return(cross_valid(data, pars, conditions, likelihood_B0_theta, growth_curve_B0_theta, NLL))

    }
}


results_df <- data.frame()
for (model in models) {
    for (param_set_name in names(parameter_sets)) {

        print(model)
        print(param_set_name)

        pars <- setNames(rep(0.0001, length(parameter_sets[[param_set_name]])), parameter_sets[[param_set_name]])

        pars[c("theta", "k")] <- c(1, -log(1 - (mean(data$agbd) / mean(data$nearest_mature_biomass))))

        conditions <- list(
            function(pars) pars["k"] < 0,
            function(pars) pars["theta"] > 10,
            function(pars) pars["theta"] < 0
        )

        for (nll_flag in c(FALSE, TRUE)) {

            result <- run_model(model, pars, data, conditions, NLL = nll_flag)

            # Create a row for the results
            row <- data.frame(
                model = model,
                param_set_name = param_set_name,
                method = ifelse(nll_flag, "NLL", "NLS"),
                r2_mean = result$r2_mean,
                r2_sd = result$r2_sd,
                r2_unseen = result$r2_unseen
            )

            # Add parameters
            best_pars <- result$best_pars
            for (par_name in names(best_pars)) {
                row[[par_name]] <- best_pars[[par_name]]
            }

            # Append the row to the results dataframe
            results_df <- bind_rows(results_df, row)
        }
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
