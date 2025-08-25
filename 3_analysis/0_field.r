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

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------- Field Data Cleaning ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data <- read.csv("./0_data/groa_field/field_predictors.csv")
field_data <- subset(field_data, biome == 1)
field_data <- field_data %>%
    rename(
        biomass = field_biom,
        asymptote = nearest_mature
    ) %>%
    mutate(age = floor(age + 0.5)) %>%
    drop_na()

# # Assuming field_data already cleaned as shown earlier
# field_summary <- field_data %>%
#     group_by(site_id, age) %>%
#     summarise(
#         biomass = median(biomass, na.rm = TRUE),
#         .groups = "drop"
#     )

site_env <- field_data %>%
    select(-age, -biomass, -plot_id) %>%
    distinct(site_id, .keep_all = TRUE)

# field_summary <- field_summary %>%
    # left_join(site_env, by = "site_id")


field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Train Model on Satellite Data -------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

satellite_data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# remove the columns in data of the form ecoreg_xxx that have the value xxx not present in field_non_repeats$ecoreg
remove_unused_dummies <- function(data, field_non_repeats, prefix) {
    cols <- grep(paste0("^", prefix, "_"), colnames(data), value = TRUE)
    valid <- paste0(prefix, "_", unique(field_non_repeats[[prefix]]))
    cols_to_remove <- setdiff(cols, valid)
    data %>% select(-all_of(cols_to_remove))
}

satellite_data <- satellite_data %>%
    remove_unused_dummies(field_data, "ecoreg") %>%
    remove_unused_dummies(field_data, "topography")

sat_data_scaled <- normalize_independently(satellite_data)
train_stats <- sat_data_scaled$train_stats
sat_data_scaled <- sat_data_scaled$train_data

tst <- subset(sat_data_scaled, age == 1)
mean(tst$biomass)


# n = 20
# # sample 100 rows per age
# sat_data_scaled <- sat_data_scaled %>%
#     group_by(age) %>%
#     slice_sample(n = n) %>%
#     ungroup()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Get theta from field data --------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# - We use only age to find the shape parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data_scaled <- apply_min_max_scaling(field_data, train_stats)
common_cols <- intersect(colnames(field_data_scaled), colnames(sat_data_scaled))
field_data_scaled <- field_data_scaled[, common_cols]
sat_data_scaled <- sat_data_scaled[, common_cols]
sat_data_scaled$satellite <- 1
field_data_scaled$satellite <- 0
unite_field_satellite <- rbind(field_data_scaled, sat_data_scaled)

unite_field_satellite <- field_data_scaled

# field_data_scaled <- field_data_scaled %>%
#     select(biomass, asymptote, age) %>%
#     mutate(satellite = 0)

# sat_data_scaled2 <- sat_data_scaled %>%
#     select(biomass, age) %>%
#     mutate(satellite = 1, asymptote = 210)
# unite_field_satellite2 <- rbind(field_data_scaled2, sat_data_scaled2)

# keep only rows with biomass < 300
unite_field_satellite <- unite_field_satellite %>%
    filter(biomass < 400)


init_params <- c(
    "k0" = 0.01,  # Initial guess for k0
    "theta" = 2  # Initial guess for theta
    # "lag" = 2  # Initial guess for lag
)

# how variable is AIC if I force theta at 1.75 with 100 data points?
# fit theta and force data = 1.75 for the same dataset (sample it 3x)

# and get R2

# source("2_modelling/2_modelling.r")

# init_params <- c(
#     "k0" = 0.01, # Initial guess for k0
#     "lag" = 2.5 # Initial guess for lag
# )

# source("2_modelling/2_modelling.r")


model <- run_optim(unite_field_satellite, init_params, conditions)
model

unite_field_satellite1 <- mutate(unite_field_satellite, age = ifelse(satellite == 1, age + model[["par"]][["lag"]], age))

pred <- growth_curve(model[["par"]], data = unite_field_satellite, lag = model[["par"]][["lag"]])

r2 <- calc_r2(unite_field_satellite1, pred)

df <- unite_field_satellite1
df$pred <- pred


label_text <- paste0(
    # "lag = ", model[["par"]][["lag"]],
    ", theta = ", model[["par"]][["theta"]],
    ", r2 = ", r2,
    ", ", n, " per age"
)

plot <- ggplot(df, aes(x = age, y = biomass)) +
    geom_point() +
    geom_point(aes(y = pred), color = "red") +
    annotate("text", x = Inf, y = Inf, label = label_text, hjust = 1.1, vjust = 2, size = 5) +
    theme_minimal()

plot

# Save to file
ggsave(paste0("./0_results/", n, ".png"),
    plot = plot,
    bg = "white"
)









# CSV with the R2 and the theta value
write.csv(data.frame(
    theta = theta,
    lag = lag
), file = "./0_results/theta_lag.csv", row.names = FALSE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(sat_data_scaled))[["all_mean_climate"]],
    sat_data_scaled
)

sat_fit_result <- run_optim(sat_data_scaled, init_params, conditions)
sat_pred_biomass <- growth_curve(sat_fit_result$par, data = sat_data_scaled, lag = sat_fit_result$par["lag"])
r2_satellite <- calc_r2(sat_data_scaled, sat_pred_biomass)
r2_satellite

field_pred_biomass <- growth_curve(field_fit_result$par, data = field_data_scaled)

r2_field <- calc_r2(field_data_scaled, field_pred_biomass)
r2_field

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Exporting results ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



# Assuming pred = predicted AGB, obs = norm_data$biomass
df <- data.frame(
    Predicted = field_pred_biomass,
    Observed = field_data_scaled$biomass
)

ext <- ggplot(df, aes(x = Predicted, y = Observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 2) +
    labs(
        x = "Predicted Biomass (Mg/ha)",
        y = "Observed Biomass (Mg/ha)"
    ) +
    coord_cartesian(expand = FALSE) +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", size = 28, family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position = "none"
    )

# Save to file
ggsave("./0_results/figures/extended/predicted_vs_observed_field.png",
    plot = ext
)






# Histogram of field ages

ext <- ggplot(field_data, aes(x = age)) +
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

# Save to file
ggsave("./0_results/figures/extended/field_age_histogram.png",
    plot = ext, width = 1800, height = 1400, res = 300
)
