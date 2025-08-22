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
    mutate(age = floor(age + 0.5))
    
# Compute age frequencies
age_freq <- field_data %>%
    count(age, name = "freq")

# Add frequency to df
field_data <- field_data %>%
    left_join(age_freq, by = "age") %>%
    mutate(weight = 1 / freq)

# For each site, sample one measurement, with higher chances for rarer ages
field_data <- field_data %>%
    group_by(site_id) %>%
    slice_sample(n = 1, weight_by = weight) %>%
    ungroup() %>%
    select(-biome, -lat, -lon, -site_id, -plot_id) %>%
        mutate(across(any_of(categorical), as.factor))

field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)
# remove points with biomass > 400
field_data <- field_data %>%
    filter(biomass <= 200)


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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Get theta from field data --------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# - We use only age to find the shape parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


field_data_scaled <- apply_min_max_scaling(field_data, train_stats)

# unite field and satellite data
common_cols <- intersect(colnames(field_data_scaled), colnames(sat_data_scaled))
field_data_scaled <- field_data_scaled[, common_cols]
sat_data_scaled1 <- sat_data_scaled[, common_cols]
sat_data_scaled1$satellite <- 1
field_data_scaled$satellite <- 0
unite_field_satellite <- rbind(field_data_scaled, sat_data_scaled1)

init_params <- c(
    k0 = 0.05, # baseline growth rate constant
    theta = 1, # shape parameter
    lag = 2.5
)

growth_curve <- function(pars, data, lag = 0) {
    fit_data_pars <- setdiff(names(pars), c(non_data_pars, "age"))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    k <- rep(pars[["k0"]], nrow(data))

    data <- data %>%
        mutate(age = if_else(satellite == 1, age + lag, age))

    age <- data[["age"]]

    if (length(fit_data_pars) > 0) {
        k <- (k + rowSums(sapply(fit_data_pars, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE))) * (age)
    } else {
        k <- k * age
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)
    k[which(k > 7)] <- 7

    if ("theta" %in% names(pars)) {
        theta <- pars[["theta"]]
    } else {
        theta <- 0.7787
    }

    return((data[["asymptote"]]) * (1 - exp(-k))^theta)
}


field_fit_result <- run_optim(unite_field_satellite, init_params, conditions)


theta <- field_fit_result[["par"]]["theta"]


unite_field_satellite <- unite_field_satellite %>%
    mutate(age = if_else(satellite == 1, age + field_fit_result[["par"]]["lag"], age))


ages <- 1:40
k <- field_fit_result[["par"]]["k0"] * ages

df_growth <- data.frame(
    Age = ages,
    Growth = 234 * (1 - exp(-k))^theta
)

ggplot(df_growth, aes(x = Age, y = Growth)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_point(
        data = unite_field_satellite,
        aes(x = age, y = biomass),
        color = "darkred", size = 2, alpha = 0.7
    ) +
    labs(
        x = "Age (years)",
        y = "Biomass",
        title = paste("Growth Curve (Theta =", round(theta, 2), ")")
    ) +
    theme_minimal()

theta

# run this 5 times and get mean and sd of theta and lag





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(satellite_data))[["all_mean_climate"]],
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

# CSV with the R2 and the theta value
write.csv(data.frame(
    r2_field = r2_field,
    r2_satellite = r2_satellite,
    theta = theta
), file = "./0_results/0_field_results.csv", row.names = FALSE)





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
