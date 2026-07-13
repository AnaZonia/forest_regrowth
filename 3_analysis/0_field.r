# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#           Field Data Analysis and Model Validation
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_error_propagation.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------- Field Data Cleaning ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data <- read.csv("./0_data/groa_field/field_predictors.csv")
field_data <- subset(field_data, biome == 1)


data <- import_data("grid_10k_amazon_uncertainty_propagation", biome = 1, n_samples = 30000, asymptote = "nearest_mature", categorical = categorical)
norm_out <- normalize_independently(data)$train_stats
norm_out <- norm_out[1:20, ]

field_data <- field_data %>%
    rename(
        biomass = field_biom,
        asymptote = nearest_mature
    ) %>%
    mutate(age = floor(age + 0.5)) %>%
    drop_na()

# get median biomass per age per site

field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

field_data <- field_data %>%
    filter(biomass < 400)

field_data_scaled <- apply_min_max_scaling(field_data, norm_out)


error_prop <- readRDS("./0_results/error_prop.rds")
pars <- error_prop$pars

pars$topography_21 <- NULL

r2_list <- numeric(1000)

for (index in 1:nrow(pars)) {
    # index = 1
    pred <- growth_curve(pars[index,], data = field_data_scaled)

    r2 <- calc_r2(field_data_scaled, pred)
    r2_list[index] <- r2
}

results <- data.frame(
    mean_r2 = mean(r2_list),
    sd_r2 = sd(r2_list)
)



write.csv(results, file = "./0_results/field_r2.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Exporting results ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Assuming pred = predicted AGB, obs = norm_data$biomass
df <- data.frame(
    Predicted = pred,
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Histogram of field ages ----------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


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
