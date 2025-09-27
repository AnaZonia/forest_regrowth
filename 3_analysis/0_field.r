# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#           Field Data Analysis and Model Validation
#
#                 Ana Avila - September 2025
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

field_data <- field_data %>%
    filter(biomass < 400)


apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}


data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 150000)
indices <- sample(c(1:5), nrow(data), replace = TRUE)

r2_list <- numeric(5)
for (index in 1:5) {
    train_data <- data[indices == index, ]
    norm_data <- normalize_independently(train_data)
    train_stats <- norm_data$train_stats
    norm_data <- norm_data$train_data

    pars_init <- find_combination_pars(
        basic_pars = basic_pars_options[["lag"]],
        data_pars = data_pars_options(colnames(norm_data))[["all"]],
        norm_data
    )

    model <- run_optim(norm_data, pars_init[[1]], conditions)

    field_data_scaled <- apply_min_max_scaling(field_data, train_stats)

    # add column in field_data_scaled for columns present in model$par but not in field_data_scaled, fill with 0
    for (col in (names(model$par)[!names(model$par) %in% c("k0", "lag")])) {
        if (!(col %in% colnames(field_data_scaled))) {
            field_data_scaled[[col]] <- 0
        }
    }

    pred <- growth_curve(model$par, data = field_data_scaled)

    r2 <- calc_r2(field_data_scaled, pred)
    r2_list[index] <- r2
}

results <- data.frame(
    mean_r2 = mean(r2_list),
    sd_r2 = sd(r2_list)
)
write.csv(results, file = "./0_results/0_field_r2.csv", row.names = FALSE)

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
