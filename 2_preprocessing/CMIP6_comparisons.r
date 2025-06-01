# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)


# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection.R")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


secondary_CMIP6 <- import_data("grid_10k_amazon_secondary_CMIP6", biome = 1, n_samples = 10000, asymptote = "quarter_biomass")

terraclim_pars <- c("soil", "temp", "vpd", "aet", "pdsi", "def", "srad")
cmip6_pars <- c("musc", "nsat", "pr", "sdsr", "nssh")

climatic_pars <- c(terraclim_pars, cmip6_pars)

pattern <- paste0("^mean_(", paste(climatic_pars, collapse = "|"), ")$")

secondary_CMIP6_selected <- secondary_CMIP6 %>%
    select(c("age", "asymptote", "biomass", matches(pattern)))

norm_data <- normalize_independently(secondary_CMIP6_selected)
norm_data <- norm_data$train_data
# how many rows have any NA values

# data_pars <- grep(pattern, colnames(secondary_CMIP6_selected), value = TRUE)
# data_pars <- paste0("mean_", c(terraclim_pars, cmip6_pars[-c(4, 5)]))
data_pars <- paste0("mean_", c(terraclim_pars[-c(6, 7)], cmip6_pars))
# data_pars <- paste0("mean_", cmip6_pars[-c(4, 5)])
# data_pars <- paste0("mean_", terraclim_pars)
data_pars
basic_pars <- basic_pars_options[["lag"]]
init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)





climatic_pars <- c("pr", "nssh", "musc", "sdsr", "nsat")

vars <- c(
    "precipitation",
    "near_surface_specific_humidity",
    "moisture_in_upper_portion_of_soil_column",
    "surface_downwelling_shortwave_radiation",
    "near_surface_air_temperature"
)

name_map <- setNames(vars, climatic_pars)

long_data <- secondary_CMIP6_1 %>%
    select(matches(paste0("(", paste(climatic_pars, collapse = "|"), ")_\\d{4}(_\\d)?$"))) %>%
    pivot_longer(
        cols = everything(),
        names_to = "var_full",
        values_to = "value"
    ) %>%
    mutate(
        variable = str_extract(var_full, paste(climatic_pars, collapse = "|")),
        year = as.integer(str_extract(var_full, "\\d{4}")),
        scenario = case_when(
            str_detect(var_full, "_2$") ~ "ssp 585",
            str_detect(var_full, "_1$") ~ "ssp 245",
            TRUE ~ "ssp 126"
        ),
        variable_full = name_map[variable]
    )

trend_data <- long_data %>%
    group_by(variable_full, year, scenario) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

for (var_name in unique(trend_data$variable_full)) {
    p <- trend_data %>%
        filter(variable_full == var_name) %>%
        ggplot(aes(x = year, y = mean_value, color = scenario)) +
        geom_line(size = 1) +
        labs(
            title = paste("Future climate projections for", var_name),
            x = "Year",
            y = "Mean value",
            color = "Scenario"
        ) +
        theme_minimal(base_size = 14) +
        theme(aspect.ratio = 1 / 3) +
        scale_color_brewer(palette = "Dark2")

    ggsave(
        filename = paste0("0_results/figures/", var_name, ".jpeg"),
        plot = p,
        width = 9, height = 3, units = "in", dpi = 300
    )
}
