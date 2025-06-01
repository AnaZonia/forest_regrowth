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

secondary <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000)

vars <- c("srad", "soil", "temp", "vpd", "aet", "def", "pr", "pdsi")

long_data <- secondary %>%
    pivot_longer(
        cols = matches(paste0("^(", paste(vars, collapse = "|"), ")_\\d+$")),
        names_to = c("variable", "year"),
        names_sep = "_",
        values_to = "value"
    ) %>%
    mutate(year = as.integer(year))

trend_data <- long_data %>%
    group_by(variable, year) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

ggplot(trend_data, aes(x = year, y = mean_value)) +
    geom_line(color = "steelblue") +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Climate variable trends", y = "Mean value", x = "Year") +
    theme_minimal()





secondary_CMIP6_1 <- import_data("grid_10k_amazon_secondary_CMIP6", biome = 1, n_samples = 10000)

secondary_CMIP6 <- secondary_CMIP6 %>%
    select(-matches("(_1|_2)$"))

climatic_pars <- c("pr", "nssh", "musc", "sdsr", "nsat")
climatic_pars <- c("srad", "temp", "def", "vpd", "pr", "pdsi", "aet")


# For each variable, compute row-wise historical and future means
for (v in vars) {
    # Identify columns for historical and future years
    historical_cols <- grep(paste0("^", v, "_(19\\d{2}|20(0\\d|1[0-9]))(_\\d)?$"), names(secondary_CMIP6), value = TRUE)
    future_cols <- grep(paste0("^", v, "_(202\\d|203\\d|204\\d)(_\\d)?$"), names(secondary_CMIP6), value = TRUE)

    # Compute row-wise means and add new columns
    secondary_CMIP6[[paste0("mean_", v, "_historical")]] <- rowMeans(secondary_CMIP6[, historical_cols], na.rm = TRUE)
    secondary_CMIP6[[paste0("mean_", v, "_future")]] <- rowMeans(secondary_CMIP6[, future_cols], na.rm = TRUE)
}

norm_data <- normalize_independently(secondary_CMIP6)
norm_data <- norm_data$train_data


data_pars = colnames(norm_data)[grepl(paste(c(paste0("mean_", climatic_pars, "_historical")), collapse = "|"), colnames(norm_data))]
basic_pars <- basic_pars_options[["lag"]]
init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)






secondary_CMIP6_clean <- secondary_CMIP6 %>%
    select(-matches("(_1|_2)$"))


vars <- c("pr", "nssh", "musc", "sdsr", "nsat")

# Pivot longer and parse variable, year, scenario
long_data <- secondary_CMIP6 %>%
    select(-matches("(_1$|_2$)")) %>%
    bind_cols(
        secondary_CMIP6 %>%
            select(matches(paste0("(", paste(vars, collapse = "|"), ")_\\d{4}(_\\d)?$")))
    ) %>%
    pivot_longer(
        cols = everything(),
        names_to = "var_full",
        values_to = "value"
    ) %>%
    mutate(
        variable = str_extract(var_full, paste(vars, collapse = "|")),
        year = as.integer(str_extract(var_full, "\\d{4}")),
        scenario = case_when(
            str_detect(var_full, "_2$") ~ "Scenario 3",
            str_detect(var_full, "_1$") ~ "Scenario 2",
            TRUE ~ "Scenario 1"
        )
    )

# Average over samples
trend_data <- long_data %>%
    group_by(variable, year, scenario) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Plot with 5 panels and color by scenario
ggplot(trend_data, aes(x = year, y = mean_value, color = scenario)) +
    geom_line(size = 1) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Future climate projections (3 scenarios)", x = "Year", y = "Mean value") +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2")






