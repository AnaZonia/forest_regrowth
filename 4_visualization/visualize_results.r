
library(tidyverse)
library(ggplot2)
# library(reshape2)

data <- read.csv("./0_data/unified_fc_cleaned.csv") %>%
    filter(biome == 1)
data$biomass <- data$ESA_CCI_2020

# Calculate median of `pred` by age group and define asymptote
median_biomass_per_age <- data %>%
    group_by(age) %>%
    summarize(median_biomass = median(biomass, na.rm = TRUE))

ggplot(median_biomass_per_age, aes(x = age, y = median_biomass)) +
    geom_line() +
    labs(
        title = "Median Biomass per Age",
        x = "Age",
        y = "Median Biomass"
    ) +
    theme_minimal()


asymptote <- median(data$nearest_mature_biomass, na.rm = TRUE)

# Extend age range and join with `median_pred_per_age`
ages_full <- data.frame(age = 1:130)
median_pred_per_age$age <- median_pred_per_age$age + 34
median_pred_per_age_extended <- ages_full %>%
    left_join(median_pred_per_age, by = "age")

# Model fitting using `optim` and Chapman-Richards growth model
pars <- c(k = 1, theta = 1)
simple_curve <- function(pars, data, conditions) {
    curve <- asymptote * (1 - exp(-pars[["k"]] * data[["age"]]))^pars[["theta"]]
    result <- sum((curve - data$median_pred)^2)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    }
    if (is.na(result) || result == 0) {
        return(-Inf)
    }
    return(result)
}

model <- optim(pars, simple_curve, data = median_pred_per_age_extended[1:68, ], conditions = conditions)

# Predict values with fitted parameters
median_pred_per_age_extended$predicted_median <- asymptote * (1 - exp(-model$par[["k"]] * median_pred_per_age_extended[["age"]]))^model$par[["theta"]]

# Additional AGBD Calculation and Plotting
median_agbd_per_age <- data %>%
    group_by(age) %>%
    summarize(median_agbd = median(agbd, na.rm = TRUE))
median_agbd_per_age$age <- median_agbd_per_age$age + 34
median_pred_per_age_extended <- median_pred_per_age_extended %>%
    left_join(median_agbd_per_age, by = "age")

# Plotting
ggplot(median_pred_per_age_extended, aes(x = age)) +
    geom_line(aes(y = median_pred), color = "blue", linetype = "dashed", size = 1, na.rm = TRUE) +
    geom_line(aes(y = predicted_median), color = "red", size = 1) +
    geom_line(aes(y = median_agbd), color = "black", linetype = "dashed", size = 1) +
    labs(
        title = "Median `pred` per Age with Extended Chapman-Richards Growth Fit",
        x = "Age",
        y = "Median `pred`"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("Original Median" = "blue", "Fitted Curve" = "red"))



# plot heatmaps of r squared values to show all the comparisons
# across models and configurations
