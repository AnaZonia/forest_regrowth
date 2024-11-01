
library(tidyverse)
library(ggplot2)
library(reshape2)

tst <- read.csv("./0_data/non_aggregated_5yr.csv")
tst <- tst %>% filter(biome == 4)

# Extract year columns for variables that vary with year
year_columns <- grep("_(\\d{4})$", names(tst), value = TRUE)
mean_values <- colMeans(tst[, year_columns], na.rm = TRUE)
print(mean_values)
# Convert mean_values into a data frame
mean_values_df <- data.frame(variable = names(mean_values), mean_value = mean_values)

# Extract year and variable names
mean_values_df$year <- as.numeric(sub(".*_(\\d{4})$", "\\1", mean_values_df$variable))
mean_values_df$variable <- sub("_(\\d{4})$", "", mean_values_df$variable)

# Reshape the data to have one row per year and one column per variable
mean_values_wide <- reshape(mean_values_df, idvar = "year", timevar = "variable", direction = "wide")

# Rename columns to remove "mean_value." prefix for clarity
colnames(mean_values_wide) <- sub("mean_value\\.", "", colnames(mean_values_wide))

# Sort by year to ensure rows are in chronological order
mean_values_wide <- mean_values_wide[order(mean_values_wide$year), ]

# Get the list of column names (excluding 'year')
variables <- colnames(mean_values_wide)[colnames(mean_values_wide) != "year"]

dev.off()
par(mfrow = c(3, 3)) # Adjust c(row, col) to match the number of variables

# Loop through each variable and plot
for (var in variables) {
    formula <- as.formula(paste(var, "~ year"))
    model <- lm(formula, data = mean_values_wide)

    # Plot the data with a regression line
    plot(mean_values_wide$year, mean_values_wide[[var]],
        main = paste(var, "vs Year"),
        xlab = "Year", ylab = var, pch = 19
    )
    abline(model, col = "blue", lwd = 2)
}

# check if there is significant annual variation in the data


### Curve fit on medians and on the entire dataset
# plot curves from the non-clustered data and the 10 fits


# Calculate median of `pred` by age group and define asymptote
median_pred_per_age <- data %>%
    group_by(age) %>%
    summarize(median_pred = median(pred, na.rm = TRUE))
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


