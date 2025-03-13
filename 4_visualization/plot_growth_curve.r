# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Load data
new_model_data <- read.csv("mapbiomas_amaz_all_results_NEW.csv")
# new_model_data <- read.csv("./0_results/lagged_nolag_unified_data.csv")

# Compute mean values and standard deviations per age
mean_biomass_data <- new_model_data %>%
    group_by(age) %>%
    summarise(
        mean_pred = mean(pred_lagged, na.rm = TRUE),
        sd_pred = sd(pred_lagged, na.rm = TRUE),
        mean_pred_nolag = mean(pred_nolag, na.rm = TRUE),
        sd_pred_nolag = sd(pred_nolag, na.rm = TRUE),
        mean_biomass = mean(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    ) %>%
    mutate(age_nolag = age + 20) %>% # Shift the ages for mean_pred_nolag
    pivot_longer(
        cols = c(mean_pred, mean_pred_nolag, mean_biomass),
        names_to = "biomass_type",
        values_to = "value"
    ) %>%
    mutate(
        sd = case_when(
            biomass_type == "mean_pred" ~ sd_pred,
            biomass_type == "mean_pred_nolag" ~ sd_pred_nolag,
            biomass_type == "mean_biomass" ~ sd_biomass
        ),
        plot_age = if_else(biomass_type != "mean_pred", age_nolag, age) # Apply shift only for no lag
    )

# Define colors (consistent with previous plot)
colors <- c(
    "mean_pred" = "green",
    "mean_pred_nolag" = "red",
    "mean_biomass" = "blue"
)

# Plot the curves with updated title and legend position
p <- ggplot(mean_biomass_data, aes(x = plot_age, y = value, color = biomass_type)) +
    geom_line(size = 1.5) + # Main lines
    geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = biomass_type), alpha = 0.2, color = NA) + # SD ribbons
    scale_color_manual(values = colors, labels = c("Predicted", "Predicted (No Lag, Shifted)", "Observed Biomass")) +
    scale_fill_manual(values = colors, guide = "none") + # Match ribbon colors but remove legend
    labs(
        title = "Mean Biomass Predictions and Observations",
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)",
        color = "Legend"
    ) +
    theme_minimal(base_size = 20) + # Larger font size
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), # Centered title
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = c(0.25, 0.75), # Top left legend
        legend.background = element_rect(fill = "white", color = "black", size = 1), # Add background box
        aspect.ratio = 1 / 2 # Makes the plot twice as wide as it is tall
    )


# Assuming `mean_biomass_data` already exists and is being used for the main curves

# Load the second dataset and modify as needed
field_biomass <- read.csv("~/Documents/data/field_biomass_with_biome.csv")
field_biomass$field_biom <- field_biomass$field_biom * 0.5

# Aggregate field biomass data by age
aggregate_biomass <- function(data, age_col, biomass_col, interval = 1) {
    data %>%
        mutate(age_interval = floor({{ age_col }} + 0.5)) %>% # Group into integer intervals
        group_by(age_interval) %>%
        summarise(mean_biomass = mean({{ biomass_col }}, na.rm = TRUE)) %>%
        rename(age = age_interval)
}

field_aggregated <- aggregate_biomass(field_biomass, field_age, field_biom)

# Combine the field data with the mean biomass data (assuming both are by age)
combined_data <- left_join(mean_biomass_data, field_aggregated, by = "age")

# Set colors for the plot
colors <- c("blue", "green", "red") # Customize as needed



# Plot the curves with updated title and legend position
p <- ggplot(mean_biomass_data, aes(x = plot_age, y = value, color = biomass_type)) +
    geom_line(size = 1.5) + # Main lines
    geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = biomass_type), alpha = 0.2, color = NA) + # SD ribbons
    scale_color_manual(values = colors, labels = c("Predicted", "Predicted (No Lag, Shifted)", "Observed Biomass")) +
    scale_fill_manual(values = colors, guide = "none") + # Match ribbon colors but remove legend
    labs(
        title = "Mean Biomass Predictions and Observations",
        x = "Forest Age (years)",
        y = "Biomass (MgC/ha)",
        color = "Legend"
    ) +
    # Add the field biomass data (as points or lines)
    geom_point(data = field_aggregated, aes(x = age, y = mean_biomass), color = "purple", size = 2, alpha = 0.7) + # Field data points
    theme_minimal(base_size = 20) + # Larger font size
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), # Centered title
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = c(0.25, 0.75), # Top left legend
        legend.background = element_rect(fill = "white", color = "black", size = 1), # Add background box
        aspect.ratio = 1 / 2 # Makes the plot twice as wide as it is tall
    )

# Print the plot
print(p)

png("biomass_with_scatter.png", width = 10, height = 5, units = "in", res = 900)

dev.off()

# Compute mean values per age
mean_biomass_data <- new_model_data %>%
    group_by(age) %>%
    summarise(
        mean_pred = mean(pred, na.rm = TRUE),
        mean_pred_nolag = mean(pred_nolag, na.rm = TRUE),
        mean_heinrich_biomass_2020 = mean(heinrich_biomass_2020, na.rm = TRUE),
        mean_biomass = mean(biomass, na.rm = TRUE)
    ) %>%
    pivot_longer(
        cols = c(mean_pred, mean_pred_nolag, mean_heinrich_biomass_2020, mean_biomass),
        names_to = "biomass_type",
        values_to = "value"
    )

# Plot the mean curves
ggplot(mean_biomass_data, aes(x = age, y = value, color = biomass_type)) +
    geom_line(size = 1) +
    labs(
        title = "Mean Biomass Predictions and Observations",
        x = "Age",
        y = "Biomass",
        color = "Legend"
    ) +
    scale_color_manual(values = c(
        "mean_pred" = "blue",
        "mean_pred_nolag" = "green",
        "mean_heinrich_biomass_2020" = "red",
        "mean_biomass" = "purple"
    )) +
    theme_minimal()

# Load your data
unified <- read.csv("~/Documents/data/mapbiomas_heinrich_field_2.csv")
unified <- unified %>% rename(age = mapbiomas_age)
write.csv(unified, "~/Documents/data/mapbiomas_heinrich_2.csv")

names(unified)

field_biomass <- read.csv("~/Documents/data/field_biomass_with_biome.csv")
field_biomass$field_biom <- (field_biomass$field_biom*0.5)
field_biomass <- subset(field_biomass, biome == 1)


heinrich_2017 <- read.csv("~/Documents/data/heinrich_2017.csv")
heinrich_2017$ESA_CCI_2017 <- (heinrich_2017$ESA_CCI_2017)

# Calculate mean biomass and mean heinrich_biomass per mapbiomas_age
mean_biomass_2017 <- heinrich_2017 %>%
    group_by(silva_age) %>%
    summarise(
        mean_2017_biomass = median(ESA_CCI_2017, na.rm = TRUE),
        mean_heinrich_2017_biomass = mean(heinrich_biomass_2017, na.rm = TRUE)
    )  %>%
    pivot_longer(
        cols = c(mean_2017_biomass, mean_heinrich_2017_biomass),
        names_to = "biomass_type",
        values_to = "value"
    )

# Calculate mean biomass and mean heinrich_biomass per mapbiomas_age
mean_biomass <- unified %>%
    group_by(silva_age) %>%
    summarise(
        mean_2020_biomass = median(biomass, na.rm = TRUE),
        mean_heinrich_2020_biomass = mean(heinrich_biomass_2020, na.rm = TRUE)
    )

# mean_biomass$mean_biomass <- mean_biomass$mean_biomass-55

# Combine the datasets for plotting
combined_data <- mean_biomass %>%
    pivot_longer(
        cols = c(mean_biomass, mean_heinrich_biomass),
        names_to = "biomass_type",
        values_to = "value"
    ) %>%
    rename(age = silva_age) %>%
    left_join(mean_field_biomass, by = c("age" = "field_age"))


# Compute mean and standard deviation per age
mean_biomass_data <- new_model_data %>%
    group_by(age) %>%
    summarise(
        mean_pred = mean(pred, na.rm = TRUE),
        sd_pred = sd(pred, na.rm = TRUE),
        mean_pred_nolag = mean(pred_nolag, na.rm = TRUE),
        sd_pred_nolag = sd(pred_nolag, na.rm = TRUE),
        mean_heinrich_biomass_2020 = mean(heinrich_biomass_2020, na.rm = TRUE),
        sd_heinrich_biomass_2020 = sd(heinrich_biomass_2020, na.rm = TRUE),
        mean_biomass = mean(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    ) %>%
    pivot_longer(
        cols = starts_with("mean_"),
        names_to = "biomass_type",
        values_to = "mean_value"
    ) %>%
    mutate(sd_value = case_when(
        biomass_type == "mean_pred" ~ sd_pred,
        biomass_type == "mean_pred_nolag" ~ sd_pred_nolag,
        biomass_type == "mean_heinrich_biomass_2020" ~ sd_heinrich_biomass_2020,
        biomass_type == "mean_biomass" ~ sd_biomass,
        TRUE ~ NA_real_
    ))

# Plot the mean curves with standard deviation as ribbons
ggplot(mean_biomass_data, aes(x = age, y = mean_value, color = biomass_type)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = biomass_type), alpha = 0.2) +
    geom_point(data = field_biomass, aes(x = field_age, y = field_biom), color = "black", alpha = 0.7, size = 2) +  # Scatterplot for field data
    labs(
        title = "Mean Biomass Predictions with Standard Deviation",
        x = "Age",
        y = "Biomass",
        color = "Legend",
        fill = "Legend"
    ) +
    scale_color_manual(values = c(
        "mean_pred" = "blue",
        "mean_pred_nolag" = "green",
        "mean_heinrich_biomass_2020" = "red",
        "mean_biomass" = "purple"
    )) +
    scale_fill_manual(values = c(
        "mean_pred" = "blue",
        "mean_pred_nolag" = "green",
        "mean_heinrich_biomass_2020" = "red",
        "mean_biomass" = "purple"
    )) +
    theme_minimal()



new_model_data <- read.csv("~/Documents/data/mapbiomas_ESA_SD.csv")

# Calculate mean biomass and mean heinrich_biomass per mapbiomas_age
mean_biomass <- new_model_data %>%
    group_by(mapbiomas_age) %>%
    summarise(
        mean_AGB = mean(AGB, na.rm = TRUE),
        mean_SD = mean(SD, na.rm = TRUE)
    )

ggplot(mean_biomass, aes(x = mapbiomas_age, y = mean_AGB)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_AGB - mean_SD, ymax = mean_AGB + mean_SD), alpha = 0.2) +
    labs(
        title = "Mean Biomass Predictions with Standard Deviation",
        x = "Age",
        y = "Biomass"
    ) +
    theme_minimal()

tst = subset(new_model_data, mapbiomas_age == 2)

hist(tst$SD)
mean(tst$SD)

# Filter forests with age greater than 25
filtered_data <- subset(new_model_data, mapbiomas_age > 25)

# Split data into two groups based on SD (50% higher and 50% lower SD)
high_sd_group <- subset(filtered_data, SD > median(SD, na.rm = TRUE))
low_sd_group <- subset(filtered_data, SD <= median(SD, na.rm = TRUE))

# Compute mean biomass for each group
mean_high_sd_biomass <- mean(high_sd_group$AGB, na.rm = TRUE)
mean_low_sd_biomass <- mean(low_sd_group$AGB, na.rm = TRUE)

mean_high_sd_biomass
mean_low_sd_biomass

filtered_data <- subset(new_model_data, SD < 50)

# Calculate mean biomass and mean heinrich_biomass per mapbiomas_age
mean_biomass <- filtered_data %>%
    group_by(mapbiomas_age) %>%
    summarise(
        mean_AGB = mean(AGB, na.rm = TRUE),
        mean_SD = mean(SD, na.rm = TRUE)
    )

ggplot(mean_biomass, aes(x = mapbiomas_age, y = mean_AGB)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_AGB - mean_SD, ymax = mean_AGB + mean_SD), alpha = 0.2) +
    labs(
        title = "Mean Biomass Predictions with Standard Deviation",
        x = "Age",
        y = "Biomass"
    ) +
    theme_minimal()

