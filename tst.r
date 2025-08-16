library(tidyverse)
setwd("/home/anazonia/Documents/forest_regrowth")

# Read datasets
UTM_10k_files <- list.files("./0_data/UTM_stratified_10k", pattern = "\\.csv$", full.names = TRUE)
df_10k <- map_dfr(UTM_10k_files, read_csv) %>% filter(!is.na(biomass))

UTM_15k_files <- list.files("./0_data/UTM_stratified_15k", pattern = "\\.csv$", full.names = TRUE)
df_15k <- map_dfr(UTM_15k_files, read_csv) %>% filter(!is.na(biomass))

UTM_20k_files <- list.files("./0_data/UTM_stratified_20k", pattern = "\\.csv$", full.names = TRUE)
df_20k <- map_dfr(UTM_20k_files, read_csv) %>% filter(!is.na(biomass))

UTM_5k_files <- list.files("./0_data/UTM_stratified_5k", pattern = "\\.csv$", full.names = TRUE)
df_5k <- map_dfr(UTM_5k_files, read_csv) %>% filter(!is.na(biomass))

UTM_100k_files <- list.files("./0_data/UTM_100k", pattern = "\\.csv$", full.names = TRUE)
df_100k <- map_dfr(UTM_100k_files, read_csv) %>% filter(!is.na(biomass))


# Determine common sample size
n_common <- min(nrow(df_5k), nrow(df_10k), nrow(df_15k), nrow(df_20k), nrow(df_100k))

# Shuffle and truncate
set.seed(123)
df_5k <- df_5k %>% sample_n(n_common)
df_10k <- df_10k %>% sample_n(n_common)
df_15k <- df_15k %>% sample_n(n_common)
df_20k <- df_20k %>% sample_n(n_common)
df_100k <- df_100k %>% sample_n(n_common)

# Running mean function
running_mean <- function(x) cumsum(x) / seq_along(x)

# Combine running mean data for plotting
df_plot <- bind_rows(
    tibble(n = 1:n_common, mean_biomass = running_mean(df_5k$biomass), dataset = "UTM_5k"),
    tibble(n = 1:n_common, mean_biomass = running_mean(df_10k$biomass), dataset = "UTM_10k"),
    tibble(n = 1:n_common, mean_biomass = running_mean(df_15k$biomass), dataset = "UTM_15k"),
    tibble(n = 1:n_common, mean_biomass = running_mean(df_20k$biomass), dataset = "UTM_20k"),
    tibble(n = 1:n_common, mean_biomass = running_mean(df_100k$biomass), dataset = "UTM_100k")
)

# Compute final means
mean_vals <- tibble(
    dataset = c("UTM_5K", "UTM_10k", "UTM_15k", "UTM_20k", "UTM_100k"),
    mean_biomass = c(mean(df_5k$biomass), mean(df_10k$biomass), mean(df_15k$biomass), mean(df_20k$biomass), mean(df_100k$biomass))
)

# Plot with dashed lines for means
ggplot(df_plot, aes(x = n, y = mean_biomass, color = dataset)) +
    geom_line() +
    geom_hline(
        data = mean_vals, aes(yintercept = mean_biomass, color = dataset),
        linetype = "dashed", show.legend = FALSE
    ) +
    labs(
        x = "Number of samples", y = "Running mean biomass",
        title = "Convergence of sample mean (CLT effect)"
    ) +
    theme_minimal()


library(dplyr)

# List all CSV files in the folder
files <- list.files("./0_data/UTM_100k", pattern = "\\.csv$", full.names = TRUE)

# Loop over each file
for (file in files) {
    # Read the CSV
    df <- read.csv(file)

    # Remove .geo and system:index columns, if present
    df <- df %>% select(-any_of(c(".geo", "system.index")))

    # Overwrite the original CSV (or use a different folder for safety)
    write.csv(df, file, row.names = FALSE)
}

