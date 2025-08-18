library(tidyverse)

setwd("/home/anazonia/Documents/forest_regrowth")

# List of folders
folders <- c(
    "UTM_stratified_5k",
    # "grid_10k_amazon_secondary",
    "UTM_stratified_15k",
    "UTM_stratified_20k"
    # "UTM_100k",
    # "UTM_100k_2"
)

# Function to read, combine, and filter biomass
read_biomass <- function(folder) {
    files <- list.files(file.path("./0_data", folder), pattern = "\\.csv$", full.names = TRUE)
    df <- map_dfr(files, read_csv) %>% filter(!is.na(biomass))
    df$dataset <- folder
    df
}

# Read all datasets into a named list
df_list <- map(folders, read_biomass) %>% set_names(folders)

# Determine common sample size
n_common <- min(map_int(df_list, nrow))

# Shuffle and truncate all datasets
set.seed(123)
df_list <- map(df_list, ~ sample_n(.x, n_common))

# Running mean
running_mean <- function(x) cumsum(x) / seq_along(x)

# Combine for plotting
df_plot <- imap_dfr(df_list, ~ tibble(
    n = 1:n_common,
    mean_biomass = running_mean(.x$biomass),
    dataset = .y
))

# Compute final means
mean_vals <- imap_dfr(df_list, ~ tibble(
    dataset = .y,
    mean_biomass = mean(.x$biomass)
))

# Plot
ggplot(df_plot, aes(x = n, y = mean_biomass, color = dataset)) +
    geom_line() +
    geom_hline(
        data = mean_vals,
        aes(yintercept = mean_biomass, color = dataset),
        linetype = "dashed", show.legend = FALSE
    ) +
    labs(
        x = "Number of samples", y = "Running mean biomass",
        title = "Convergence of sample mean (CLT effect)"
    ) +
    theme_minimal()





# # Read all CSVs, drop system:index and .geo
# UTM_100k_files <- list.files("./0_data/UTM_100k", pattern = "\\.csv$", full.names = TRUE)

# df_100k <- map_dfr(UTM_100k_files, ~ read_csv(.x) %>%
#     select(-`system.index`, -`.geo`))

# # Check columns
# glimpse(df_100k)
