# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#       Plot the Model with lag and Model without lag models
#
#                 Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(tidyverse)

csv_files <- list.files(paste0("./0_data/tmf_biomass"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()


df <- df[!is.na(df$tmf),]


df_summary <- df %>%
    group_by(tmf) %>%
    summarise(
        mean_biomass = mean(biomass, na.rm = TRUE)
    ) %>%
    ungroup()

head(df_summary)

ggplot(df_summary, aes(x = tmf, y = mean_biomass)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
        title = "Mean Biomass vs Time since last disturbance",
        x = "Time since last disturbance (years)",
        y = "Mean Biomass (Mg/ha)"
    ) +
    theme_minimal()
