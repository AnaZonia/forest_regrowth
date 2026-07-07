# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#       Get average biomass per year with EU TMF data
#
#                 Ana Avila - November 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(tidyverse)

csv_files <- list.files(paste0("./0_data/tmf_biomass"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()


df <- df[!is.na(df$tmf),]

# Ensure df_summary includes ymin/ymax for error bars
df_summary <- df %>%
    group_by(tmf) %>%
    summarise(
        mean_biomass = mean(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    ) %>%
    mutate(
        ymin = mean_biomass - sd_biomass,
        ymax = mean_biomass + sd_biomass
    ) %>%
    ungroup()


ext <- ggplot(df_summary, aes(x = tmf, y = mean_biomass)) +
    geom_line(color = "black", size = 1.2) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
        alpha = 0.2) +
    labs(
        x = "Forest Age estimated by EU TMF data (years)",
        y = expression("Mean Biomass (Mg " * ha^{
            -1
        } * ")")
    ) +
    theme_pubr(base_size = 20) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.8),
        plot.title = element_text(face = "bold")
    )




ext


# Save to file
ggsave("./0_results/figures/extended/tmf_average_biomass_per_age.png",
    plot = ext
)
