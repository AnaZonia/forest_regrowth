

library(tidyverse)



csv_files <- list.files(paste0("./0_data/grid_1k_amazon_secondary_allpixels"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()


edge <- df %>% filter(edge == 0)
non_edge <- df %>% filter(edge == 1)


edge_0 <- df %>%
    filter(edge == 0) %>%
    slice_sample(n = 20000)

edge_1 <- df %>%
    filter(edge == 1) %>%
    slice_sample(n = 20000)

sampled_df <- bind_rows(edge_0, edge_1)

p <- ggplot(sampled_df, aes(x = biomass, color = factor(edge))) +
    geom_histogram(
        aes(y = after_stat(density), fill = factor(edge)),
        position = "identity",
        alpha = 0.4,
        bins = 50,
        color = "white",
        linewidth = 0.3
    ) +
    geom_density(size = 1.1, fill = NA) +
    scale_fill_manual(
        values = c("0" = "#4575b4", "1" = "#d73027"),
        breaks = c("0", "1"),
        labels = c("Edges", "Interior"),
        name = ""
    ) +
    scale_color_manual(
        values = c("0" = "#4575b4", "1" = "#d73027"),
        breaks = c("0", "1"),
        labels = c("Edges", "Interior"),
        name = ""
    ) +
    labs(
        x = "Biomass (Mg/ha)",
        y = "Density"
    ) +
    coord_cartesian(expand = FALSE) +
    theme(
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", size = 24, family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position  = c(0.8, 0.8),  # upper-right inside plot
        legend.justification = c("right", "top"),
        legend.title = element_text(size = 18, family = "Helvetica"),
        legend.text = element_text(size = 16, family = "Helvetica")
    )

ggsave(
    filename = "0_results/figures/extended/edge_biomass.jpeg",
    plot = p,
    width = 15,
    height = 8,
    units = "in",
    dpi = 300
)
