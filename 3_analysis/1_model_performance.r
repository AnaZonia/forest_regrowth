# ------------------------------------------------- #
# Figure - Model Performance Section
# Plot the importance of each category of variables
# ------------------------------------------------- #

library(ggplot2)
library(tidyverse)
library(RColorBrewer)

R2_asymptote <- read.csv("./0_results/R2_asymptotes.csv")
# Filter for the variables of interest
R2_asymptote <- R2_asymptote %>%
    filter(asymptote %in% c("nearest_mature", "quarter_biomass", "full_amazon"))

p <- ggplot(
    R2_asymptote %>% mutate(asymptote = reorder(asymptote, mean_r2)),
    aes(x = asymptote, y = mean_r2)
) +
    geom_bar(stat = "identity", fill = "#003f5c") +
    coord_flip() +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )

p
