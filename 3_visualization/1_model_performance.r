# ------------------------------------------------- #
# Figure - Model Performance Section
# Plot the importance of each category of variables
# ------------------------------------------------- #



library(ggplot2)
library(tidyverse)


importance_full_amazon <- read.csv("./0_results/importance_full_amazon.csv")
importance_nearest_mature <- read.csv("./0_results/importance_nearest_mature.csv")
importance_nearest_mature

importance_full_amazon$group <- "Amazon"
importance_nearest_mature$group <- "NearestMature"

all_data <- bind_rows(importance_full_amazon, importance_nearest_mature) %>%
    mutate(importance_scaled = importance_pct * r2 / 100) %>%
    filter(importance_scaled > 0)

# ---- Plot ----

ggplot(all_data, aes(x = group, y = importance_scaled, fill = variable)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = expression(R^2 ~ contribution), fill = "Variable") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
        legend.position = "right",
        text = element_text(size = 12)
    )

library(viridis) # install.packages("viridis") if needed

ggplot(all_data, aes(x = group, y = importance_scaled, fill = variable)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = expression(R^2 ~ contribution), fill = "Variable") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_viridis(discrete = TRUE, option = "D") + # "D" is a good option
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "right",
        text = element_text(size = 12),
        panel.grid.major.x = element_blank(), # cleaner vertical grid
        panel.grid.minor = element_blank()
    )
