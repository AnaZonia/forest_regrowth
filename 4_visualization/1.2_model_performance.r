# ------------------------------------------------- #
# Figure - Model Performance Section
# Plot the importance of each category of variables
# ------------------------------------------------- #

library(ggplot2)
library(tidyverse)
library(RColorBrewer)



importance_data <- data.frame(
    category = c("age_1", "age_4", "age_n"),
    percentage = c(2, 10, 18)
)


p <- ggplot(importance_data, aes(x = category, y = percentage)) +
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

