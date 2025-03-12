library(tidyverse)
library(ggplot2)

# Read and rename columns in each data frame
# data <- read.csv("aggregated_amaz_all_results.csv")
# data <- read.csv("./0_results/results/non_aggregated_all_pred.csv")
data <- read.csv("./0_results/lagged_nolag_unified_data.csv")

# Prepare plot data
plot_data <- data %>%
    group_by(age) %>%
    summarize(across(
        c(biomass, starts_with("pred")),
        list(
            median = ~ median(., na.rm = TRUE),
            sd = ~ sd(., na.rm = TRUE)
        )
     ))# %>%
    # mutate(
    #     age_nolag = age + 30,
    #     age_pred_all = age + 1
    # )

# Bartlow color palette from HCL
bartlow_colors <- c(
    # "Estimated biomass out of observed range" <- hcl(h = 130, c = 70, l = 50), # Green
    "Observed Biomass" = hcl(h = 255, c = 70, l = 50), # Blue
    "Predicted Biomass" = hcl(h = 10, c = 70, l = 50) # Red

)


# Create the plot
p <- ggplot(plot_data) +
    # Biomass
    geom_ribbon(
        aes(
            x = age,
            ymin = biomass_median - biomass_sd,
            ymax = biomass_median + biomass_sd,
            fill = "Observed Biomass"
        ),
        alpha = 0.2
    ) +
    geom_line(aes(x = age, y = biomass_median, color = "Observed Biomass"),
        size = 1.2
    ) +

    # Pred_all_nolag
    geom_ribbon(
        aes(
            x = age,
            ymin = predictions_median - predictions_sd,
            ymax = predictions_median + predictions_sd,
            fill = "Predicted Biomass"
        ),
        alpha = 0.2
    ) +
    geom_line(aes(x = age, y = predictions_median, color = "Predicted Biomass"),
        size = 1.2, linetype = "dashed"
    ) +

    # Customize the plot
    labs(
        title = "Biomass and Prediction Trends by Adjusted Age",
        subtitle = "Median values with shaded areas for ±1 SD",
        x = "Adjusted Age (years)",
        y = "Biomass / Prediction",
        color = "Legend",
        fill = "Legend"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA, size = 1),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.05, 0.95), # Top-left corner
        legend.justification = c(0, 1), # Align legend box to top-left
        legend.background = element_rect(fill = "white", color = "grey80", size = 0.5),
        legend.box.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size = 14), # Increase legend text size
        legend.title = element_text(size = 16, face = "bold"), # Increase legend title size
        legend.key.size = unit(1.5, "lines") # Increase the size of the legend keys
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values = bartlow_colors) +
    scale_fill_manual(values = bartlow_colors) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed", "dotted"))))

print(p)

# # Create the plot
# p <- ggplot(plot_data) +
#     # Biomass
#     geom_ribbon(
#         aes(
#             x = age_nolag,
#             ymin = biomass_median - biomass_sd,
#             ymax = biomass_median + biomass_sd,
#             fill = "Observed Biomass"
#         ),
#         alpha = 0.2
#     ) +
#     geom_line(aes(x = age_nolag, y = biomass_median, color = "Observed Biomass"),
#         size = 1.2
#     ) +

#     # Pred_all_nolag
#     geom_ribbon(
#         aes(
#             x = age_nolag,
#             ymin = pred_all_nolag_median - pred_all_nolag_sd,
#             ymax = pred_all_nolag_median + pred_all_nolag_sd,
#             fill = "Predicted Biomass"
#         ),
#         alpha = 0.2
#     ) +
#     geom_line(aes(x = age_nolag, y = pred_all_nolag_median, color = "Predicted Biomass"),
#         size = 1.2, linetype = "dashed"
#     ) +

#     # Pred_all
#     geom_ribbon(
#         aes(
#             x = age_pred_all,
#             ymin = pred_all_median - pred_all_sd,
#             ymax = pred_all_median + pred_all_sd,
#             fill = "Estimated biomass out of observed range"
#         ),
#         alpha = 0.2
#     ) +
#     geom_line(aes(x = age_pred_all, y = pred_all_median, color = "Estimated biomass out of observed range"),
#         size = 1.2, linetype = "dotted"
#     ) +

#     # Customize the plot
#     labs(
#         title = "Biomass and Prediction Trends by Adjusted Age",
#         subtitle = "Median values with shaded areas for ±1 SD",
#         x = "Adjusted Age (years)",
#         y = "Biomass / Prediction",
#         color = "Legend",
#         fill = "Legend"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#         plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         plot.subtitle = element_text(size = 16, hjust = 0.5),
#         axis.title = element_text(face = "bold"),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "grey80", fill = NA, size = 1),
#         plot.background = element_rect(fill = "white", color = NA),
#         legend.position = c(0.05, 0.95), # Top-left corner
#         legend.justification = c(0, 1), # Align legend box to top-left
#         legend.background = element_rect(fill = "white", color = "grey80", size = 0.5),
#         legend.box.margin = margin(5, 5, 5, 5),
#         legend.text = element_text(size = 14), # Increase legend text size
#         legend.title = element_text(size = 16, face = "bold"), # Increase legend title size
#         legend.key.size = unit(1.5, "lines") # Increase the size of the legend keys
#     ) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     scale_color_manual(values = bartlow_colors) +
#     scale_fill_manual(values = bartlow_colors) +
#     guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed", "dotted"))))

png("biomass_prediction_trends.png", width = 12, height = 8, units = "in", res = 900)

print(p)

dev.off()
