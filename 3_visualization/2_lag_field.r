# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(patchwork) # For legend extraction

# Set global options
options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 20))

# ---------------------------- Data Loading ----------------------------
field_data <- read.csv("0_data/groa_field/aggregated_field_biomass.csv")
lag_data <- read.csv("0_results/pred_vs_obs_amazon_lag.csv")
intercept_data <- read.csv("0_results/pred_vs_obs_amazon_intercept.csv")

# Calculate the age lag
lag <- lag_data$corrected_age[1] - lag_data$uncorrected_age[1]

# ---------------------------- Data Preparation ----------------------------
prepare_summary_data <- function(data, group_var, pred_var, obs_var = NULL) {
    summary <- data %>%
        group_by({{ group_var }}) %>%
        summarise(
            mean_pred = median({{ pred_var }}, na.rm = TRUE),
            sd_pred = sd({{ pred_var }}, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        rename(age = {{ group_var }})

    if (!rlang::quo_is_null(enquo(obs_var))) {
        obs_summary <- data %>%
            group_by({{ group_var }}) %>%
            summarise(
                mean_obs = median({{ obs_var }}, na.rm = TRUE),
                sd_obs = sd({{ obs_var }}, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            rename(age = {{ group_var }})

        summary <- full_join(summary, obs_summary, by = "age")
    }

    return(summary)
}

# Prepare all datasets
uncorrected_summary <- prepare_summary_data(
  lag_data, corrected_age, pred_uncorrected, obs
)

lag_corrected_summary <- prepare_summary_data(
  lag_data, uncorrected_age, pred_corrected
) %>% filter(age <= (lag + 1))

lag_future <- prepare_summary_data(
    lag_data, age_future, pred_future
)

intercept_summary <- prepare_summary_data(
  intercept_data, age, pred_corrected
)

intercept_future <- prepare_summary_data(
    intercept_data, age_future, pred_future
)

# First, combine all data including future predictions
all_pred_data <- uncorrected_summary %>%
    full_join(lag_corrected_summary, by = "age", suffix = c("_uncorrected", "_lag_corrected")) %>%
    full_join(intercept_summary, by = "age") %>%
    rename(
        mean_pred_intercept = mean_pred,
        sd_pred_intercept = sd_pred
    ) %>%
    full_join(
        lag_future %>% rename(
            mean_pred_lag_future = mean_pred,
            sd_pred_lag_future = sd_pred
        ),
        by = "age"
    ) %>%
    full_join(
        intercept_future %>% rename(
            mean_pred_intercept_future = mean_pred,
            sd_pred_intercept_future = sd_pred
        ),
        by = "age"
    )

# Update plot colors and linetypes to include future predictions
plot_colors <- c(
    "Predicted (lag-corrected)" = "blue",
    "Predicted (uncorrected)" = "blue",
    "Predicted (intercept-only)" = "green4",
    "Remote Sensing" = "red",
    "Field Measurements" = "black",
    "Future (lag-corrected)" = "blue",
    "Future (intercept-only)" = "green4"
)

linetypes <- c(
    "Predicted (lag-corrected)" = "dotted",
    "Predicted (uncorrected)" = "solid",
    "Predicted (intercept-only)" = "dashed",
    "Remote Sensing" = "solid",
    "Field Measurements" = "solid",
    "Future (lag-corrected)" = "solid",
    "Future (intercept-only)" = "dashed"
)

# Modify the plotting function
create_biomass_plot <- function(data, field_data, lag, output_width = 20, output_height = 8) {
    p <- ggplot(data, aes(x = age)) +
        # Remote sensing data
        geom_line(aes(y = mean_obs, color = "Remote Sensing"), size = 1) +
        geom_ribbon(
            aes(ymin = mean_obs - sd_obs, ymax = mean_obs + sd_obs, fill = "Remote Sensing"),
            alpha = 0.2, color = NA
        ) +

        # Predicted data
        geom_line(aes(
            y = mean_pred_lag_corrected, color = "Predicted (lag-corrected)",
            linetype = "Predicted (lag-corrected)"
        ), size = 1) +
        geom_ribbon(
            aes(
                ymin = mean_pred_lag_corrected - sd_pred_lag_corrected,
                ymax = mean_pred_lag_corrected + sd_pred_lag_corrected,
                fill = "Predicted (lag-corrected)"
            ),
            alpha = 0.2, color = NA
        ) +
        geom_line(aes(
            y = mean_pred_uncorrected, color = "Predicted (uncorrected)",
            linetype = "Predicted (uncorrected)"
        ), size = 1) +
        geom_ribbon(
            aes(
                ymin = mean_pred_uncorrected - sd_pred_uncorrected,
                ymax = mean_pred_uncorrected + sd_pred_uncorrected,
                fill = "Predicted (uncorrected)"
            ),
            alpha = 0.2, color = NA
        ) +
        geom_line(aes(
            y = mean_pred_intercept, color = "Predicted (intercept-only)",
            linetype = "Predicted (intercept-only)"
        ), size = 1) +
        geom_ribbon(
            aes(
                ymin = mean_pred_intercept - sd_pred_intercept,
                ymax = mean_pred_intercept + sd_pred_intercept,
                fill = "Predicted (intercept-only)"
            ),
            alpha = 0.2, color = NA
        ) +

        # Field data points
        geom_point(
            data = field_data,
            aes(x = age, y = mean_biomass, color = "Field Measurements"),
            size = 3, alpha = 0.7
        ) +

        # Vertical line for age lag
        geom_vline(
            xintercept = (lag + 1),
            linetype = "dotted", color = "black", size = 1
        ) + 
        geom_text(
            aes(x = (lag + 1), y = 300, label = paste(lag, "year lag")),
            hjust = 0.5, color = "black", size = 5
        ) +
        
        # Add future prediction lines
        geom_line(
            aes(
                y = mean_pred_lag_future,
                color = "Future (lag-corrected)",
                linetype = "Future (lag-corrected)"
            ),
            size = 1, na.rm = TRUE
        ) +
        geom_ribbon(
            aes(
                ymin = mean_pred_lag_future - sd_pred_lag_future,
                ymax = mean_pred_lag_future + sd_pred_lag_future,
                fill = "Future (lag-corrected)"
            ),
            alpha = 0.15, color = NA, na.rm = TRUE
        ) +
        geom_line(
            aes(
                y = mean_pred_intercept_future,
                color = "Future (intercept-only)",
                linetype = "Future (intercept-only)"
            ),
            size = 1, na.rm = TRUE
        ) +
        geom_ribbon(
            aes(
                ymin = mean_pred_intercept_future - sd_pred_intercept_future,
                ymax = mean_pred_intercept_future + sd_pred_intercept_future,
                fill = "Future (intercept-only)"
            ),
            alpha = 0.15, color = NA, na.rm = TRUE
        ) +

        # Scale definitions
        scale_color_manual(values = plot_colors, name = NULL) +
        scale_fill_manual(values = plot_colors, name = NULL) +
        scale_linetype_manual(values = linetypes, name = NULL) +
        scale_y_continuous(limits = c(0, 310), expand = expansion(mult = c(0, 0.05))) +

        # Labels and theme
        labs(
            x = "Forest Age (years)",
            y = "Biomass (Mg/ha)"
        ) +
        theme(
            aspect.ratio = 0.5,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.title = element_text(color = "black", family = "Helvetica"),
            axis.text = element_text(color = "black", size = 14, family = "Helvetica")
        )

    return(p)
}

# Create and save plot
biomass_plot <- create_biomass_plot(all_pred_data, field_data, lag)


# ---------------------------- Save Outputs ----------------------------
# Save main plot
ggsave(
  filename = "0_results/figures/lag_field_biomass.jpeg",
  plot = biomass_plot,
  width = 15,
  height = 8,
  units = "in",
  dpi = 300
)

# Function to extract and save legend
save_legend_separately <- function(plot, filename, width = 10, height = 2) {
  # Extract legend
  legend <- cowplot::get_legend(plot)
  
  # Create a blank plot with just the legend
  legend_plot <- ggpubr::as_ggplot(legend)
  
  # Save the legend
  ggsave(
    filename = filename,
    plot = legend_plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

# Save legend separately
save_legend_separately(
  biomass_plot,
  "0_results/figures/biomass_plot_legend.jpeg"
)
