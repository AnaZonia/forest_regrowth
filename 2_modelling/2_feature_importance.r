library(ggplot2)
library(tidyverse)


# Source other scripts
source("2_R_scripts/1_modelling.r")



# Add this function to your main.R file or a new visualization.R file
plot_variable_importance <- function(importance_df, title = "Variable Importance", use_pct = TRUE) {
  # Determine which column to use
  y_col <- if(use_pct) "importance_pct" else "importance"
  y_lab <- if(use_pct) "Relative Importance (%)" else "Importance"
  
  # Convert variable to factor with levels in descending importance order
  importance_df$variable <- factor(importance_df$variable, 
                                   levels = importance_df$variable[order(-importance_df[[y_col]])])
  
  # Create plot
  p <- ggplot(importance_df, aes(x = variable, y = !!sym(y_col))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = title,
         y = y_lab) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Highlight specific variables if desired
  highlight_vars <- c("num_fires", "sur_cover")
  if (any(highlight_vars %in% importance_df$variable)) {
    highlight_df <- importance_df[importance_df$variable %in% highlight_vars, ]
    p <- p + geom_bar(data = highlight_df, 
                      aes(x = variable, y = !!sym(y_col)),
                      stat = "identity", 
                      fill = "darkred")
  }
  
  return(p)
}








model <- run_optim(data, pars_init, conditions)

# Calculate variable importance using different methods
coef_importance <- analyze_variable_importance(model, data_pars)
perm_importance <- calculate_permutation_importance(model, data, data_pars)
loocv_importance <- calculate_loocv_importance(data, basic_pars, data_pars, conditions)

# Create plots
p1 <- plot_variable_importance(coef_importance, 
                            "Variable Importance (Coefficient Magnitude)")
p2 <- plot_variable_importance(perm_importance, 
                            "Variable Importance (Permutation)")
p3 <- plot_variable_importance(loocv_importance, 
                            "Variable Importance (Leave-One-Out)")

# Save plots
ggsave("results/importance_coefficients.png", p1, width = 10, height = 6)
ggsave("results/importance_permutation.png", p2, width = 10, height = 6)
ggsave("results/importance_loocv.png", p3, width = 10, height = 6)




# Compare importance across methods
plot_importance_comparison <- function(importance_results) {
# Combine all methods into one dataframe
importance_combined <- data.frame()

for (method_name in names(importance_results)) {
method_df <- importance_results[[method_name]]
method_df$method <- method_name
importance_combined <- rbind(importance_combined, method_df)
}

# Plot
ggplot(importance_combined, aes(x = variable, y = importance_pct, fill = method)) +
geom_bar(stat = "identity", position = "dodge") +
labs(title = "Variable Importance Comparison Across Methods",
        x = "Variables",
        y = "Relative Importance (%)") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_brewer(palette = "Set1")
}

# Create an importance heatmap
plot_importance_heatmap <- function(importance_results) {
# Combine all methods into one dataframe
importance_combined <- data.frame()

for (method_name in names(importance_results)) {
method_df <- importance_results[[method_name]]
importance_combined <- rbind(
    importance_combined,
    data.frame(
    variable = method_df$variable,
    method = method_name,
    importance = method_df$importance_pct
    )
)
}

# Plot
ggplot(importance_combined, aes(x = variable, y = method, fill = importance)) +
geom_tile() +
scale_fill_gradient(low = "white", high = "steelblue") +
labs(title = "Variable Importance Heatmap",
        x = "Variables",
        y = "Method",
        fill = "Importance (%)") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
}