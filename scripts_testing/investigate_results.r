################################################################################
#                                                                              #
#                 Comparing and investigating data                             #
#                                                                              #
################################################################################

library(corrplot)
library(tidyverse)

source("fit_1_import_data.r")
source("fit_3_run_model.r")
source("fit_2_functions.r")
find_combination_pars(iterations_optim)
data <- dataframes[[2]][[3]]

dataframes_lm <- lapply(datafiles, import_climatic_data, normalize = TRUE, convert_to_dummy = TRUE)
dataframes_lm <- prepare_dataframes(dataframes_lm, c(1, 4))
data_lm <- dataframes_lm[[2]][[3]]
colnames(data_lm)

lm_cv_output <- cross_valid(data_lm, run_lm, c("cwd", "mean_prec", "mean_si"))

pars_iter <- data_pars[[4]]
pars_iter
numeric_cols <- pars_iter[!grepl("prec|si|agbd|biome|distance|soil|ecoreg|LU|B0|k0|theta", pars_iter)]
numeric_cols
tst <- c(numeric_cols, "nearest_mature")
numeric_cols
lm_cv_output <- cross_valid(data_lm, run_lm, tst)

optim_cv_output <- cross_valid(data, run_optim, pars_iter, conditions)

# Remove non-numeric columns
numeric_df <- data_lm %>% select_if(is.numeric)

# Calculate correlation matrix
cor_matrix <- cor(numeric_df)

# Check for zero variance variables
zero_var_cols <- sapply(numeric_df, function(x) var(x, na.rm = TRUE) == 0)
if (any(zero_var_cols)) {
    print("Zero variance columns:")
    print(names(numeric_df)[zero_var_cols])
    data_clean <- numeric_df[, !zero_var_cols]
}
data_clean <- data_clean[!grepl("prec|si|biome|distance", colnames(data_clean))]

cor_matrix <- cor(data_clean, use = "complete.obs")

# Plot the correlation matrix
corrplot(cor_matrix,
    method = "color", type = "upper", order = "original",
    tl.col = "black", tl.srt = 45, addCoef.col = "black",
    number.cex = 0.7, cl.pos = "n"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Results ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


top_5_per_data_pars <- results_all %>%
    group_by(data_pars) %>%
    top_n(5, rsq) %>%
    arrange(data_pars, desc(rsq)) %>%
    ungroup()


mean_rsq_per_data_name <- results_lm %>%
    # filter(model_type == "optim")%>%
    group_by(data_pars) %>%
    summarise(mean_rsq = mean(rsq, na.rm = TRUE)) %>%
    arrange(desc(mean_rsq))

anova_result <- aov(rsq ~ data_name, data = results_all)

print(summary(anova_result))

# Assuming results_all is your dataframe
mean_rsq_per_data_name <- tst %>%
    filter(biome == 1) %>%
    group_by(data_pars) %>%
    summarise(mean_rsq = mean(rsq, na.rm = TRUE)) %>%
    arrange(desc(mean_rsq))

best_rsq_df <- lm_fit %>%
    # group_by(biome, data_name, data_pars) %>%
    slice_max(rsq, with_ties = FALSE) # %>%
# arrange(biome, data_pars, data_name)
best_rsq_df[, c(1:8)]
best_rsq_df
