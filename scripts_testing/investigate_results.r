################################################################################
#                                                                              #
#                 Comparing and investigating data                             #
#                                                                              #
################################################################################

library(corrplot)
library(tidyverse)

source("fit_3_run_model.r")
source("fit_2_functions.r")
lm_cv_output <- cross_valid(data_lm, run_lm, c("ecoreg"))

find_combination_pars(iterations_optim)

data <- dataframes_lm[[2]][[3]]
colnames(data_lm)

lm_cv_output <- cross_valid(data_lm, run_lm, c("ecoreg"))

indices <- sample(c(1:5), nrow(data), replace = TRUE)
test <- data[indices == 1, ]
train_data <- data[!indices == 1, ]
for (col in names(train_data)) {
    if (is.factor(train_data[[col]])) {
        levels_train <- levels(train_data[[col]])
        print(levels_train)
        test <- test[test[[col]] %in% levels_train, ]
        levels_test <- levels(test[[col]])
        print(levels_test)
    }
}


source("fit_1_import_data.r")




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
