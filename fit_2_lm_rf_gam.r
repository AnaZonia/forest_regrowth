library(mgcv)
library(randomForest)
library(tidyverse)

source("fit_1_import_data.r")

#------------------- Main Functions -------------------#

run_rf_gam_lm <- function(data, pars_categ, pars_smooth, model) {
    set.seed(123)
    train_indices <- sample(1:nrow(data), size = floor(0.7 * nrow(data)))
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]

    # Fit a GAM model
    if (model == "gam") {
        formula <- as.formula(paste(
            "agbd ~",
            paste(sapply(pars_smooth, function(x) paste0("s(", x, ")")), collapse = " + "),
            "+",
            paste(pars_categ, collapse = " + ")
        ))
        model <- gam(formula, data = train_data)
    } else {
        rf_lm_pars <- c(pars_smooth, pars_categ)
        rf_lm_formula <- as.formula(paste("agbd ~", paste(rf_lm_pars, collapse = " + ")))
        if (model == "lm") {
            model <- lm(rf_lm_formula, data = train_data)
        } else {
            model <- randomForest(rf_lm_formula,
                data = train_data,
                ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE,
                keep.forest = TRUE, oob.score = TRUE
            )
        }
    }

    pred <- predict(model, newdata = test_data)
    rsq <- cor(test_data$agbd, pred)^2
    print(rsq)

    return(results <- list(
        model = model,
        pred = pred,
        test_agbd = test_data$agbd,
        rsq = rsq
    ))
}



intervals <- list(
    "5y",
    "10y",
    "15y",
    "all"
)

datafiles <- lapply(intervals, function(file) {
    paste0("./data/", file, "_LULC_dist_amaz_500.csv")
})

dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# prop_amaz <- lapply(dataframes, function(df) {
#     nrow(df %>% filter(biome == 1)) / nrow(df)
# })
# prop_amaz
# 67% of them are in the Amazon anyways.

pars_categ <- c("indig", "protec", names(data)[str_detect(names(data), "LU")])
# pars_categ <- c()
pars_smooth <- c(
    "num_fires_before_regrowth", "num_fires_after_regrowth",
    "ts_fire_before_regrowth", "ts_fire_after_regrowth", "fallow" # , "mean_si", "mean_prec", "cwd"
    # ,"lulc_sum_9", "lulc_sum_15", "lulc_sum_41", "lulc_sum_21"
)


lm_df <- data.frame()
for (i in seq_along(dataframes)) {
    data <- dataframes[[i]]

    val <- run_rf_gam_lm(data, pars_categ, pars_smooth, "lm")

    new_row <- summary(val$model)$coefficients[-1, 1] # -1 to remove (Intercept)
    new_row <- as.data.frame(t(new_row))

    new_row$model_name <- paste0(intervals[[i]])
    new_row$model_type <- "lm"
    new_row$rsq <- val$rsq

    # Reorder columns to have model_name first and rsq second
    new_row <- new_row %>% select(model_name, model_type, rsq, everything())

    lm_df <- bind_rows(lm_df, new_row)
}
head(lm_df)

# Plot predictions vs observed
plot(pred_gam, data$agbd)
abline(0, 1, col = "red")

lulc_sum_columns <- names(data)[str_detect(names(data), "lulc")]
# lulc_sum_columns <- lulc_sum_columns[!grepl("lulc_sum_(48|47|40)$", lulc_sum_columns)]
lulc_sum_columns

lulc_unique_values <- lapply(lulc_sum_columns, function(col) unique(data[[col]]))
names(lulc_unique_values) <- lulc_sum_columns
# Print the results
for (col in lulc_sum_columns) {
    cat(col, ":\n")
    print(lulc_unique_values[[col]])
    cat("\n")
}



