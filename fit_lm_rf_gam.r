library(mgcv)
library(randomForest)
library(tidyverse)

source("1_import_data.r")

data <- import_climatic_data("./data/15y_LULC.csv", normalize = TRUE)
data_raw <- read.csv("./data/15y_LULC.csv")
head(data)
hist(data_raw$agbd)
mean(data2_raw$agbd, na.rm = TRUE)
data2 <- import_climatic_data("./data/15y_LULC_mat_30m.csv", normalize = TRUE)
data2_raw <- read.csv("./data/15y_LULC_mat_30m.csv")
hist(data2_raw$num_fires_before_regrowth)



lulc_unique_values <- lapply(lulc_columns, function(col) unique(data[[col]]))
names(lulc_unique_values) <- lulc_columns
# Print the results
for (col in lulc_columns) {
    cat(col, ":\n")
    print(lulc_unique_values[[col]])
    cat("\n")
}

pars_categ <- c("indig", "protec", names(data)[str_detect(names(data), "LU")])
pars_smooth <- c(
    "num_fires_before_regrowth", "fallow", "sur_cover" # "mean_si", "mean_prec", "cwd"
    # lulc_sum_columns # "lulc_sum_15"
)

run_rf_gam_lm <- function(data, pars_categ, pars_smooth, model) {
    set.seed(123)
    train_indices <- sample(1:nrow(data), size = floor(0.7 * nrow(data)))
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]

    lulc_sum_columns <- names(data)[str_detect(names(data), "lulc")]
    lulc_sum_columns <- lulc_sum_columns[!grepl("lulc_sum_(48|47|40)$", lulc_sum_columns)]

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

val <- run_rf_gam_lm(data, pars_categ, pars_smooth, "lm")

# Plot predictions vs observed
plot(pred_gam, data$agbd)
abline(0, 1, col = "red")



data <- read.csv("./data/mat_for_biomass_and_distance.csv")
data <- data[, 2:3]

# why in here it is showing only a few points so far away? makes no sense.
data <- subset(data, distance < 5000)
plot(data$distance, data$mature_biomass)
