####################################################################
############### Predict mature forest biomass from environmental predictors ##################
####################################################################
# Ana Avila - July 2024
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

library(ggplot2)
library(terra)
library(tidyverse) # for stringr operations
library(mlr) # for createDummyFeatures

#------------------ SWITCHES ------------------#

fit_with_all_parameters <- TRUE

#------------------ BODY ------------------#

data <- read.csv("data/mature_biomass_climate_categ.csv")
# the data comes with yearly columns of seasonality and precipitation. Since for the
# mature forest prediction we only want yearly values, we will calculate the mean of each climatic predictor.
patterns <- c("si_", "prec_")
means <- sapply(patterns, function(pat) rowMeans(data[, grep(pat, names(data))], na.rm = TRUE))
colnames(means) <- c("mean_si", "mean_prec")
data <- cbind(data, means)

# remove unnecessary columns
data <- data[, -grep("prec_|si_|biome|geo|system.index", names(data))]

# turn categorical variables into dummy variables
categorical <- c("ecoreg", "soil")
data[categorical] <- lapply(data[categorical], as.factor)
data <- createDummyFeatures(data, cols = categorical)
data <- data %>%
    rename(agbd = b1, cwd = b1_1)

# Normalize numeric columns (0-1)
# Store max and min values to transform back the predctions
min_agbd <- min(data$agbd, na.rm = TRUE)
max_agbd <- max(data$agbd, na.rm = TRUE)
# transform all numeric columns
numeric_cols <- c("mean_si", "mean_prec", "cwd", "agbd")
data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

if (fit_with_all_parameters) {
    pars_chosen <- names(data)[!names(data) %in% "agbd"]
} else {
    pars_chosen <- c("cwd", "protec", "mean_si", "mean_prec", "indig")
}

pars <- setNames(rep(0.1, length(pars_chosen)), pars_chosen)
pars_categ <- names(data)[!names(data) %in% numeric_cols]

mat_biomass_function <- function(pars, data, pars_chosen) {
    pred_agbd <- data[[1]] * 0
    for (i in seq_along(pars_chosen)) {
        pred_agbd <- pred_agbd + pars[[pars_chosen[i]]] * data[[pars_chosen[i]]]
    }
    return(pred_agbd)
}

likelihood <- function(pars, data, pars_chosen) {
    return(sum((mat_biomass_function(pars, data, pars_chosen) - data$agbd)^2))
}

o <- optim(pars,
    likelihood,
    data = data,
    pars_chosen = pars_chosen)
o

pred <- mat_biomass_function(o$par, data, pars_chosen)
plot(pred, data$agbd)
abline(0, 1)

# Calculate the mean of observed data
mean_agbd <- mean(data$agbd, na.rm = TRUE)

# Calculate Total Sum of Squares (TSS)
TSS <- sum((data$agbd - mean_agbd)^2, na.rm = TRUE)

# Calculate Residual Sum of Squares (RSS)
RSS <- sum((data$agbd - pred)^2, na.rm = TRUE)

# Calculate R-squared
R_squared <- 1 - (RSS / TSS)

# Print R-squared value
print(R_squared)


# Load the mgcv package
library(mgcv)

# Fit a GAM model
# Assuming linear relationships for simplicity; use s(variable) for non-linear terms
gam_model <- gam(agbd ~ s(cwd) + s(mean_si) + s(mean_prec), data = data)

# Construct the formula dynamically
formula <- as.formula(paste("agbd ~ s(cwd) + s(mean_si) + s(mean_prec) + indig + protec +", paste(pars_categ, collapse = " + ")))

# Fit the GAM model with the dynamically created formula
gam_model <- gam(formula, data = data)

# Predict using the GAM model
pred_gam <- predict(gam_model, newdata = data)



summary(gam_model)
predictors <- all.vars(formula(gam_model))
predictors
# Plot predictions vs observed
plot(pred_gam, data$agbd)
abline(0, 1)

# Calculate R-squared for the GAM model
mean_agbd <- mean(data$agbd, na.rm = TRUE)
TSS <- sum((data$agbd - mean_agbd)^2, na.rm = TRUE)
RSS <- sum((data$agbd - pred_gam)^2, na.rm = TRUE)
R_squared_gam <- 1 - (RSS / TSS)
print(R_squared_gam)


# Install and load the randomForest package
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
library(randomForest)

# Assuming 'data' is your dataset and 'agbd' is the response variable
# Split data into training and testing sets
set.seed(123) # For reproducibility
train_indices <- sample(1:nrow(data), size = floor(0.7 * nrow(data)))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Fit a Random Forest model on the training data
rf_model <- randomForest(agbd ~ ., data = train_data, ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE, keep.forest = TRUE, oob.score = TRUE)

# Print model summary including OOB error
print(rf_model)
print(rf_model$oob.error)

# Predict using the Random Forest model on the test data
pred_rf_test <- predict(rf_model, newdata = test_data)

# Calculate R-squared for the Random Forest model on the test data
mean_agbd_test <- mean(test_data$agbd, na.rm = TRUE)
TSS_test <- sum((test_data$agbd - mean_agbd_test)^2, na.rm = TRUE)
RSS_test <- sum((test_data$agbd - pred_rf_test)^2, na.rm = TRUE)
R_squared_rf_test <- 1 - (RSS_test / TSS_test)
print(R_squared_rf_test)

# Plot predictions vs observed for the test data
plot(pred_rf_test, test_data$agbd)
abline(0, 1)