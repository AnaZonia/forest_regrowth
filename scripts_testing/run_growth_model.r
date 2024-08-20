library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)

# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")
datafiles_amaz <- paste0("./data/amaz_", intervals, ".csv")
dataframes_amaz <- lapply(datafiles_amaz, import_climatic_data, normalize = TRUE)
datafiles_countrywide <- paste0("./data/countrywide_", intervals, ".csv")
dataframes_countrywide <- lapply(datafiles_countrywide, import_climatic_data, normalize = TRUE)

run_growth_model <- function(data, test_data, pars) {

    conditions <- list(
        'pars["theta"] > 10'
        ,'pars["theta"] < 0'
        ,'pars["B0"] < 0'
        # ,'pars["age"] < 0'
        # ,'pars["age"] > 5'
    )

    growth_curve <- function(pars, data) {
        # k <- c(pars[["age"]] * data[["age"]] + pars[["cwd"]] * data[["cwd"]])
        k <- rowSums(sapply(lu_pars, function(par) pars[[par]] * data[[par]] * data[["age"]]))

        pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]]
    }

# head(growth_curve(pars, data))
    likelihood <- function(pars, data, conditions) {
        result <- sum((growth_curve(pars, data) - data$agbd)^2)

        if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
            return(Inf)
        } else if (is.na(result) || result == 0) {
            return(Inf)
        } else {
            return(result)
        }
    }

    o <- optim(pars, likelihood, data = data, conditions = conditions)

    pred <- growth_curve(o$par, test_data)
    # rsq <- cor(data$agbd, pred)^2 #rsq is NA when pred is all the same value

    Model <- lm(test_data$agbd ~ pred)
    Residuals <- summary(Model)$residuals
    SumResSquared <- sum(Residuals^2)
    TotalSumSquares <- sum((test_data$agbd - mean(test_data$agbd))^2)
    rsq <- 1 - (SumResSquared / TotalSumSquares)
    rsq

    return(list(
        optimized_parameters = o$par,
        rsq = rsq,
        predictions = pred,
        optimization_result = o
    ))
}

run_lm <- function(train_data, test_data, pars) {
    lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

    model <- lm(lm_formula, data = train_data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Check for rank deficiency, and remove variables if necessary
    # (usually ones with too few unique values)
    aliased_vars <- summary(model)$aliased

    if (any(aliased_vars)) {
        problematic_vars <- names(aliased_vars)[aliased_vars]
        print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

        # Remove problematic variables from the formula
        pars <- pars[!pars %in% problematic_vars]
        lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))
        model <- lm(lm_formula, data = train_data)
    }

    Model <- lm(test_data$agbd ~ predict(model, newdata = test_data))
    Residuals <- summary(Model)$residuals
    SumResSquared <- sum(Residuals^2)
    TotalSumSquares <- sum((test_data$agbd - mean(test_data$agbd))^2)
    rsq <- 1 - (SumResSquared / TotalSumSquares)

    return(list(
        pars = t(summary(model)$coefficients[-1, 1, drop = FALSE]), # -1 to remove (Intercept),
        rsq = rsq
    ))
}

pars_fit <- readRDS("./data/amaz_ideal_par_combination.rds")

lu_pars <- names(pars_fit[[81]])[!names(pars_fit[[81]]) %in% c("theta", "B0_exp", "k0")]

# pars <- c(setNames(
#     rep(0, length(lu_pars)),
#     c(lu_pars)
# ))
# pars <- c(B0 = mean(dataframes_amaz[[i]][["agbd"]]), theta = 1, pars) # age = 0, cwd = 0)

data <- dataframes_amaz[[1]]
pars <- pars_fit[[81]]

r_squared_values_opt <- c()
r_squared_values_lm <- c()

indices <- sample(c(1:5), nrow(data), replace = TRUE)

for (i in 1:5) {
    print(i)
    # Define the test and train sets
    test_data <- data[indices == i, ]
    train_data <- data[!indices == i, ]

    modellm <- run_lm(train_data, test_data, lu_pars)
    print(modellm$rsq)
    print(modellm$pars)

    modelopt <- run_growth_model(train_data, test_data, pars)
    print(modelopt$rsq)
    print(modelopt$optimized_parameters)

}


mean_r_squared <- mean(r_squared_values)
sd_r_squared <- sd(r_squared_values)

# for (i in seq_along(dataframes_amaz)) {
#   print("----------------------------------------------------")
#   print(i)

#   # Without asymptote (using mature_biomass)
#   result1 <- run_growth_model(dataframes_amaz[[i]], c(B0 = mean(dataframes_amaz[[i]][["agbd"]]), theta = 1, age = 0, cwd = 0))
#   print(paste("R-squared amaz:", result1$r_squared))

#   # # With asymptote
#   result2 <- run_growth_model(dataframes_countrywide[[i]], c(B0 = 40, theta = 5, age = 0, cwd = 0))
#   print(paste("R-squared countrywide_amaz_subset:", result2$r_squared))
# }
