# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#              Cross-Validation and R-Squared
#
#                   Ana Avila - August 2025
#
#     Evaluates the model performance using 5-fold cross-validation.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------ Calculate R-squared -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Computes the coefficient of determination (R²).
#'
#' @param data A dataframe containing observed response values (`biomass` column).
#' @param pred A numeric vector of predicted values from the model.
#'
#' @return A numeric value: the R-squared indicating goodness-of-fit.

calc_r2 <- function(data, pred) {
    obs_pred <- lm(data$biomass ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$biomass - mean(data$biomass))^2)
    r2 <- 1 - (sum_res_squared / total_sum_squares)

    return(r2)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- 5-Fold Cross-Validation ------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Each fold is used once as a test set while training on
#' the remaining folds. Computes R² values for each fold.
#'
#' @param data A dataframe containing the full dataset to be split.
#' @param basic_pars List of basic parameters to pass to the model.
#' @param data_pars Vector of predictor names to include in the model.
#' @param conditions Additional conditions to pass to optim.
#'
#' @return A numeric vector of length 5 containing the R² for each fold.
#'
#' @details
#' - Each fold is randomly assigned using equal probability.
#' - Training and test sets are normalized independently, with the test set
#'   scaled according to the training set's min/max values.
#' - The model is trained using `run_optim` and evaluated using `calc_r2`.




cross_validate <- function(data, basic_pars, data_pars, conditions) {

    indices <- sample(c(1:5), nrow(data), replace = TRUE)
    data$pred_cv <- NA
    data$pred_final <- NA
    r2_list <- numeric(5)

    for (index in 1:5) {
        # Define the test and train sets
        test_data <- data[indices == index, -grep("pred", names(data))]
        train_data <- data[indices != index, -grep("pred", names(data))]

        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        # Function to perform direct optimization
        pars_init <- find_combination_pars(basic_pars, data_pars, train_data)

        # Run the model function on the training set and evaluate on the test set
        model <- run_optim(train_data, pars_init, conditions)

        pred_cv <- growth_curve(model$par, test_data, lag = model$par["lag"])

        # save the predicted values of each iteration of the cross validation.
        data$pred_cv[indices == index] <- pred_cv
        r2 <- calc_r2(data[indices == index, ], pred_cv)
        r2_list[index] <- r2
        print(r2)
    }

    return(r2_list)
}
