


run_growth_model <- function(data, initial_pars) {
    conditions <- list(
        'pars["theta"] > 10',
        'pars["theta"] < 0',
        'pars["B0"] < 0'
    )

    use_asymptote <- length(initial_pars) > 3

    if (use_asymptote) {
        conditions <- c(conditions, 'pars["B0"] > pars["A"]')
    }

    growth_curve <- if (use_asymptote) {
        function(pars, data) {
            pars[["B0"]] + (pars[["A"]] - pars[["B0"]]) * (1 - exp(-pars[["age"]] * data$age))^pars[["theta"]]
        }
    } else {
        function(pars, data) {
            pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-pars[["age"]]))^pars[["theta"]]
        }
    }

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

    o <- optim(initial_pars, likelihood, data = data, conditions = conditions)

    pred <- growth_curve(o$par, data)
    rsq <- cor(data$agbd, pred)^2

    return(list(
        optimized_parameters = o$par,
        r_squared = rsq,
        predictions = pred,
        optimization_result = o
    ))
}


# # Usage
# for (i in seq_along(dataframes)) {
#   print("----------------------------------------------------")
#   print(names_dataframes[i])

#   # Without asymptote (using mature_biomass)
#   result1 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, age = 0.1))
#   print(paste("R-squared (fixed asymptote, fit growth rate):", result1$r_squared))

#   # # With asymptote
#   # result2 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, k = 0.1, A = 100))
#   # print(paste("R-squared (fit asymptote, rate fit from age):", result2$r_squared))
# }
