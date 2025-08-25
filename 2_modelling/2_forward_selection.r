# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                  Forward Selection for Optim
#
#                    Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------- Identify Optimal Parameter Combination ------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Identify Optimal Parameter Combination
#'
#' Iteratively fits parameter combinations with `run_optim()` and selects the one
#' that minimizes the Akaike Information Criterion (AIC), excluding parameters that
#' do not improve the model.
#'
#' @param basic_pars Character vector. The basic parameters (lag, k0, B0).
#' @param data_pars Character vector. The data parameters (predictors).
#' @param data Data frame. The dataset to be used for fitting the model.
#'
#' @return A list. The best parameter combination identified.
#'
#' @details
#' - Categorical variables are handled by grouping their dummy variables together.



find_combination_pars <- function(basic_pars, data_pars, data) {

    # Initialize parameter vector with data parameters
    all_pars <- c(setNames(
        rep(0, length(data_pars)),
        c(data_pars)
    ))

    all_pars[["k0"]] <- 1

    if ("lag" %in% basic_pars) {
        all_pars[["lag"]] <- 2.5
    }

    if ("theta" %in% basic_pars) {
        all_pars[["theta"]] <- 1
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Helper function to group dummy variables
    for (var in categorical) {
        dummy_indices <- grep(paste0(var, "_"), data_pars)
        if (length(dummy_indices) > 0) {
            data_pars <- c(data_pars[-dummy_indices], var)
        }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initialize the best model with basic parameters
    remaining <- 1:length(data_pars)
    taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

    # best model list
    best <- list(AIC = 0)
    best[["par"]] <- all_pars[names(all_pars) %in% basic_pars]

    base_row <- all_pars
    base_row[names(all_pars)] <- NA
    base_row <- c(RSS = 0, base_row)

    should_continue <- TRUE
    # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
    if (length(data_pars) == 0) {
        return(best$par)
    } else {
        for (i in 1:length(data_pars)) {
            if (!should_continue) break

            iter_df <- foreach(j = remaining[-taken]) %dopar% {
                # check for categorical variables (to be included as a group)
                if (data_pars[j] %in% c(categorical)) {
                    inipar <- c(best$par, all_pars[grep(paste0(data_pars[j], "_"), names(all_pars))])
                } else {
                    # as starting point, take the best values from last time
                    inipar <- c(best$par, all_pars[data_pars[j]])
                }

                model <- run_optim(data, inipar, conditions)
                iter_row <- base_row
                iter_row[names(inipar)] <- model$par
                iter_row["RSS"] <- model$value

                return(iter_row)
            }

            iter_df <- as.data.frame(do.call(rbind, iter_df))
            best_model <- which.min(iter_df$RSS)
            best_model_AIC <- 2 * (i + length(best$par)) + nrow(data) * log(iter_df$RSS[best_model] / nrow(data))

            if (best$AIC == 0 | best_model_AIC < best$AIC) {
                best$AIC <- best_model_AIC
                best$par <- iter_df[best_model, names(all_pars)]
                best$par <- Filter(function(x) !is.na(x), best$par)
                taken <- which(sapply(data_pars, function(x) any(grepl(x, names(best$par)))))
                print(paste0(i, " parameters included: ", toString(data_pars[taken])))
            } else {
                not_taken <- data_pars[!data_pars %in% data_pars[taken]]
                print(paste("No improvement. Exiting loop. Parameters not taken:", toString(not_taken)))
                should_continue <- FALSE
            }
        }

        return(best$par)
    }
}


