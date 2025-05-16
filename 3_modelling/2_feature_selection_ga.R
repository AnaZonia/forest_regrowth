# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------- Identify Optimal Parameter Combination -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   This function identifies the optimal combination of parameters for a given dataset
#   by iteratively fitting parameter combinations with run_optim and selecting the one that minimizes
#   the Akaike Information Criterion (AIC).
#
# Arguments:
#   iterations       : A dataframe where each row contains the information needed to perform
#                      one iteration of the parameter selection process, including the land use history
#                      interval, data parameter set, basic parameter set, and biome.
#
# Returns:
#   ideal_par_combination : A list where each element contains the best parameter combination
#                           identified for each iteration in the `iterations` dataframe.
#
# Notes:
#   - Categorical variables are handled by grouping their dummy variables together during optimization.
#   - The function writes the results of each iteration to an RDS file for future use.
# External Functions:
#   run_optim()


find_combination_pars <- function(basic_pars, data_pars, data) {
    # basic_pars = c(basic_pars_options[["lag"]], "theta")
    # data <- norm_field
    # data_pars = c("num_fires", "dist")

    # Initialize parameter vector with data parameters
    all_pars <- c(setNames(
        rep(0, length(data_pars)),
        c(data_pars)
    ))

    all_pars[["k0"]] <- 1

    if ("lag" %in% basic_pars) {
        all_pars[["lag"]] <- 2.5
    } else {
        all_pars[["B0"]] <- mean(data[["biomass"]])
    }

    if ("theta" %in% basic_pars) {
        all_pars[["theta"]] <- 1
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Helper function to group dummy variables
    group_dummies <- function(vars, data_pars) {
        for (var in vars) {
            dummy_indices <- grep(paste0(var, "_"), data_pars)
            if (length(dummy_indices) > 0) {
                data_pars <- c(data_pars[-dummy_indices], var)
            }
        }
        return(data_pars)
    }

    # Process categorical variables
    data_pars <- group_dummies(categorical, data_pars)
    # Process climatic variables
    data_pars <- group_dummies(climatic_pars, data_pars)


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

                # check for categorical variables or yearly climatic variables (to be included as a group)
                if (data_pars[j] %in% c(categorical, climatic_pars)) {
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
            print(best_model_AIC)

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





# OPTIM GA ----------------------------------------------------------------------

# ini_par <- data.frame(B0 = 80, k0 = 1, theta = 2)

# for (j in 2:ncore){
#     ini_par[j, ] <- c(
#         ini_par[1, "B0"] * (1.5 * runif(1) + .5),
#         ini_par[1, "k0"] * (1.5 * runif(1) + .5),
#         ini_par[1, "theta"] * (1.5 * runif(1) + 0.5)
#     )
# }

# conditions <- list('pars[["theta"]] > 10', 'pars[["theta"]] < 0', 'pars[["k0"]] < 0')

# optim_ga(ini_par, norm_data)

# --------------------------------------------------------------------------------



# Genetic Algorithm integrated with optim to reduce the local minima problem by Brian Leung - March 2025.

# Draws, for each individual, a label 1, 2, or 3 with probabilities given by p_change (e.g. 50% mutation, 25% crossover, 25% co‑dominance).

# row 1 is the best one, and the rest are the offspring

change <- function(par, p_change) {
    tp <- sample(1:3, ncore, replace = TRUE, prob = p_change) 
    for (i in 2:ncore) {
        if (tp[i] == 1) { # mutation
            par[i, ] <- par[i, ] * (.8 + .4 * runif(ncol(par)))
        } else if (tp[i] == 2) { # cross over
            s <- sample(1:ncore, 1) # find other parent
            t <- sample(1:ncol(par), 1) # choose 1 trait to swap
            par[i, t] <- par[s, t]
            t <- sample(1:ncol(par), 1) # also mutate one other trait, so that make sure don't have exact duplicates
            par[i, t] <- par[i, t] * (.8 + .4 * runif(1))
        } else { # co-dominance
            s <- sample(1:ncore, 1) # find other parent
            par[i, ] <- (par[i, ] + par[s, ]) / 2
            t <- sample(1:ncol(par), 1) # also mutate one other trait, so that make sure don't have exact duplicates
            par[i, t] <- par[i, t] * (.8 + .4 * runif(1))
        }
    }
    return(par)
}

optim_ga <- function(par, norm_data = norm_data, control = list(), ngen = 50, maxit = 50, p_change = c(.5, .25, .25)) # p_change - mutation, cross over, co-dominance - issue with none, is that all might be identical 
{
#	operations of ga - create population, select who survives based on "fitness". Has mutation and cross-over (which in this case includes independent assortment). Can also mix the values of both (e.g, take the midpoint).
# keep the best performer. Choose the other ones based on their RSSs. 
	
    # can be passed into optim to control the number of iterations per optim run
    control$maxit = maxit
	
    options(cores = ncore)

    for (gen in 1:ngen) {
        mem_optim <- foreach(i = 1:ncore) %dopar% {
            # for (i in 1:ncore) {
            tmp_optim <- run_optim(norm_data, par[1, ], conditions = conditions)
            print(unlist(c(lk = tmp_optim$val, tmp_optim$par)))
        }

        mem_optim <- as.data.frame(do.call(rbind, mem_optim))
        # now compare - keep the best one, allow the others to reproduce

        min <- mem_optim[which.min(mem_optim$lk), ]
        par[1, ] <- min[, -1]
        # for each individual, decide which change or whether to do it
        # which individuals survive?
        # p is the probability of each individual being selected for the next generation
        p <- exp(-(mem_optim[, 1] - min$lk))
        # s is the selected individuals according to the probability p
        s <- sample(1:ncore, ncore - 1, replace = T, prob = p)
        print(mem_optim)[s, -1]
        par[-1, ] <- mem_optim[s, -1] # take the parameter outcomes of each of the optims (remove the RSS values)
        par <- change(par, p_change)
        print(Sys.time())
        print(c(gen, min$lk))
        print(par)
    }

	return(min)
}

