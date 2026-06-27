

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


optim_ga <- function(par, norm_data = norm_data, control = list(), ngen = 50, maxit = 50, p_change = c(.5, .25, .25)) {
    # 	operations of ga - create population, select who survives based on "fitness". Has mutation and cross-over (which in this case includes independent assortment). Can also mix the values of both (e.g, take the midpoint).
    # keep the best performer. Choose the other ones based on their RSSs.

    # can be passed into optim to control the number of iterations per optim run
    control$maxit <- maxit

    options(cores = ncore)

    for (gen in 1:ngen) {
        gen = 1

        # run optim for each individual in the population and store the resulting parameters and RSS
        mem_optim <- foreach(i = 1:ncore, .combine = rbind) %dopar% {
        # for (i in 1:ncore) {
            tmp_optim <- run_optim(norm_data, par[i, ], conditions = conditions)
            print(unlist(c(lk = tmp_optim$val, tmp_optim$par)))
        }
        mem_optim <- as.data.frame(mem_optim)

        # Keep the row with the lowest sum of squares, and reproduce the other rows
        best_row <- mem_optim[which.min(mem_optim$lk), ]
        par[1, ] <- best_row[, -1]

        # which individuals survive?
        # p is the probability of each individual being selected for the next generation.
        # delta shows the difference between the best likelihood and each individual
        # tau avoids situations where the probability is exactly zero for all individuals except the top likelihood one
        delta <- mem_optim[, 1] - min(mem_optim[, 1])
        tau <- sd(mem_optim[, 1]) + 1e-8
        p <- exp(-delta / tau)
        p <- exp(-(mem_optim[, 1] - best$lk))
        
        # s is the selected individuals according to the probability p
        s <- sample(1:ncore, ncore - 1, replace = T, prob = p)
        
        par[-1, ] <- mem_optim[s, -1] # take the parameter outcomes of each of the optims (remove the RSS values)
        
        par <- change(par, p_change)
        
        print(Sys.time())
        print(c(gen, min$lk))
        print(par)
    }

    return(min)
}


