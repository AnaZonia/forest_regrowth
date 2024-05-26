library(rstan)
library(tidybayes)
library(fitdistrplus)
library(ggplot2)

land_use_5_years <- read.csv("data/unified_data_5_years.csv")
land_use_10_years <- read.csv("data/unified_data_10_years.csv")
land_use_15_years <- read.csv("data/unified_data_15_years.csv")
land_use <- read.csv("data/unified_data.csv")

# 10k takes about 10 min
tst <- aggregate(agbd ~ age, land_use_10_years, median)
plot(tst$age, tst$agbd)

data_aggregated <- list(age = tst$age, agbd = tst$agbd, n = nrow(tst))
data_regular <- list(age = land_use_10_years$age,
                     agbd = land_use_10_years$agbd, n = nrow(land_use_10_years))

data <- land_use$age
fit_gamma <- fitdist(land_use$age, distr = "gamma", method = "mle")
fit_weibull <- fitdist(land_use$age, distr = "weibull", method = "mle")
fit_normal <- fitdist(land_use$age, distr = "norm", method = "mle")
fit_gamma$aic
fit_weibull$aic
fit_normal$aic

iter <- 2000
warmup <- 1000
chains <- 2

fit_medians <- stan(file = "age_agbd.stan", data = data_aggregated,
                    iter = iter, warmup = warmup,
                    chains = chains, cores = 4,
                    init = list(list(age = 1, theta = 0.5, log_A = 4, B0 = 60),
                                list(age = 0.01, theta = 1, log_A = 6, B0 = 1)),
                    control = list(max_treedepth = 12))
print(fit_medians)

# prior check
fit_medians %>%
  gather_draws(mean_agbd[i], ndraws = 15) |>
  mutate(age = tst$age[i]) |>
  ggplot(aes(x = age, y  = .value)) +
  geom_line() +
  facet_wrap(~.draw)


traceplot(fit_medians, par = c("age_par", "A", "B0", "theta", "sigma"))

fit_regular.flat <- stan(file = "age_agbd.stan", data = data_regular, iter = iter, warmup = warmup,
                       chains = chains, cores = 4,
                      init = list(list(age = 1, theta = 0.5, A = 180, B0 = 60),
                                  list(age = 0.01, theta = 1, A = 120, B0 = 1)),
                       control = list(max_treedepth = 12))
print(fit_regular.flat)

data_cwd = list(age = land_use_10_years$age, agbd =  land_use_10_years$agbd, b1 =  land_use_10_years$all, n = nrow(land_use_10_years))

fit_regular_cwd <- stan(file = "age_agbd_cwd.stan", data = data_cwd,
                        iter = iter, warmup = warmup,
                        chains = chains, cores = 4,
                        init = list(list(age = 1, theta = 0.5, A = 180, B0 = 60, b1 = 0.01),
                                      list(age = 0.01, theta = 1, A = 120, B0 = 1, b1 = 1)),
                        control = list(max_treedepth = 12))

print(fit_regular_cwd)


traceplot(fit_medians, par = c("age", "A", "B0", "theta", "sigma"))
plot(fit_medians, par = c("age", "A", "B0", "theta", "sigma"))

G_3par <- function(pars, data) {
  pars["B0"] + pars["A"] * (1 - exp(-pars["age"]*data$age))^pars["theta"]
}

pars <- c(B0 = 40, A = 80, theta = 1.5, age = 1) # intercept, asymptote, shape term, standard deviation

# Nonlinear Least Squares
nls <- function(pars, data, G) {
  if (pars["age"] < 0 || pars["age"] > 5 ||
    pars["theta"] > 10 || pars["theta"] < 0 || pars["B0"] < 0) {
    return(-Inf)
  }
  result = sum((G(pars, data) - data$agbd)^2)
  ifelse(result  =  =  0, -Inf, result)
}

o = optim(pars, fn = nls, data = data, G = G_3par)
o

pred = G_3par(data)
plot(data$agbds, pred)
abline(c(0,1))


# ---------------- compare models


data = list(ages = import$age, agbds = import$agbd, last_LU = import$last_LU, b1 = import$b1,
            #num_fires_before_regrowth = import$num_fires_before_regrowth, fallow = import$fallow, 
            n = nrow(import))

iter = 2000
warmup = 1000
chains = 2

fit.stan <- stan(file = "model.stan", data = data, iter = iter, warmup = warmup,
                 chains = chains, cores = 4,
                 init = list(list(age = 1, b1 = 1, last_LU = 1, 
                                  #num_fires_before_regrowth = 1, fallow = 1,
                                  theta = 0.5, A = 80, B0 = 5),
                             list(age = 0.01, b1 = 0.01, last_LU = 0.01,
                                  #num_fires_before_regrowth = 0.01, fallow = 0.01,
                                  theta = 1, A = 100, B0 = 1)),
                 control =  list(max_treedepth = 12))

data = list(ages = import$age, agbds = import$agbd, last_LU = import$last_LU, b1 = import$b1, n = nrow(import))
            
fit.stan <- stan(file = "model.stan", data = data, iter = iter, warmup = warmup,
                 chains = 1, cores = 4,
                 init = list(list(age = 0.5, b1 = 0.5,
                                  theta = 1, A = 80, B0 = 1)),
                 control =  list(max_treedepth = 12))

print(fit.stan)
# traceplot(fit.stan, par = c("age", "A", "B0", "theta", "sigma"))
# plot(fit.stan, par = c("age", "A", "B0", "theta", "sigma"))
# 
# ext_fit <- extract(fit.stan)
# for (val in names(ext_fit)){
#   print(mean(ext_fit[]))
# }
# 
# G_3par <- function(data) {
#   mean(ext_fit["B0"][[1]]) + mean(ext_fit["A"][[1]]) * (1 - exp(-mean(ext_fit["age"][[1]])*data$age))^mean(ext_fit["theta"][[1]])
# }
# 
# pred = G_3par(data)
# plot(data$agbds, pred)
# abline(c(0,1))

library(stats)

