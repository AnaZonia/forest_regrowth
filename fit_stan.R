library(rstan)
library(tidybayes)
library(fitdistrplus)
library(ggplot2)

data_5 <- read.csv("data/unified_data_5_years.csv")
data_10 <- read.csv("data/unified_data_10_years.csv")
data_15 <- read.csv("data/unified_data_15_years.csv")
all_data <- read.csv("data/unified_data.csv")

# 10k takes about 10 min
tst <- aggregate(agbd ~ age, data_10, median)
plot(tst$age, tst$agbd)

data_aggregated <- list(age = tst$age, agbd = tst$agbd, n = nrow(tst))

data_regular <- list(
  age = data_10$age,
  agbd = data_10$agbd,
  n = nrow(data_10)
)

# fit_gamma <- fitdist(all_data$age, distr = "gamma", method = "mle")
# fit_weibull <- fitdist(all_data$age, distr = "weibull", method = "mle")
# fit_normal <- fitdist(all_data$age, distr = "norm", method = "mle")
# fit_gamma$aic
# fit_weibull$aic
# fit_normal$aic


iter <- 2000
warmup <- 1000
chains <- 2

fit_medians <- stan(
  file = "age_agbd.stan", data = data_aggregated,
  iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  init = list(
    list(age = 1, theta = 6, log_A = 4, B0 = 60),
    list(age = 0.8, theta = 4, log_A = 6, B0 = 30)
  ),
  control = list(max_treedepth = 12)
)
print(fit_medians)

# write a sequence from 1 to 10
test <- function() {
  # prior check
  fit_medians %>%
    gather_draws(mean_agbd[i], ndraws = 15) |>
    mutate(age = tst$age[i]) |>
    ggplot(aes(x = age, y = .value)) +
    geom_line() +
    facet_wrap(~.draw)
}

traceplot(fit_medians, par = c("age_par", "A", "B0", "theta", "sigma"))


fit_regular_flat <- stan(
  file = "age_agbd.stan", data = data_regular, iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  init = list(
    list(age = 1, theta = 0.5, A = 180, B0 = 60),
    list(age = 0.01, theta = 1, A = 120, B0 = 1)
  ),
  control = list(max_treedepth = 12)
)
print(fit_regular.flat)

data_cwd <- list(
  age = data_10$age, agbd = data_10$agbd,
  b1 = data_10$all, n = nrow(data_10)
)

fit_regular_cwd <- stan(
  file = "age_agbd_cwd.stan", data = data_cwd,
  iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  init = list(
    list(age = 1, theta = 0.5, A = 180, B0 = 60, b1 = 0.01),
    list(age = 0.01, theta = 1, A = 120, B0 = 1, b1 = 1)
  ),
  control = list(max_treedepth = 12)
)

print(fit_regular_cwd)


traceplot(fit_medians, par = c("age", "A", "B0", "theta", "sigma"))
plot(fit_medians, par = c("age", "A", "B0", "theta", "sigma"))

# ---------------- compare models


data <- list(
  ages = import$age, agbds = import$agbd, last_LU = import$last_LU, b1 = import$b1,
  # num_fires_before_regrowth = import$num_fires_before_regrowth, fallow = import$fallow,
  n = nrow(import)
)

iter <- 2000
warmup <- 1000
chains <- 2

fit.stan <- stan(
  file = "model.stan", data = data, iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  init = list(
    list(age = 1, b1 = 1, last_LU = 1, theta = 0.5, A = 80, B0 = 5),
    list(age = 0.01, b1 = 0.01, last_LU = 0.01, theta = 1, A = 100, B0 = 1)
  ),
  control = list(max_treedepth = 12)
)

data <- list(
  ages = import$age, agbds = import$agbd,
  last_LU = import$last_LU, b1 = import$b1, n = nrow(import)
)

fit.stan <- stan(
  file = "model.stan", data = data, iter = iter, warmup = warmup,
  chains = 1, cores = 4,
  init = list(list(
    age = 0.5, b1 = 0.5,
    theta = 1, A = 80, B0 = 1
  )),
  control = list(max_treedepth = 12)
)

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

