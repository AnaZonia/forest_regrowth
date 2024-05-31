library(rstan)
library(tidybayes)
library(fitdistrplus)
library(ggplot2)

data_5 <- read.csv("data/unified_data_5_years.csv")
data_10 <- read.csv("data/unified_data_10_years.csv")
data_15 <- read.csv("data/unified_data_15_years.csv")
all_data <- read.csv("data/unified_data.csv")

#--------------------- Switches

# fit with mean mature forest biomass in the asymptote and see if it works

#-----------------------

# 10k takes about 10 min
tst <- aggregate(agbd ~ age, all_data, median)
data_aggregated <- list(age = tst$age, agbd = tst$agbd, n = nrow(tst))
plot(data_aggregated$age, data_aggregated$agbd)

# data_regular <- list(
#   age = data_10$age,
#   agbd = data_10$agbd,
#   n = nrow(data_10)
# )

iter <- 20000
warmup <- 5000
chains <- 2
# init <- list(
#   list(age_par = 1, theta = 6, B0 = 60, sigma = 0.3),
#   list(age_par = 0.8, theta = 4, B0 = 30, sigma = 0.1)
# )

init <- list(
  list(age_par = 1, theta = 6, B0 = 60, sigma = 0.3, A = 100),
  list(age_par = 0.8, theta = 4, B0 = 30, sigma = 0.1, A = 80)
)

fit_medians <- stan(
  file = "age_agbd.stan", data = data_aggregated,
  iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  init = init,
  control = list(max_treedepth = 12)
)
print(fit_medians)

vals = 30 + 160 * (1 - exp(-(0.01 * tst$age)))

dnorm(
  x = tst$agbd - vals, mean = 0,
  sd = 0.3, log = TRUE
)


# traceplot(fit_medians, par = c("age_par", "A", "B0", "theta", "sigma"))

# plot(fit_medians, par = c("age", "A", "B0", "theta", "sigma"))




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
