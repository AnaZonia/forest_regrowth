library(rstan)
library(tidybayes)
library(fitdistrplus)
library(ggplot2)
library(shinystan)

data_5 <- read.csv("data/unified_data_5_years.csv")
data_10 <- read.csv("data/unified_data_10_years.csv")
data_15 <- read.csv("data/unified_data_15_years.csv")
all_data <- read.csv("data/unified_data.csv")

#--------------------- Switches

# fit with mean mature forest biomass in the asymptote and see if it works

#-----------------------
# install.packages("httpgd")
# 10k takes about 10 min
tst <- aggregate(agbd ~ age, all_data, median)
data_aggregated <- list(age = tst$age, agbd = tst$agbd, n = nrow(tst))
plot(data_aggregated$age, data_aggregated$agbd)

# Plot the data
vals <- 30 + 160 * (1 - exp(-(0.01 * c(1:35))))^0.3 + rnorm(n = 35, sd = 5)
plot(c(1:35), vals)
fake_data <- list(age=c(1:35), agbd=vals, n=35)

# data_regular <- list(
#   age = data_10$age,
#   agbd = data_10$agbd,
#   n = nrow(data_10)
# )

iter <- 2000
warmup <- 1000
chains <- 2

# nelder mead in stan?
# transforming things so that all follows more or less a normal around zero is ideal

fit_medians2 <- stan(
  file = "age_agbd.stan", data = fake_data,
  iter = iter, warmup = warmup,
  chains = chains, cores = 4,
  control = list(max_treedepth = 12)
)
print(fit_medians2)

launch_shinystan(fit_medians2)

curve(dnorm(x, mean = 40, sd = 10), from = 0, to = 80)

curve(dexp(x, rate = 0.1), from = 0, to = 100)

curve(dbeta(x, shape1 = 3, shape2 = 5),
  from = 0, to = 1,
  main = "Beta Distribution", xlab = "x", ylab = "Density"
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



tst <- read.csv("./data/all_land_use.csv")
nrow(tst)
hist(tst$age)
