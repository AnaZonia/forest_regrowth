library(cmdstanr)

setwd("C:/Users/anaca/Desktop/forest_regrowth/")

import <- read.csv("data/unified_data_10_years.csv")

iter <- 2000
warmup <- 1000

data <- list(ages = import$age, agbds = import$agbd,
             last_LU = import$last_LU, b1 = import$b1, n = nrow(import))

mod_bern <- cmdstan_model(stan_file = "model.stan")
fit_bern <- mod_bern$sample(data = data, iter_sampling = iter,
                            iter_warmup = warmup, chains = 2, cores = 4)