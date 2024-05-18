library(rstan)

setwd("C:/Users/anaca/Desktop/forest_regrowth/")

data <- read.csv('./data/amazon_df_sample_10million.csv')

tst <- aggregate(age~agbd, data, median)


data = list(ages=tst$age, agbds=tst$agbd, n=nrow(tst))

iter = 2000
warmup=1000
chains=2

fit.stan <- stan(file="model.stan", data = data, iter = iter, warmup = warmup,
                 chains = chains, cores = 3,
                 init = list(list(age = 1, theta = 0.5, A = 80, B0 = 5),
                             list(age = 0.01, theta = 1, A = 100, B0 = 1)),
                 control=list(max_treedepth=12))