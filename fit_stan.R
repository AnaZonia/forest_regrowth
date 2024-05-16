library(rstan)

setwd("C:/Users/anaca/Desktop/forest_regrowth/")

data <- read.csv('./data/amazon_df_sample_10million.csv')

tst <- aggregate(age~agbd, data, median)


data = list(df=tst, predictors=ncol(tst), pixels=nrow(tst))

iter = 2000
warmup=1000
chains=2


fit.stan <- stan(file=model.stan, data=data,iter=iter, warmup = warmup, chains=chains, cores=3,
                 init = list(list()))