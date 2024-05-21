library(rstan)

setwd("C:/Users/anaca/Desktop/forest_regrowth/")

land_use_5_years <- read.csv('data/unified_data_5_years.csv') # 44751 rows
land_use_10_years <- read.csv('data/unified_data_10_years.csv') # 22952 rows
land_use_15_years <- read.csv('data/unified_data_15_years.csv') # 19845 rows
land_use <- read.csv('data/unified_data.csv')

# 10k takes about 10 min
tst <- aggregate(agbd~age, land_use_10_years, median)
# plot(tst$age, tst$agbd)

data = list(age=tst$age, agbd=tst$agbd, n=nrow(tst))

iter = 2000
warmup = 1000
chains = 2

fit_medians <- stan(file="model.stan", data = data, iter = iter, warmup = warmup,
                 chains = chains, cores = 4,
                 init = list(list(age = 1, theta = 0.5, A = 80, B0 = 5),
                             list(age = 0.01, theta = 1, A = 100, B0 = 1)),
                 control=list(max_treedepth=12))

print(fit_medians)
traceplot(fit_medians, par=c("age", "A", "B0", "theta", "sigma"))
plot(fit_medians, par=c("age", "A", "B0", "theta", "sigma"))

ext_fit <- extract(fit.stan)
for (val in names(ext_fit)){
  print(mean(ext_fit[]))
}

G_3par <- function(data) {
  mean(ext_fit['B0'][[1]]) + mean(ext_fit['A'][[1]]) * (1 - exp(-mean(ext_fit['age'][[1]])*data$age))^mean(ext_fit['theta'][[1]])
}

pred = G_3par(data)
plot(data$agbds, pred)
abline(c(0,1))


# ---------------- compare models


data = list(ages=import$age, agbds=import$agbd, last_LU=import$last_LU, b1 = import$b1,
            #num_fires_before_regrowth=import$num_fires_before_regrowth, fallow=import$fallow, 
            n=nrow(import))

iter = 2000
warmup = 1000
chains = 2

fit.stan <- stan(file="model.stan", data = data, iter = iter, warmup = warmup,
                 chains = chains, cores = 4,
                 init = list(list(age = 1, b1 = 1, last_LU = 1, 
                                  #num_fires_before_regrowth = 1, fallow = 1,
                                  theta = 0.5, A = 80, B0 = 5),
                             list(age = 0.01, b1 = 0.01, last_LU = 0.01,
                                  #num_fires_before_regrowth = 0.01, fallow = 0.01,
                                  theta = 1, A = 100, B0 = 1)),
                 control=list(max_treedepth=12))

data = list(ages=import$age, agbds=import$agbd, last_LU=import$last_LU, b1 = import$b1, n=nrow(import))
            
fit.stan <- stan(file="model.stan", data = data, iter = iter, warmup = warmup,
                 chains = 1, cores = 4,
                 init = list(list(age = 0.5, b1 = 0.5,
                                  theta = 1, A = 80, B0 = 1)),
                 control=list(max_treedepth=12))

print(fit.stan)
# traceplot(fit.stan, par=c("age", "A", "B0", "theta", "sigma"))
# plot(fit.stan, par=c("age", "A", "B0", "theta", "sigma"))
# 
# ext_fit <- extract(fit.stan)
# for (val in names(ext_fit)){
#   print(mean(ext_fit[]))
# }
# 
# G_3par <- function(data) {
#   mean(ext_fit['B0'][[1]]) + mean(ext_fit['A'][[1]]) * (1 - exp(-mean(ext_fit['age'][[1]])*data$age))^mean(ext_fit['theta'][[1]])
# }
# 
# pred = G_3par(data)
# plot(data$agbds, pred)
# abline(c(0,1))


library(stats)

cor(import2)
