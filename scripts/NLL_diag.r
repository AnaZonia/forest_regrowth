# ----------------------------
# Investigando qual foi o problema com o classificador de negative log-likelihood.
#------------------------------

library(ggplot2)

data_init <- readRDS('./forest_regrowth/data/santoro_cwd_fire_LU.rds')
data <- data_init[,1:2]
head(data)
nrow(data)

remove_outliers <- function(data){
  means <- aggregate(agbd ~ age, data, median)
  colnames(means) <- c('age', 'mean')
  sd <- aggregate(data$agbd, by = list(data$age), FUN = sd)
  colnames(sd) <- c('age', 'sd')
  data <- merge(data, means, by = 'age')
  data <- merge(data, sd, by = 'age')
  data[abs(data$agbd-data$mean) < 0.25*data$sd, ]
}

NLL = function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 1){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  if (pars['sd'] < 0){
    return(-Inf) 
  }
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}


nls <- function(pars, data) {
  result <- sum((G(pars, data) - data$agbd)^2)
  ifelse(result == 0, -Inf, result)
}


G <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

data_reduced <- remove_outliers(data)

pars <- c(A = 100, age = 0.5, theta = 1.5)

o_nls <- optim(pars, nls, data=data_reduced)
o_nls

# weirdly, here the default seems to outperform BFGS. Not sure why this is different
# from the example in exploring_heinrich

plot(data_reduced$age, data_reduced$agbd)
lines(unique(G(o_nls$par, data = data_reduced)), col = 'red')

pars <- c(A = 100, age = 0.5, theta = 1.5, sd = 0.05)

o <- optim(pars, NLL, data=data_reduced)
o
lines(unique(G(o$par, data = data_reduced)), col = 'blue')

result <- G(o$par, data = data_reduced)
rss_nll <- sum((data_reduced$agbd - result)^2)

# o_nls$value < rss_nll, so log likelihood is not finding the best parameters once more.

################### bare bones example, with heinrich data.
# this is mcwd_0 from heinrich
data_example <- data.frame(age = 1:32, agbd = c(19.0, 0.0, 5.0, 5.5, 14.0, 12.5, 15.0, 21.0, 28.0, 31.0, 26.0, 34.5, 51.5, 43.5, 52.5, 47.0, 42.5, 57.0, 48.0, 44.0, 52.5, 76.0, 59.0, 
79.0, 62.0, 77.5, 49.5, 58.5, 77.0, 55.5, 77.0, 89.5))

pars <- c(A = 100, age = 0.5, theta = 1.5)

o_nls <- optim(pars, nls, data=data_example, method = 'BFGS')
o_nls # o_nls$value = 2612.017

#####################
# note that
o_nls <- optim(pars, nls, data=data_example)
o_nls # o_nls$value = 2628.707
####################

pars <- c(A = 100, age = 0.5, theta = 1.5, sd = 0.5)
o_nll <- optim(pars, NLL, data=data_example)
o_nll # o_nll$value = 2628.707
sum((data_example$agbd - G(o_nll$par, data_example))^2) #2638.79

#so, it's not reaching the correct value. why?

# inputting parameters from nls to NLL:
o <- optim(c(o_nls$par, sd = 0.5), NLL, data=data_example)
o

sum((data_example$agbd - G(o$par, data = data_example))^2)
# 2612.054 - slightly larger than the value gotten from nls, but it's negligible.

# trying "repuffing" NLL:

for(i in 1:30){
  o <- optim(o$par, NLL, data=data_example)
  print(o$value)
  print(o$par)
}
sum((data_example$agbd - G(o$par, data = data_example))^2)
# 2612.003 - absolute minimum! repuffing worked!

# what if we try to run it only once with sd fixed?

NLL = function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 1){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = 1, log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

pars <- c(A = 100, age = 0.5, theta = 1.5)

o <- optim(pars, NLL, data=data_example)
o
sum((data_example$agbd - G(o$par, data = data_example))^2)
#[1] 2612.003 - absolute minimum as well. the issue is the standard deviation perhaps?

###############################################
# now trying that with my own data:
pars <- c(A = 100, age = 0.5, theta = 1.5, sd = 0.05)

o <- optim(pars, NLL, data=data_reduced)
o
NLL(c(o_nls$par, sd = 11.55545802), data_reduced)

for(i in 1:30){
  o <- optim(o$par, NLL, data=data_reduced)
  print(o$value)
  print(o$par)
}

values <- c()
for(i in 1:30){
  new_par <- o$par + c(runif(1,-5,5), runif(1,0.1,0.5), runif(1,-0.1,1), runif(1,-1,1))
  o <- optim(new_par, NLL, data=data_reduced)
  values <- c(values, o$value)
  print(o$value)
  print(o$par)
}

# okay - so repuffing won't work with my own data, as it goes even further away.

############################    SD     ###############################3

# Now, the issue comes when sd is fit. What if it's held fixed?

NLL = function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 1){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = 1, log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

pars <- c(A = 100, age = 0.5, theta = 1.5)

o <- optim(pars, NLL, data=data_reduced)
o
sum((data_reduced$agbd - G(o$par, data = data_reduced))^2)

# okay - so now we get a better value with NLL if sd is held fixed.

# let me try with a small variation for sd:

NLL = function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 1){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  if (pars['sd'] < 0){
    return(-Inf) 
  }
  if (pars['sd'] > 10){
    return(-Inf) 
  }
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

pars <- c(A = 100, age = 0.5, theta = 1.5, sd = 0.5)

o <- optim(pars, NLL, data=data_reduced)
o
sum((data_reduced$agbd - G(o$par, data = data_reduced))^2) #62618695
lines(unique(G(o$par, data = data_reduced)), col = 'green')
