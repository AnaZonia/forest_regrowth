# investigating Heinrich data with optim vs. linear model
library(ggplot2)
#Set the working directory will need to change for each user
setwd("/home/aavila/forest_regrowth/heinrich_poorter/Fig1_data_input/Fig1_data_input/")

########### Load the point data ####################

# Creating a name of data and reading the csv table
dados2<-read.csv2("MCWD_assessment_may2020_CHIRPSv1.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2

dados2 <- dados2[,c(1,3,5)]
colnames(dados2) <- c('age', 'agbd', 'threshold')
mcwd_0<-dados2[dados2$threshold == 0,]
mcwd_180<-dados2[dados2$threshold == -180,]
mcwd_277<-dados2[dados2$threshold == -277,]
mcwd_350<-dados2[dados2$threshold == -350,] 

data <- mcwd_0
# write.csv(data, "/home/aavila/forest_regrowth/heinrich_data.csv")

# here's the result of her fit of the same model with nls:

# Nonlinear regression model
#   model: Corrected_AGB ~ 133 * (1 - exp(-theta2 * age))^theta3
#    data: mcwd_0
#  theta2  theta3 
# 0.02992 1.11620 
#  residual sum-of-squares: 2695

# Number of iterations to convergence: 5 
# Achieved convergence tolerance: 1.19e-06

# when we try her data with optim and the curve:
pars = c(A = 133, age = 0.02992, theta = 1.1162)
pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407) #, sd = 0.05)

G <- function(pars) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

G(pars)

##############

pars = list(A = seq(0, 150, 10), age = seq(0.005, 0.1, 0.005), theta = seq(0.5, 2, 0.05))

tst <- expand.grid(pars)

G <- function(pars) {
  pars[[1]] * (1 - exp(-pars[[2]]*data$age))^pars[[3]]
}

nls <- function(pars) {
  sum((G(pars) - data$agbd)^2)
}

leastsquares <- apply(tst, 1, nls)
negloglik <- apply(tst, 1, NLL)
compare <- cbind(tst, leastsquares, negloglik)
plot(compare$A, compare$negloglik)
plot(compare$age, compare$negloglik)
plot(compare$theta, compare$negloglik)
plot(compare$negloglik)
range(compare$leastsquares[9000:length(compare$leastsquares)])


library(tidyverse)
median_values <- compare %>%
  group_by(A) %>%
  summarize(median_nll = median(negloglik, na.rm = TRUE))
plot(median_values$A, median_values$median_nll)

compare[which.min(compare$negloglik),]

A0 <- compare[compare$A==120,]
plot(A0$age, A0$negloglik)
plot(A0$theta, A0$negloglik)

age0 <- compare[compare$age==0.08,]
plot(age0$A, age0$negloglik)
plot(age0$theta, age0$negloglik)

persp(A0$age, A0$theta, A0$negloglik)

o = optim(par = pars, fn = nls)
o

plot(data$age, data$agbd)
curve( 87.07455636 * (1 - exp(-0.07435007*x))^1.69029407, add = TRUE)

# that's starting from the values given by NLS.
# if we start with values a bit off, it doesn't work anymore:
# starting with:

pars = c(A = 100, age = 1, theta = 1)
o = optim(par = pars, fn = nls)
o

# now the likelihood is 18604, rather than 2612.
# starting with age=0.5, it works again.
# so I guess we gotta try a bunch of different starting values.

# trying NLL with her functional form (no intercept, with a delay term)

#let's investigate why NLL breaks, making zeroes.
pars = c(A = 174.839886, age = 0.301995, theta = 0.532769, sd = 0.310642)

G <- function(pars) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

data$agbd - G(pars)

dnorm(x = data$agbd - G(pars), mean = 0, sd = pars['sd'])

data['agbd'] - G
norm.pdf(data['agbd'] - G)

NLL = function(pars) {
  # Negative log-likelihood
  print(-sum(dnorm(x = data$agbd - G(pars), mean = 0, sd = 0.05, log = TRUE), na.rm = TRUE))
}

#pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407) #, sd = 0.05)

o = optim(par = pars, fn = NLL)
o

G(o$par)

sum((G(o$par) - data$agbd)^2)

  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']

-sum(dnorm(x = data$agbd - G(o$par), mean = 0, sd = pars['sd'], log = TRUE))






# try with random forest as well out of curiosity.