# investigating Heinrich data with optim vs. linear model
library(ggplot2)

########### Load the point data ####################

# Creating a name of data and reading the csv table
dados2<-read.csv2("./forest_regrowth/heinrich_poorter/Fig1_data_input/Fig1_data_input/MCWD_assessment_may2020_CHIRPSv1.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2

dados2 <- dados2[,c(1,3,5)]
colnames(dados2) <- c('age', 'agbd', 'threshold')
mcwd_0<-dados2[dados2$threshold == 0,]
mcwd_180<-dados2[dados2$threshold == -180,]
mcwd_277<-dados2[dados2$threshold == -277,]
mcwd_350<-dados2[dados2$threshold == -350,] 

plot_compare <- function(title, pars, data){
  par(cex.lab = 1.5)
  plot(data$age, data$agbd, main = title, xlab = "Age", ylab = "Biomass")
  lines(G(pars, data), col = "red", lwd = 2, legend.text = "Heinrich")
  lines(G(o$par, data), col = "blue", lwd = 2, legend.text = "Optim")
  legend("topright", legend = c("Heinrich", "Optim"), col = c("red", "blue"), lwd = 4, cex = 1.2)
}
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
pars_0 = c(A = 133, age = 0.02992, theta = 1.1162) # this is the value she found with nls function
# and with RSS = 2695.

G <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

# Indeed, starting out with her parameters as initial parameters we get better results:

nls <- function(pars, data) {
  result = sum((G(pars, data) - data$agbd)^2)
  ifelse(result == 0, -Inf, result)
}

o = optim(par = pars_0, fn = nls, data = mcwd_0)
o

# And this leaves us with
# pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407)
# with RSS 2612.003
# when running with method='BFGS', we get no deviation from initial values
# and RSS is 2695.041.

# now, the issue is that optim default method seems very sensitive to initial conditions.
# if we start with values a bit off, it doesn't work anymore:
# starting with:

pars = c(A = 100, age = 1, theta = 1)
o = optim(par = pars, fn = nls, data = mcwd_0)
o

# now the RSS is 18604, rather than 2612.
# starting with age=0.5, and the default method, we have once more acceptable predictions.

pars = c(A = 50, age = 2, theta = 0.5)
o = optim(par = pars, fn = nls, data = mcwd_0)
o

plot_compare('Low water stress', pars_0, mcwd_0)

# now, things get weird when we start comparing other CWD intervals.

############# to visualize

dfs <- list(mcwd_0, mcwd_180, mcwd_277, mcwd_350)
colors <- c("red", "blue", "green", "purple")
plot(NULL, xlim = range(mcwd_0$age), ylim = range(mcwd_0$agbd), 
     xlab = "X-axis Label", ylab = "Y-axis Label", main = "Combined Data")

for (i in 1:4) {
  points(dfs[[i]]$age, dfs[[i]]$agbd, col = colors[i], pch = i) # Change 'points' to 'lines' if needed
}

legend("topright", legend = c("1st quartile", "2nd quartile", "3rd quartile", "4th quartile"), col = colors, pch = 1:4)

#------------------

lines(G(pars_0, mcwd_0), col = "red", lwd = 2, legend.text = "Heinrich_0")
lines(G(o$par, mcwd_0), col = "red", lwd = 2, legend.text = "Optim_0", lty = 2)

pars_180 = c(A = 119, age = 0.02766, theta = 1.02515) # RSS 2002
pars <- c(A = 100, age = 0.5, theta = 1.5)
o = optim(par = pars, fn = nls, data = mcwd_180, method='BFGS')
o # RSS 1642
 
lines(G(pars_180, data = mcwd_180), col = "blue", lwd = 2, legend.text = "Heinrich_180")
lines(G(o$par, data = mcwd_180), col = "blue", lwd = 2, legend.text = "Optim_180", lty = 2)

pars_277 = c(A = 97, age = 0.0244, theta = 0.8920) # RSS 362
pars <- c(A = 100, age = 0.5, theta = 1.5)
o = optim(par = pars, fn = nls, data = mcwd_277, method='BFGS')
o # RSS 290.6
 
lines(G(pars_277, data = mcwd_277), col = "green", lwd = 2, legend.text = "Heinrich_277")
lines(G(o$par, data = mcwd_277), col = "green", lwd = 2, legend.text = "Optim_277", lty = 2)

pars_350 = c(A = 88.5, age = 0.0162, theta = 0.8019) # RSS 379.9
pars <- c(A = 100, age = 0.5, theta = 1.5)
o = optim(par = pars, fn = nls, data = mcwd_350, method='BFGS')
o # RSS 279.3
 
lines(G(pars_350, data = mcwd_350), col = "purple", lwd = 2, legend.text = "Heinrich_350")
lines(G(o$par, data = mcwd_350), col = "purple", lwd = 2, legend.text = "Optim_350", lty = 2)

##############################################################

# I then decided to try NLL and see what happens.
# Ideally we'd find pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407)
# with RSS of 2612 or better.

# starting from the expected ideal:

pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407, sd = 0.5) # this is the value she found with nls function

NLL = function(pars, data) {
  # Negative log-likelihood
  print(-sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE))
}

o = optim(par = pars, fn = NLL, data = mcwd_0)
o

# this gives zero likelihood.
# pars = c(A = 91.428284, age = -8.633106, theta = 6.044022, sd = 4.853728 ) # using the output that gave zero likelihood to investigate
# dnorm(x = mcwd_0$agbd - G(pars, mcwd_0), mean = 0, sd = pars['sd'])
# it is finding zeroes when it tests negative values for age.

# try to correct to return -Inf for NaN,
NLL = function(pars, data) {
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

o = optim(par = pars, fn = NLL, data = mcwd_0)
o

# perfect - initial values were maintained, with sd = 9, and NLL = 115.84.
# testing with different initial values
pars = c(A = 100, age = 1, theta = 1, sd = 0.5)
o = optim(par = pars, fn = NLL, data = mcwd_0)
o

# a different value was found, with NLL = 147.25 and
#       A      age    theta       sd 
# 43.93502 34.03816 12.01732 24.07951

# okay - so local minima are an issue. I can try to restrain the initial values to a more reasonable range.
pars = c(A = 100, age = 1, theta = 1, sd = 0.5)

G <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

NLL = function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 10){
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

o = optim(par = pars, fn = NLL, data = mcwd_0)
o
# now it converges to another local minimum:
# $par
#            A          age        theta           sd 
# 140.33915251   0.02707856   1.08620319   9.21236535 

# $value
# [1] 116.4026

head(G(o$par, data = mcwd_0))
plot(mcwd_0$age, mcwd_0$agbd)
lines(G(o$par, data = mcwd_0), col = "blue", lwd = 2)

# I tried running the same with ax (bayesian optimization) on Python, also
# got stuck on local minima. the best it could do was NLL = 120, so worse than optim.


# G <- function(pars, data) {
#   k = pars['age']*data$age+pars['pasture']*data$pasture+pars['other_annual']*data$other_annual
#   pars['A'] * (1 - exp(-k))^pars['theta']
# }

# pars <- c(A = 100, age = 0.5, theta = 1.5, pasture = 0.05, other_annual = 0.05)

# o <- optim(pars, nls, data=data_LU, method = 'BFGS')
# o

# looking into land use change

# data <- data_init[data_init$last_LU %in% c(15, 41, 48),]
# data$last_LU <- factor(data$last_LU)
# dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
# names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
# data <- cbind(data, dummy_LU)

# k = o$par['age']*1+o$par['pasture']*1#+o$par['other_annual']*data_LU$other_annual
# k
# result = (1 - exp(-k))^o$par['theta']
# result
# range(result)
# var <- o$par['theta']
# (o$par['age']+o$par['pasture'])^var
# (pars['age']+pars['pasture'])^pars['theta']
