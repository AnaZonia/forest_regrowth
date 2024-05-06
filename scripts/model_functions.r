####################################################################
##################### Functions ######################
####################################################################


# ----------------------------------
# Growth Equations
# ----------------------------------

G_3par <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

G_4par <- function(pars, data) {
  pars['B0'] + pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

# ----------------------------------
# Fitting Functions
# ----------------------------------

# Negative log-likelihood
NLL = function(pars, data, G) {
  if (pars['age'] < 0 || pars['age'] > 10 ||
      pars['theta'] > 10 || pars['theta'] < 0 ||
      pars['sd'] < 0 || pars['sd'] > 10){
    return(-Inf)
  }
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0,
                sd = pars['sd'], log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

# Nonlinear Least Squares
nls <- function(pars, data, G) {
  if (#pars['age'] < 0 || pars['age'] > 5 ||
      pars['theta'] > 10 || pars['theta'] < 0 || pars['B0'] < 0) {
    return(-Inf)
  }
  result = sum((G(pars, data) - data$agbd)^2)
  ifelse(result == 0, -Inf, result)
}

# repuff NLL output
repuff <- function(pars, data, repeats) {
  for(i in 1:repeats){
    o <- optim(o$par, NLL, data)
    print(o$value)
    print(o$par)
  }
  return('RSS: ' + sum((data$agbd - G(o$par, data))^2))
}
