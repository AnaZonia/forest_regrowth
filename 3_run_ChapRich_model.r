################# passing into the model ##################

# pars[2]*tavg+pars[3]*prec+
# central_df$tavg, central_df$prec, 
#  0.001, 0.0001, 

columns = c(other_perennial, other_annual, pasture, soy, forest_age, num_fires, last_burn)

G = function(pars, columns) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  for (i in 1:length(columns)){
    k = pars[2]*other_perennial+pars[3]*pasture+pars[4]*soy+pars[5]*other_annual + pars[6]*other_annual + pars[7]*num_fires^(pars[8]*last_burn)
  }
  # Prediction of the model
  Gmax * (1 - exp(-k*(forest_age)))
}

####################
#finding the right parameters
par0 = c(50, 0.0001,  0.0001,  0.001, 0.001,  0.001, 0.1, 0.1)

Gpred = G(par0, central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

# dat$tavg, dat$prec, 

NLL = function(pars, dat) {
  # Values prediced by the model
  if(pars[9] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars, dat$other_perennial, dat$other_annual, dat$pasture, dat$soy, dat$forest_age, dat$num_fires, dat$last_burn)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = dat$agbd, mean = Gpred, sd = pars[9], log = TRUE))
  #print(pars)
  return(fun)
}

par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.1, 0.1)
NLL(par0, dat=central_df)

o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "CG")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG")
for (i in meth0){
  o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

#  central_df$tavg, central_df$prec, 

pred = G(o$par[1:9], central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

plot(central_df$agbd, pred, abline(0,1))

