#Genetic Algorithm integrated with optim to reduce the local minima problem by Brian Leung - March 2025.

library(doParallel)
library(foreach)

# ncore = 25

registerDoParallel()


growth_curve <- function(par, data, lag = 0) {

    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(par), c(non_data_pars, climatic_pars))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k

    k <- rep(par[["k0"]], nrow(data)) * data[["age"]]

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    return(par[["B0"]] + (data[["nearest_biomass"]] - par[["B0"]]) * (1 - exp(-k))^par[["theta"]])
}


# growth_curve(ini_par[1,], norm_data)

wrapper<-function(par,le_func, data, conditions)
{
    if ("B0" %in% names(par)) {
        conditions <- c(conditions, list('par[["B0"]] < 0'))
    } else {
        conditions <- c(conditions, list('par[["lag"]] < 0'))
    }

	return(optim(par, le_func, data=data, conditions = conditions))
}

likelihood <- function(par, data, conditions) {

    if ("lag" %in% names(par)) {
        growth_curves <- growth_curve(par, data, par[["lag"]])
        residuals <- growth_curves - data$biomass
        result <- mean(residuals^2)
    } else {
        result <- mean((growth_curve(par, data) - data$biomass)^2)
    }

    print(par[["B0"]])
    print(result)
    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
    return(-Inf)
    } else if (is.na(result) || result == 0) {
    return(-Inf)
    } else {
    return(result)
    }

}



change=function(par,p_change)
{
	tp=sample(1:3,ncore,replace=T,prob=p_change)
	for(i in 2:ncore)
	{
		if(tp[i]==1){ #mutation
			par[i,]=par[i,] * (.8+ .4*runif(ncol(par)))
		}else if(tp[i]==2){ #cross over
			s=sample(1:ncore,1)#find other parent
			t=sample(1:ncol(par),1) #choose 1 trait to swap
			par[i,t]=par[s,t]
			t=sample(1:ncol(par),1) #also mutate one other trait, so that make sure don't have exact duplicates
			par[i,t]=par[i,t] * (.8+ .4*runif(1))

		}else{ #co-dominance
			s=sample(1:ncore,1)#find other parent
			par[i,]=(par[i,]+par[s,])/2
			t=sample(1:ncol(par),1) #also mutate one other trait, so that make sure don't have exact duplicates
			par[i,t]=par[i,t] * (.8+ .4*runif(1))
		}
	}
	return(par)

}



optim_gm<-function(par,le_func, df=NULL,wrap_optim=wrapper, control=list(), ngen=50,maxit=50, ncore=15, p_change=c(.5,.25,.25)) #p_change - mutation,cross over, co-dominance - issue with none, is that all might be identical 
{
#	operations of ga - create population, select who survives based on "fitness". Has mutation and cross-over (which in this case includes independent assortment). Can also mix the values of both (e.g, take the midpoint).
#keep the best performer. Choose the other ones based on their likelihoods. 
	control$maxit=maxit
	options(cores=ncore)
		
	for(gen in 1:ngen)
	{
        for (i in 1:ncore)
            # mem_optim=foreach(i = 1:ncore) %dopar% #simplify
            {
                tmp_optim=wrap_optim(par[i,],le_func, control, df, conditions)
                # return(unlist(c(lk=tmp_optim$val,tmp_optim$par)))
                print(unlist(c(lk=tmp_optim$val,tmp_optim$par)))
            }
		mem_optim=as.data.frame(do.call(rbind,mem_optim))
		#now compare - keep the best one, allow the others to reproduce
#		par2=par
		min=mem_optim[which.min(mem_optim$lk),]
		par[1,]=min[,-1]
		#for each individual, decide which change or whether to do it
		#which individuals survive?
		p=exp(-(mem_optim[,1]-min$lk))
		s=sample(1:ncore,ncore-1,replace=T,prob=p) #these are the new individuals
        print(mem_optim)[s,-1]
		par[-1,]=mem_optim[s,-1] #take the parameter outcomes of each of the optims (remove the likelihood values)
		par=change(par,p_change)
		print(Sys.time())
		print(c(gen,min$lk))
        print(par)
	}
	return(min)
}

ncore = 15

ini_par = data.frame(B0 = 80, k0 = 1, theta = 2)

for(j in 2:ncore)
{
    ini_par[j,]=c(ini_par[1,"B0"]*(1.5*runif(1)+.5), ini_par[1,"k0"]*(1.5*runif(1)+.5), ini_par[1, "theta"]*(1.5*runif(1)+0.5)) #the null model does not have b (which is the parameter of interest in our case)
}
ini_par

# optim_gm(ini_par, likelihood, norm_data, wrapper)
conditions <- list('par[["theta"]] > 10', 'par[["theta"]] < 0', 'par[["k0"]] < 0')

for (i in 1:ncore)
    # mem_optim=foreach(i = 1:ncore) %dopar% #simplify
    {
        print(i)
        i = 1
        ini_par[1,] = data.frame(B0 = 120, k0 = 1, theta = 1)

        if ("B0" %in% names(ini_par[i,])) {
            conditions_iter <- c(conditions, list('par[["B0"]] < 0'))
        } else {
            conditions_iter <- c(conditions, list('par[["lag"]] < 0'))
        }
        print(conditions_iter)

	# return(optim(par, le_func, data=data, conditions = conditions))
        tmp_optim = optim(ini_par[i,], likelihood, data=norm_data, conditions=conditions_iter)
        # return(unlist(c(lk=tmp_optim$val,tmp_optim$par)))
        # print(tmp_optim)
    }# head(norm_data$biomass)


#*** running the optim_gm
#**	set initialization parameters 

#	ini_par=data.frame(V1 , V2)


#** iterate through remaining cores allocated to project
#** set DIFFERENT starting values 

#	for(j in 2:ncore)
#	{
#		
#		ini_par[j,]=c(ini_par[1,"V1"]*(1.5*runif(1)+.5), ini_par[1,"V2"]*(1.5*runif(1)+.5)) #the null model does not have b (which is the parameter of interest in our case)
#	}



#** Note that if are comparing this against a previous model (e.g., simpler model), can set the best optim run parameters as one of the starting values

#** pass it parameter matrix to fit, the evaluation function, and a number of optional variables

#** optional variables include the data, the wrapper function around optim, control variables used in optim, number of generations, number of iterations for optim, number of cores to use, and the probability of the different "evolutionary" operations

#	ret=optim_gm(ini_par,le_func, df,wrapper,control, ngen,maxit) 





