####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################
library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")

########################
# SWITCHES
data_prec_temp = TRUE
scaled = FALSE

minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

if (data_prec_temp == TRUE){
  data <- readRDS('santoro_ESA_alldata.rds')
  data <- data[data$last_LU %in% c(15, 41, 48),]
  data$last_LU <- factor(data$last_LU)
  dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
  names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
  data_raw <- data[,-7]
  if(scaled==TRUE){
    data_scaled <- cbind(agbd=data_raw$agbd, scale(data_raw[,2:ncol(data_raw)]))
    data <- cbind(data_scaled, dummy_LU)
  }else{
    data_maxmin <- cbind(agbd=data_raw$agbd, minMax(data_raw[,2:ncol(data_raw)]))
    data <- cbind(data_maxmin, dummy_LU)
  }
}else{
  data <- readRDS('total_aet_cwd.rds')
}

######################################################################
#################        passing into the model     ##################
######################################################################

sds <- aggregate(agbd ~ age, data, sd)
means <- aggregate(agbd ~ age, data, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'agbd', 'sd')
data <- sum_stats

fit <- lm(sum_stats$agbd ~ sum_stats$age)
summary(fit)

ggplot(sum_stats, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))

# Predict values using the model
predicted_values <- predict(fit, as.data.frame(data$age))

mean((data$agbd - predicted_values)^2)
mean((data$agbd - pred)^2)


G <- function(pars) {
  E = pars['prec'] * data$age # + pars['temp'] * data$temp
  #LU = pars['total_fires'] * data$total_fires + pars['ts_fire'] * data$ts_fire + pars['pasture'] * data$pasture + 
  #    pars['other_perennial'] * data$other_perennial + pars['other_annual'] * data$other_annual 
  k = E
  pars['B_0'] + pars['A'] * (1 - exp(-k))
}


#pars = c(B_0 = 10, A = 100, temp =- 0.002, prec = 0.000005, total_fires = 0.05, ts_fire = 0.05, pasture = 0.05, other_perennial = 0.05, other_annual = 0.05, sd = 0.05)
pars = c(B_0 = 10, A = 100, prec = 2000, sd = 0.05)
Gpred <- G(pars)
Gpred
# there are a lot of repeated values of cwd and aet with varying agbd values. Note that the scales are very different!
# just for test, we should try this with medians too.

NLL = function(pars) {
  if(pars['sd'] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  # Negative log-likelihood 
  print(-sum(dnorm(x = data$agbd - G(pars), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE))
}

o = optim(par = pars, fn = NLL, hessian = FALSE)
o

pred = G(o$par)
unique(pred)

outcome <- data.frame(data$agbd, pred)
outcome <- round(outcome, 3)
head(outcome)

median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

plot(median_values$pred, median_values$median_agbd, abline(0,1), xlim=c(0, 100))

############################
############################
# Assuming your dataframe is called 'my_data'
# and the column you want to use for splitting is 'column_name'

# Calculate the standard deviation of the column
std_dev <- sd(data$cwd)

# Calculate the thresholds for splitting into thirds
thresholds <- c(mean(data$prec) - 0.43*std_dev,
                 mean(data$prec) + 0.43*std_dev)

data_split_low <- data[data$prec < thresholds[1], ]
data_split_mid <- data[data$prec > thresholds[1] & data$prec < thresholds[2], ]
data_split_high <- data[data$prec > thresholds[2], ]

means <- aggregate(agbd ~ age, data_split_high, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'agbd', 'sd')

fit <- lm(sum_stats$median ~ sum_stats$age)
summary(fit)

ggplot(sum_stats, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))