####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################
library(ggplot2)
library(terra)
setwd("/home/aavila/forest_regrowth/")

rasters_files <- list.files(path = "./data/drive_export", pattern = '\\.tif$', full.names = TRUE)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}


########################
# SWITCHES
# data_prec_temp = FALSE
# scaled = FALSE

# tst <- read.csv('santoro_alldata.csv')

# minMax <- function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }

# if (data_prec_temp == TRUE){
#   data <- readRDS('santoro_ESA_alldata.rds')
#   data <- data[data$last_LU %in% c(15, 41, 48),]
#   data$last_LU <- factor(data$last_LU)
#   dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
#   names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
#   data_raw <- data[,-7]
#   if(scaled==TRUE){
#     data_scaled <- cbind(agbd=data_raw$agbd, scale(data_raw[,2:ncol(data_raw)]))
#     data <- cbind(data_scaled, dummy_LU)
#   }else{
#     data_maxmin <- cbind(agbd=data_raw$agbd, minMax(data_raw[,2:ncol(data_raw)]))
#     data <- cbind(data_maxmin, dummy_LU)
#   }
# }else{
#   data <- readRDS('total_aet_cwd.rds')
# }

######################################################################
#################        passing into the model     ##################
######################################################################

sds <- aggregate(sd ~ age, data, sd)
means <- aggregate(sd ~ age, data, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'agbd', 'sd')

fit <- lm(sum_stats$agbd ~ sum_stats$age)
summary(fit)

ggplot(sum_stats, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))

fit <- lm(data$agbd ~ data$age)
summary(fit)

ggplot(data, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))

# Predict values using the model
predicted_values <- predict(fit, as.data.frame(data$age))

mean((data$agbd - predicted_values)^2)
mean((data$agbd - pred)^2)

# --- NLS

pars_0 = c(A = 133, age = 0.02992, theta = 1.1162) # this is the value she found with nls function
# and with RSS = 2695.

G <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

# Indeed, starting out with her parameters as initial parameters we get better results:
pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407, sd = 0.5) # this is the value she found with nls function
Gpred <- G(pars, sum_stats)
Gpred

nls <- function(pars, data) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 10){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  if (pars['theta'] < 0){
    return(-Inf) 
  }
  result = sum((G(pars, data) - data$agbd)^2)
  ifelse(result == 0, -Inf, result)
}


o = optim(par = pars_0, fn = nls, data = sum_stats)
o

pred <- G(o$par, sum_stats)
 
outcome <- data.frame(sum_stats$agbd, pred)
head(outcome)
outcome <- round(outcome, 3)
head(outcome)

plot(outcome$pred, outcome$data.agbd, abline(0,1)) # , xlim=c(0, 100))
mean((data$agbd - pred)^2)


# -- NLL

G <- function(pars, data) {
  pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

pars = c(A = 87.07455636, age = 0.07435007, theta = 1.69029407, sd = 0.5) # this is the value she found with nls function
Gpred <- G(pars, data)
Gpred

NLL = function(pars) {
  if (pars['age'] < 0){
    return(-Inf) 
  }
  if (pars['age'] > 10){
    return(-Inf) 
  }
  if (pars['theta'] > 10){
    return(-Inf) 
  }
  if (pars['theta'] < 0){
    return(-Inf) 
  }
  if (pars['sd'] < 0){
    return(-Inf) 
  }
  if (pars['sd'] > 10){
    return(-Inf) 
  }
  # Negative log-likelihood 
  result = -sum(dnorm(x = data$agbd - G(pars, data), mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
  ifelse(result == 0, -Inf, result)
}

o = optim(par = pars, fn = NLL, hessian = FALSE)
o

pred <- G(o$par)
 
outcome <- data.frame(data$agbd, pred)
outcome <- round(outcome, 3)
head(outcome)
plot(outcome$pred, outcome$data.agbd, abline(0,1)) # , xlim=c(0, 100))
mean((data$agbd - pred)^2)

median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

plot(median_values$pred, median_values$median_agbd, abline(0,1), xlim=c(0, 100))

############################
############################

breaks <- quantile(data$cwd, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
categories <- cut(data$cwd, breaks, labels = c("low", "mid", "high"))
df_list <- split(data, categories)

sum_stats <- aggregate(agbd ~ age, df_list[[3]], median)
colnames(sum_stats) <- c('age', 'agbd')
fit <- lm(sum_stats$agbd ~ sum_stats$age)
summary(fit)

predicted_values <- predict(fit, as.data.frame(sum_stats$age))
mean((sum_stats$agbd - predicted_values)^2)


ggplot(sum_stats, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))+
  labs(title = "Median, low cwd")+
  theme(plot.title = element_text(size = 40)) 

##############################################################

G <- function(pars, data) {
  pars['B0'] + pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

# now fit with land use and fire and see what happens.

setwd("/home/aavila/forest_regrowth")
data <- readRDS('santoro_cwd_fire_LU.rds')

pasture <- data[data$last_LU == 15,]
other_annual <- data[data$last_LU == 41,]
other_perennial <- data[data$last_LU == 48,]
dfs <- list(pasture, other_annual, other_perennial)
dfs <- lapply(dfs, aggregate, agbd ~ age, median)

colors <- c("red", "blue", "green")
plot(NULL, xlim = range(data$age), ylim = range(0:150), xlab = 'years', ylab = 'agbd',
main = 'Land Use Type')

for (i in 1:length(dfs)) {
  points(dfs[[i]]$age, dfs[[i]]$agbd, col = colors[i], pch = i) # Change 'points' to 'lines' if needed
}

legend("topright", legend = c("pasture", "other_annual", "other_perennial"), col = colors, pch = 1:3)

pars <- c(B0 = 40, A = 100, age = 0.05, theta = 1.5)
o = optim(par = pars, fn = nls, data = dfs[[1]], method='BFGS')
o
 
lines(G(o$par, data = dfs[[1]]), col = "red", lwd = 2, lty = 2)

o = optim(par = pars, fn = nls, data = dfs[[2]], method='BFGS')
o
 
lines(G(o$par, data = dfs[[2]]), col = "blue", lwd = 2, lty = 2)

o = optim(par = pars, fn = nls, data = dfs[[3]], method='BFGS')
o
 
lines(G(o$par, data = dfs[[3]]), col = "green", lwd = 2, lty = 2)

#########################################
# looking at fire

median_fire_total <- aggregate(agbd ~ fire_total, data, median)
median_last_fire <- aggregate(agbd ~ last_fire, data, median)

plot(NULL, xlim = range(data$age), ylim = range(0:150), xlab = 'years', ylab = 'agbd',
main = 'Fire')

points(median_fire_total$fire_total, median_fire_total$agbd, col = 'red', pch = 1)
points(median_last_fire$last_fire, median_last_fire$agbd, col = 'blue', pch = 2)

legend("topright", legend = c("Total fires", "Time since fire"), col = c('red', 'blue'), pch = 1:2)

G <- function(pars, data) {
  pars['B0'] + pars['last_fire']*data$last_fire
}

pars <- c(B0 = 40, last_fire = 0.05)

o = optim(par = pars, fn = nls, data = median_last_fire, method='BFGS')
o
 
lines(G(o$par, data = median_last_fire), col = "blue", lwd = 2, lty = 2)

G <- function(pars, data) {
  pars['a']*data$fire_total^2+pars['b']*data$fire_total+pars['c']
}

pars <- c(a = 2, b = 0.05, c = 1)

o = optim(par = pars, fn = nls, data = median_fire_total, method='BFGS')
o
 
lines(G(o$par, data = median_fire_total), col = "red", lwd = 2, lty = 2)

G <- function(pars, data) {
  pars['a'] + pars['b']*exp(-data$fire_total*pars['c'])
}
pars <- c(a = 2, b = 0.05, c = 0.05)

o = optim(par = pars, fn = nls, data = median_fire_total, method='BFGS')
o

G(o$par, median_fire_total)

lines(G(o$par, data = median_fire_total), col = "red", lwd = 2)

G <- function(pars, data) {
  pars['B0'] + pars['A'] * (1 - exp(-pars['age']*data$age))^pars['theta']
}

median_age_fire_total <- aggregate(fire_total ~ age, data, mean)
median_age_last_fire <- aggregate(last_fire ~ age, data, mean)

plot(median_age_last_fire$age, median_age_last_fire$last_fire)
plot(median_age_fire_total$age, median_age_fire_total$fire_total)

#######################################
# looking at total values
 
median_agbd <- aggregate(agbd ~ age, data, mean)
sd_agbd <- aggregate(agbd ~ age, data, sd) 
plot(median_agbd$age, median_agbd$agbd)
o = optim(par = pars, fn = nls, data = median_agbd, method='BFGS')
o
lines(G(o$par, data = median_agbd), col = "red", lwd = 2)
o = optim(par = o$par, fn = NLL, data = median_agbd, method='BFGS')
o
G(o$par)

merged_data <- merge(data, median_agbd, by = "age", all.x = TRUE)
merged_data <- merge(merged_data, sd_agbd, by = "age", all.x = TRUE)
colnames(merged_data) <- c('age', 'agbd', 'cwd', 'last_LU', 'fire_total', 'last_fire', 'mean', 'sd')
new_agbd <- merged_data[abs(merged_data$agbd - merged_data$mean) < 0.25*merged_data$sd, ]
plot(new_agbd$age, new_agbd$agbd)

pars = c(A = 100, age = 1, theta = 1, sd = 0.05)
o = optim(par = pars, fn = NLL, data = new_agbd)
o

G(o$par, data = new_agbd)

lines(G(o$par, data = new_agbd), col = "blue", lwd = 2)


#########################################
# looking at CWD

median_agbd <- aggregate(cwd ~ agbd, data, mean)
median_age <- aggregate(cwd ~ age, data, mean)
#plot(median_age$age, median_age$cwd)

color_gradient <- colorRampPalette(c("blue", "red"))

# Determine the colors based on age values
age_range <- range(median_age$age)
color_vector <- color_gradient(5)  # Generate 5 colors along the gradient
color_vector <- color_vector[as.numeric(cut(median_age$age, breaks = 5))]
# Create the scatter plot with the color gradient
plot(median_agbd$cwd, median_agbd$agbd, col = color_vector)

sd_agbd <- aggregate(cwd ~ agbd, data, sd) 

merged_data <- merge(data, median_agbd, by = "agbd", all.x = TRUE)
merged_data <- merge(merged_data, sd_agbd, by = "agbd", all.x = TRUE)
colnames(merged_data) <- c('age', 'agbd', 'cwd', 'last_LU', 'fire_total', 'last_fire', 'mean', 'sd')
new_agbd <- merged_data[abs(merged_data$cwd - merged_data$mean) < 0.25*merged_data$sd, ]
plot(new_agbd$cwd, new_agbd$agbd)
