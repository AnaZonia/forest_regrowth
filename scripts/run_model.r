####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
####################################################################
# Ana Avila - Dec 2023
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

library(ggplot2)
library(terra)
library(tidyverse) # for stringr operations
setwd("/home/aavila/forest_regrowth/")
source("./scripts/regrowth_functions.r")

# ----------------------------------
# Switches
# ----------------------------------

make_amazon_df = F

# ----------------------------------
# Importing data from .tif files and converting into dataframe
# ----------------------------------

# 10 k per region would be okay
# would patterns emerge with more data?
# incorporating temporal dynamics

if (make_amazon_df == T){
  files <- list.files(path = "./data/drive_export", pattern = '\\.tif$', full.names = TRUE)
  rasters <- lapply(files, rast)
  groups <- split(rasters, gsub("-[0-9]{10}-[0-9]{10}", "", files))

  count <- 0
  merge_rasters <- function(raster_list) {
    if (length(raster_list) > 1) {
      do.call(merge, raster_list)
      count <<- count + 1
      print(count)
    } else {
      raster_list[[1]]
      count <<- count + 1
      print(count)
    }
  }

  # Apply the merging function to each group
  merged_rasters <- lapply(groups, merge_rasters)

  count <- 0
  rast_to_df <- function(raster){
    print('----------------')
    count <<- count + 1
    print(count)
    coords <- crds(raster, df=FALSE, na.rm=TRUE)
    print('done with coords')
    values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
    print('done with values')
    central_df <- values_stack[complete.cases(values_stack), ]
    return(central_df)
  }

  dataframes_per_ecoregion <- lapply(merged_rasters, rast_to_df)
  numbers <- unique(gsub("^img_export_([0-9]{3}).*", "\\1", basename(files)))
  dataframes_per_ecoregion <- map2(dataframes_per_ecoregion, numbers, ~mutate(.x, ecoregion = .y))

  amazon_df <- do.call(rbind, dataframes_per_ecoregion)
  rownames(amazon_df) <- NULL
  write.csv(amazon_df, file = "amazon_df.csv")
}

# ----------------------------------
# Passing into the model
# ----------------------------------

# inputting parameters from nls to NLL
#       o <- optim(c(o_nls$par, sd = 0.5), NLL, data=data_example)
# repuffing the model with the new parameters
# fixing SD to 0.5
# 

amazon_df <- read.csv('./data/amazon_df_sample_10million.csv')
dummy_factor <- as.factor(amazon_df$ecoregion)
dummy_matrix <- model.matrix(~dummy_factor - 1)  # -1 removes intercept
amazon_df_dummy <- cbind(amazon_df, as.data.frame(dummy_matrix))

sds <- aggregate(agbd ~ age, amazon_df, sd)
means <- aggregate(agbd ~ age, amazon_df, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'agbd', 'sd')

ggplot(sum_stats, aes(x = age, y = agbd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme(text = element_text(size = 20))

fit <- lm(agbd ~ age, amazon_df)
summary(fit) # R-squared:  0.09077

model2 <- lm(agbd ~ age + cwd + ecoregion,
             data = amazon_df)
summary(model2) # R-squared: 0.1226

# Predict values using the linear model
predicted_values <- predict(fit, as.data.frame(amazon_df$age))

mean((amazon_df$agbd - predicted_values)^2) #3986.709

# ------------ NLS ------------ #

data <- amazon_df_dummy

pars = c(B0 = 10, A = 133, age = 0.02992, theta = 1.1162,
cwd = 0.001, d485 = 0.001, d508 = 0.01, d529 = 0.01)

Gpred <- G_6par(pars, data)
head(Gpred)

o = optim(pars, fn = nls, data = data, G = G_6par)
o

pred <- G_6par(o$par, data)

outcome <- data.frame(data$agbd, pred)
head(outcome)

median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

plot(median_values$pred, median_values$median_agbd, abline(0,1))

# with the medians
mean((data$agbd - pred)^2) #3961.879 with G_3par
# 3945.359 with G_4par
# 3963.008 with G_6par (including CWD)
# 3868.773 with dummy variables

# use them as variance explained
# not directly related with probability

# divide by the sum of squares of the observed values 
# complement of 1-MSE/Sum of squares

# use LM on observed vs predicted
# find if the intercept and that the LM will find is zero and the slope is 1.

# if asymptote is a constant, it should not make a difference, but if it a function
# of surrounding forest height, then it should make a difference. the asymptote should

# ------------ with dummy variable ------------ #

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
