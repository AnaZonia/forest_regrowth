library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_data_utils.r")
source("./2_R_scripts/2_model_utils.r")

set.seed(123)
re_base <- dnorm(1000)

# Import data
data <- import_data("./0_data/non_aggregated.csv", convert_to_dummy = TRUE)
data <- data %>% filter(biome == 1)
# unseen_data <- data %>% sample_n(2000)
# data <- data %>% sample_n(10000)
data_old <- subset(data, data$age == 34)
data_new <- subset(data, data$age == 1)
h1 <- hist(data_new$nearest_mature_biomass, plot = FALSE, breaks = seq(0, 500, length.out = 30))
h2 <- hist(data_old$agbd, plot = FALSE, breaks = seq(0, 500, length.out = 30))
h3 <- hist(data_new$agbd, plot = FALSE, breaks = seq(0, 500, length.out = 30))
mean_h1 <- mean(data$nearest_mature_biomass, na.rm = TRUE)
mean_h2 <- mean(data_old$agbd, na.rm = TRUE)
mean_h3 <- mean(data_new$agbd, na.rm = TRUE)

# Set up the plot area
plot(h1,
    col = rgb(0, 0, 1, 1 / 4), xlim = c(0, 500),
    main = "Forest Biomass - Amazon", xlab = "Value"
)

# Add the second histogram
plot(h2, col = rgb(1, 0, 0, 1 / 4), add = TRUE)

# Add the second histogram
plot(h3, col = rgb(0, 1, 0, 1 / 4), add = TRUE)

# Add a legend with the means included
legend("topright",
    legend = c(
        paste("Nearest Mature (Mean:", round(mean_h1, 2), ")"),
        paste("34 year old forest (Mean:", round(mean_h2, 2), ")"),
        paste("1 year old forest (Mean:", round(mean_h3, 2), ")")
    ),
    fill = c(rgb(0, 0, 1, 1 / 4), rgb(1, 0, 0, 1 / 4), rgb(0, 1, 0, 1 / 4))
)

graphics.off()


basic_pars <- c("k0", "m_base", "sd_base", "theta", "sd", "B0", "k")


curve((1 - exp(-20*x)), from = 0, to = 1, col = "red", lwd = 2, ylab = "Growth rate", xlab = "Time since disturbance")

# df <- subset(data, select = -c(agbd, biome))
# par_columns <- colnames(df)


# par_columns <- c("sur_cover")
# par_columns <- c("sur_cover", "cwd", "num_fires_before_regrowth", "mean_soil")
# pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
# pars[["B0"]] <- 0
# pars[["k"]] <- 0.5
# pars[["sd"]] <- 1

# conditions <- list(
#     function(pars) pars["B0"] < 0,
#     function(pars) pars["k"] < 0,
#     function(pars) pars["sd"] < 0
# )

# # Pass the growth_func to the likelihood function
# model <- optim(pars, function(pars) {
#     likelihood_B0_theta(pars, data, conditions, growth_curve_B0)
# })

# pred <- growth_curve_B0(model$par, data)
# calc_rsq(data, pred)


# # Run the optimization using the optim function
# model <- optim(par = pars, fn = likelihood_B0_theta, data = data, conditions = conditions)

# calc_rsq(data, growth_curve_B0_theta(model$par, data))


# par_columns <- c("sur_cover", "cwd")
# pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
# pars[["theta"]] <- 1
# pars[["k0"]] <- 0
# pars[["m_base"]] <- 0
# pars[["sd_base"]] <- 1
# pars[["sd"]] <- 1

# conditions <- list(
#     function(pars) pars["theta"] > 10,
#     function(pars) pars["theta"] < 0,
#     function(pars) pars["sd"] < 0,
#     function(pars) pars["sd_base"] < 0,
#     function(pars) pars["m_base"] < 0
# )

# # Run the optimization using the optim function
# model <- optim(pars, function(pars) {
#     likelihood_lag(pars, data, conditions)
# })
# pred = unname(growth_curve_lag(model$par, data, exp((re_base + model$par["m_base"]) * model$par["sd_base"])))
# calc_rsq(data, pred)

# Calculate histograms


# source("./2_R_scripts/2_model_utils.r")

# par_columns <- c("sur_cover")
# # par_columns <- c("sur_cover", "cwd", "num_fires_before_regrowth", "mean_soil")
# conditions <- list(
#     function(pars) pars["theta"] > 10,
#     function(pars) pars["theta"] < 0,
#     function(pars) pars["sd"] < 0
# )

# pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
# pars[["theta"]] <- 1
# pars[["k0"]] <- 0
# pars[["B0"]] <- mean(data$agbd)
# pars[["sd"]] <- 1

# # Pass the growth_func to the likelihood function
# model <- optim(pars, function(pars) {
#     likelihood_B0_theta(pars, data, conditions, growth_curve_B0_theta)
# })

# calc_rsq(data, growth_curve_B0_theta(model$par, data))
