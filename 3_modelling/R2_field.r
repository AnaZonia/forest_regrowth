# check field data

library(dplyr)
library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_cross_validate.r")
source("3_modelling/2_feature_selection_ga.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# Load data
biome <- 1
n_samples = 10000

# identify which ones are in the same site
unified_field <- read.csv("./0_data/groa_field/unified_field.csv") %>%
    # remove columns system.index and .geo
    select(-c(system.index, .geo)) %>%
    mutate(date = if_else(date < 0, NA_real_, date))

plot(unified_field$field_age, unified_field$field_biom, xlab = "Field Age", ylab = "Field Biomass", main = "Field Age vs Field Biomass")


# try it out with only 5 years old onwards to see if theta remains 1

data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)
# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data

pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = c("num_fires", "dist"), norm_data)

final_model <- run_optim(norm_data, pars_init, conditions)

pred <- growth_curve(final_model$par, data = norm_data)

r2 <- calc_r2(norm_data, pred)
r2

tst <- as.data.frame(cbind(norm_data$age, pred))
colnames(tst) <- c("age", "pred")
head(tst)
plot(tst$age, tst$pred, xlab = "Field Age", ylab = "Predicted AGB", main = "Field Age vs Predicted AGB")

plot(norm_data$age, norm_data$biomass, xlab = "Field Age", ylab = "Predicted AGB", main = "Field Age vs Predicted AGB")



plot(pred, norm_data$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)

plot(pred, norm_data$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")


field <- unified_field %>%
    # rename field_age to age
    rename(age = field_age) %>%
    # rename first to nearest_biomass
    rename(nearest_biomass = first) %>%
    # rename b1 to biomass
    rename(biomass = field_biom) %>%
    # keep only columns dist, num_fires, sur_cover, nearest_biomass, biomass and age
    select(dist, num_fires, nearest_biomass, biomass, age)

norm_field <- normalize_independently(field)$train_data
nrow(norm_field) # 582
# remove rows with NA in any column
norm_field <- norm_field[complete.cases(norm_field), ]
nrow(norm_field) # 402

pars_init <- find_combination_pars(
    basic_pars = c(basic_pars_options[["intercept"]], "theta"),
    data_pars = c("num_fires", "dist"),
    norm_field
)
final_model <- run_optim(norm_field, pars_init, conditions)
final_model$par


plot(norm_field$age, norm_field$biomass, xlab = "Field Age", ylab = "Field Biomass", main = "Field Age vs Field Biomass")
pred_agb <- growth_curve(final_model$par, data = norm_field)
points(norm_field$age, pred_agb, col = "red", pch = 19)

# get R2
r2 <- calc_r2(norm_field, pred_agb)
r2

plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)

# color the points by their age
plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB", col = tst$age)
# make the color a gradient from blue to red
library(ggplot2)
ggplot(tst, aes(x = pred_agb, y = biomass, color = age)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "Predicted AGB", y = "Observed AGB", title = "Predicted vs Observed AGB") +
    theme_minimal()

# tell me if there is a relationship between age and distance between observed from predicted AGB
# calculate the difference between observed and predicted AGB
tst$diff <- tst$biomass - pred_agb
# plot the difference against age
plot(tst$diff, tst$age, xlab = "Difference between observed and predicted AGB", ylab = "Age", main = "Difference vs Age")
summary(lm(tst$diff ~ tst$age))

