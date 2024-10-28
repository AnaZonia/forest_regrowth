library(tidyverse)
library(fastDummies)

data_good <- read_csv("./0_data/non_aggregated.csv", show_col_types = FALSE)

data_bad <- read_csv("./0_data/non_aggregated_all.csv", show_col_types = FALSE)


nrow(data_good)
# colnames(data_bad)
nrow(data_bad)

# data_good <- subset(data_good, biome == 1)
# data_bad <- subset(data_bad, biome == 1)

hist(data_good$age)
hist(data_bad$age)
table(data_good$lulc_sum_48)



colnames(data_good)


data_good <- read_csv("./0_data/non_aggregated_5yr.csv", show_col_types = FALSE)
data_good <- subset(data_good, biome == 4)

# List of categorical variable names
categorical_vars <- c("ecoreg", "last_LU", "topography")

# Convert specified columns to factors
data_good <- data_good %>%
    mutate_at(vars(categorical_vars), as.factor) %>%
    select(-biome)

# Construct the formula
all_predictors <- paste(setdiff(names(data_good), c("biome", "agbd")), collapse = " + ")
non_lu_predictors <- paste(setdiff(names(data_good)[1:20], c("biome", "agbd")), collapse = " + ")
lu_predictors <- paste(setdiff(names(data_good)[21:33], c("biome", "agbd")), collapse = " + ")

# Fit the linear model and summarize it
model <- lm(as.formula(paste("agbd ~", lu_predictors)), data = data_good)
summary(model)
model <- lm(as.formula(paste("agbd ~", non_lu_predictors)), data = data_good)
summary(model)
model <- lm(as.formula(paste("agbd ~", all_predictors)), data = data_good)
summary(model)
model <- lm(as.formula(paste("agbd ~ age")), data = data_good)
summary(model)
