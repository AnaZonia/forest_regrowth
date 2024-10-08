
library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(randomForest)


# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_import_data.r")
source("./2_R_scripts/2_functions.r")

set.seed(1)
ncores <- 20
registerDoParallel(cores = ncores)


# Import data

df <- read_csv("./0_data/eu.csv", show_col_types = FALSE)

head(df)


tst <- import_data("./0_data/eu.csv", convert_to_dummy = TRUE)
