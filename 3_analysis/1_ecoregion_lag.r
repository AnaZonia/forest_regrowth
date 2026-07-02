
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      Compare R2 of different asymptotes and land use aggregations
#
#                     Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("all_points_amazon", biome = 1, n_samples = "all", asymptote = "nearest_mature", categorical = c("topography", "last_lu"))

# regions with enough representation of older ages to have a realistic estimate of the lag
regions <- c(476, 481, 507, 508, 518)

data <- subset(data, data$ecoreg %in% regions)

df_results <- data.frame()


for (ecoregion in unique(data$ecoreg)) {
    print(ecoregion)

    df_ecoreg <- subset(data, data$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 3)

    result <- c(ecoregion, mean(cv_results[[3]]), sd(cv_results[[3]]))
    df_results <- rbind(df_results, result)
    write.csv(df_results, file = "./0_results/lag_per_ecoregion.csv")
}

df_results

# avg_biomass <- data %>%
#     filter(ecoreg %in% regions) %>%
#     group_by(ecoreg, age) %>%
#     summarise(mean_biomass = mean(biomass, na.rm = TRUE), .groups = "drop") %>%
#     mutate(grp = paste("ecoreg", ecoreg))

# ggplot(avg_biomass, aes(x = age, y = mean_biomass, color = grp)) +
#     geom_line(linewidth = 1) +
#     labs(x = "Age", y = "Average biomass", color = "Region") +
#     theme_minimal()


df_results <- read.csv("./0_results/lag_per_ecoregion.csv")

head(df_results)


growth_rates <- data.frame()

for (ecoregion in unique(data$ecoreg)) {
    print(ecoregion)

    df_ecoreg <- subset(data, data$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 2)

    pars <- cv_results[[4]][1, ]

    k <- rep(pars[["k0"]], nrow(df_ecoreg))

    k <- (k + rowSums(sapply(names(pars), function(par) {
        pars[[par]] * df_ecoreg[[par]]
    }, simplify = TRUE))) * (df_ecoreg$age)

    growth_rates <- rbind(growth_rates, c(ecoregion, k))
    # write.csv(df_results, file = "./0_results/lag_per_ecoregion.csv")
}



tst <- sapply(names(pars), function(par) {
    pars[[par]] * df_ecoreg[[par]]
})

head(tst[[1]])





