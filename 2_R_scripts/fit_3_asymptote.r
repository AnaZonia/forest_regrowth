# 1. Just age and biomass:
# - Asymptote for the climatic region
# - Asymptote for the ecoregion
# - Asymptote of the nearest neighbor

library(tidyverse)
library(ggplot2)

data <- read_csv("./0_data/unified_fc.csv") %>% na.omit()

calc_r2 <- function(data, pred) {
    obs_pred <- lm(data$biomass ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$biomass - mean(data$biomass))^2)
    r2 <- 1 - (sum_res_squared / total_sum_squares)
    return(r2)
}

growth_curve <- function(pars, data, biomass_col) {
    k <- pars[["k0"]] * data[["age"]]

    # Constrain k
    k[k < 1e-10] <- 1e-10
    k[k > 7] <- 7

    return(pars[["B0"]] + (data[[biomass_col]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
}

likelihood <- function(pars, data, conditions, biomass_col) {
    result <- mean((growth_curve(pars, data, biomass_col) - data$biomass)^2)

    # Check parameter constraints
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

run_optim <- function(train_data, pars, conditions, biomass_col) {

    model <- optim(pars, likelihood, data = train_data, conditions = conditions, biomass_col = biomass_col)
    return(model)
}

pars <- c()
pars["theta"] <- 1
pars["k0"] <- 1
pars["B0"] <- 100

conditions <- list('pars["k0"] < 0', 'pars["B0"] < 0', ' pars["theta"] > 10', 'pars["theta"] < 0')

results <- tibble(quarter = integer(), r2_quarter = numeric(), r2_ecoreg = numeric())

for (q in c(1, 2, 3, 4, NA)) {
    df <- data %>%
        filter(biome == 1) %>%
        {
            if (!is.na(q)) filter(., quarter == q) else .
        }

    df$biomass <- df$ESA_CCI_2020
    # Initialize variables to store R2 values for each biomass column
    r2_quarter <- NA
    r2_ecoreg <- NA
    r2_first <- NA

    # Loop over the biomass columns
    for (biomass_col in c("first", "quarter_biomass", "ecoreg_biomass")) {
        model <- run_optim(df, pars, conditions, biomass_col)
        pred <- growth_curve(model$par, df, biomass_col)
        r2 <- calc_r2(df, pred)

        # Store R2 values in the appropriate variable
        if (biomass_col == "quarter_biomass") {
            r2_quarter <- r2
        } else if (biomass_col == "ecoreg_biomass") {
            r2_ecoreg <- r2
        } else if (biomass_col == "first") {
            r2_first <- r2
        }
    }

    # Bind the results to the tibble
    results <- bind_rows(results, tibble(quarter = q, r2_quarter = r2_quarter, r2_ecoreg = r2_ecoreg, r2_first = r2_first))

    if (is.na(q)) {
        df_last <- df %>% mutate(pred = pred_ecoreg)
    }
}

print(results)


# Calculate median of pred by age group and define asymptote
mean_biomass_per_age <- df_last %>%
    group_by(age) %>%
    summarize(
        mean_biomass = mean(biomass, na.rm = TRUE),
        mean_pred = mean(pred, na.rm = TRUE)
    )


ggplot(mean_biomass_per_age, aes(x = age)) +
    geom_point(aes(y = mean_biomass), color = "blue", size = 2) +
    geom_line(aes(y = mean_pred), color = "red", linewidth = 1) +
    labs(
        title = "Mean Biomass per Age with Fitted Curve",
        x = "Age",
        y = "Biomass",
        caption = "Blue points: Observed mean biomass | Red line: Fitted growth curve"
    ) +
    theme_minimal()
