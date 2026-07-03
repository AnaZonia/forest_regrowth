
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



path = "all_points_amazon"
biome = 1
n_samples = "all"
asymptote = "nearest_mature"
categorical = c("topography", "last_lu")




csv_files <- list.files(paste0("./0_data/", path), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

# remove columns with all NA values
df <- df[, colSums(is.na(df)) < nrow(df)]


# remove any rows with NA values
df <- df %>% filter(rowSums(is.na(.)) == 0)

# Convert categorical to factors
df <- df %>%
    mutate(across(any_of(categorical), as.factor)) %>%
    filter(biome == biome)

# remove columns with less than 50 unique values
df <- df %>%
    group_by(across(any_of(categorical))) %>%
    ungroup() %>%
    mutate(across(any_of(categorical), droplevels))

df <- dummy_cols(df,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

asymptotes <- c("nearest_mature", "ecoreg_biomass", "quarter_biomass")

if (asymptote == "full_amazon") {
    df$asymptote <- mean(df$nearest_mature, na.rm = TRUE)
    df <- df %>% select(-any_of(c(asymptotes, "quarter", "biome")))
} else {
    # remove the columns in asymptotes that are not the designated asymptote
    # rename to asymptote the column named the same as the value of asymptote
    df <- df %>% rename(asymptote = !!sym(asymptote))
    df <- df %>% select(-any_of(c(
        asymptotes[asymptotes != asymptote],
        "quarter", "biome"
    )))
}

# remove columns with less than 50 non-zero values
df <- df %>% select(where(~ sum(. != 0) >= 50))





sampled_df <- df[sample(nrow(df), min(30000, nrow(df)), replace = FALSE), ]


grid_data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 30000, asymptote = "nearest_mature", categorical = c("topography", "last_lu", "ecoreg"))

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(grid_data))[["all"]]
data_pars
source("2_modelling/2_modelling.r")

cv_results <- cross_validate(grid_data, basic_pars, data_pars, conditions, 5)

saveRDS(cv_results, "./0_results/tst/ecoreg_dummy_lag.rds")


topog_pars <- data_pars[[]]




grid_data <- data.frame(age = grid_data$age, biomass = grid_data$biomass, type = "grid")
df_comparison <- data.frame(age = df$age, biomass = df$biomass, type = "polygons")

plot_data <- rbind(grid_data, df_comparison)

avg_biomass <- plot_data %>%
    group_by(age, type) %>%
    summarise(mean_biomass = mean(biomass, na.rm = TRUE), .groups = "drop") %>%
    mutate(grp = paste("type", type))

intercepts <- avg_biomass %>%
    group_by(type) %>%
    slice_min(order_by = age, n = 1, with_ties = FALSE)

ggplot(avg_biomass, aes(x = age, y = mean_biomass, color = type)) +
    geom_line(linewidth = 1.2) +
    geom_point(data = intercepts, size = 3) +
    geom_text(
        data = intercepts,
        aes(label = round(mean_biomass, 2)),
        vjust = -1,
        show.legend = FALSE,
        size = 5
    ) +
    labs(
        x = "Age",
        y = "Average biomass",
        color = "Type",
        title = "Average biomass by age"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold")
    )


regions <- c(476, 481, 507, 508, 518)


df_results <- data.frame()

for (ecoregion in regions) {
    print(ecoregion)

    # df_ecoreg <- subset(sampled_df, sampled_df$ecoreg == ecoregion)
    df_ecoreg <- df %>%
        filter(ecoreg == ecoregion) %>%
        slice_sample(by = age, n = 1000)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 3)

    result <- c(ecoregion, mean(cv_results[[3]]), sd(cv_results[[3]]))
    df_results <- rbind(df_results, result)
    write.csv(df_results, file = "./0_results/lag_per_ecoregion.csv")
}

df_results



# summary:
# results with random sampling 10k points (ecoregion)
# results with sampling 10k points being one per polygon of the same age (ecoregion)
# ecoregion as a dummy variable





