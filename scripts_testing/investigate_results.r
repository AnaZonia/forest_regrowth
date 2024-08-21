################################################################################
#                                                                              #
#                 Comparing and investigating data                             #
#                                                                              #
################################################################################

library(corrplot)
library(tidyverse)

source("fit_1_import_data.r")

# Define land-use history intervals to import four dataframes
intervals <- list(
    "5yr",
    "10yr",
    "15yr",
    "all"
)

datafiles <- lapply(intervals, function(file) {
    paste0("./data/amaz_", file, ".csv")
})

dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

data <- dataframes[[1]]

tst <- read.csv("./data/dist_mature_1000m_countrywide.csv")


create_correlation_plot <- function(df, interval_name) {
    # Remove non-numeric columns
    numeric_df <- df %>% select_if(is.numeric)

    # Calculate correlation matrix
    cor_matrix <- cor(numeric_df, use = "pairwise.complete.obs")

    # Create correlation plot
    png(paste0("correlation_plot_", interval_name, ".png"), width = 1200, height = 1000, res = 100)
    corrplot(cor_matrix,
        method = "color",
        type = "upper",
        order = "hclust",
        tl.col = "black",
        tl.srt = 45,
        addCoef.col = "black",
        number.cex = 0.7,
        tl.cex = 0.7,
        title = paste("Correlation Matrix -", interval_name),
        mar = c(0, 0, 1, 0)
    )
    dev.off()

    print(paste("Correlation plot for", interval_name, "saved."))
}


results_optim <- read.csv("./data/amaz_results_all.csv")
results_lm <- read.csv("./data/countrywide_results_lm.csv")
results_rf <- read.csv("./data/countrywide_results_rf.csv")

head(results_optim)

results_all <- bind_rows(results_optim, results_lm, results_rf) %>%
    arrange(data_pars, basic_pars, data_name)
write.csv(results_all, "./data/results_all.csv", row.names = FALSE)

results_all <- read.csv("./data/results_all.csv")

top_5_per_data_pars <- results_all %>%
    group_by(data_pars) %>%
    top_n(5, rsq) %>%
    arrange(data_pars, desc(rsq)) %>%
    ungroup()

print(top_5_per_data_pars)
write.csv(top_5_per_data_pars, "./data/amaz_top_5_per_data_pars.csv", row.names = FALSE)

# Order the data frame by the column rsq
tst_ordered <- results_all %>%
    filter(model_type == "lm") %>%
    arrange(desc(rsq))

tst_ordered[c(1:5), c(1:5)]

mean_rsq_per_data_name <- results_lm %>%
    # filter(model_type == "optim")%>%
    group_by(data_pars) %>%
    summarise(mean_rsq = mean(rsq, na.rm = TRUE)) %>%
        arrange(desc(mean_rsq))

print(mean_rsq_per_data_name)

anova_result <- aov(rsq ~ data_name, data = results_all)

print(summary(anova_result))


tst <- read.csv("./data/all_LULC_countrywide.csv")

# Assuming results_all is your dataframe
mean_rsq_per_data_name <- tst %>%
    filter(biome == 1)%>%
    group_by(data_pars) %>%
    summarise(mean_rsq = mean(rsq, na.rm = TRUE)) %>%
    arrange(desc(mean_rsq))


data <- read.csv("./data/amaz_results_optim.csv")

best_rsq_df <- data %>%
    group_by(data_name, data_pars) %>%
    slice_max(rsq, with_ties = FALSE) %>%
    arrange(data_pars, data_name)
best_rsq_df[, c(1:5)]


# Function to count occurrences of a substring in a file
count_occurrences <- function(file_path, substring) {
    count <- 0
    file_content <- readLines(file_path)

    for (line in file_content) {
        if (grepl(substring, line)) {
            count <- count + 1
        }
    }

    return(count)
}

total_string_count <- count_occurrences("nohup_amaz.out", "Time so far")
print(paste("Total occurrences of string:", total_string_count))




amazon <- read.csv("./data/amaz_5yr.csv")

amazon_subset <- read.csv("./data/countrywide_5yr.csv")
amazon_subset <- subset(amazon_subset, biome == 1)

range(amazon$cwd)
range(amazon_subset$cwd)

summary(lm(agbd ~ cwd, amazon))
summary(lm(agbd ~ cwd, amazon_subset))
summary(lm(agbd ~ mature_biomass, amazon))
summary(lm(agbd ~ mature_biomass, amazon_subset))

# Example character vector
# Example character vector
vec <- c(
    "theta:1.31979333088066",
    "fallow:5.30356802403215",
    "lulc_sum_41:0.129769992345671",
    "num_fires_after_regrowth:-0.514237069703973",
    "num_fires_before_regrowth:-0.0422835348018823",
    "protec:0.0624703092189218",
    "sur_cover:0.723205699759347",
    "ts_fire_after_regrowth:-0.288191882147403",
    "ts_fire_before_regrowth:-0.0358069389792399",
    "ecoreg_446:10.149607211274",
    "ecoreg_464:0.326839933603352",
    "ecoreg_465:0.0166891133152462",
    "ecoreg_466:0.375550319421151",
    "ecoreg_467:-0.075066596486651",
    "ecoreg_469:4.20822231013806",
    "ecoreg_473:0.54306256718872",
    "ecoreg_474:0.433638435318312",
    "ecoreg_476:0.380942011954866",
    "ecoreg_480:-0.0299180281881416",
    "ecoreg_481:0.193396063436363",
    "ecoreg_482:0.435423794889215",
    "ecoreg_484:0.66260358371134",
    "ecoreg_485:0.378149989980357",
    "ecoreg_490:10.3259203819715",
    "ecoreg_496:0.340517292109795",
    "ecoreg_497:0.424245959872319",
    "ecoreg_498:0.160239549461896",
    "ecoreg_503:0.21325730789524",
    "ecoreg_505:1.17566454986428",
    "ecoreg_507:0.671915895164663",
    "ecoreg_508:0.408966481117558",
    "ecoreg_511:0.0499144529129286",
    "ecoreg_518:0.244738032607572",
    "ecoreg_529:0.279182351982756",
    "ecoreg_540:0.5658506478689",
    "ecoreg_567:0.23536885714865",
    "ecoreg_570:0.0431826319683144",
    "ecoreg_584:0.335136359975363",
    "ecoreg_611:0.0133026489911424",
    "B0:87.5818691937891",
    "age:0.72828971492447"
)

# Split each element into name and value
split_vec <- strsplit(vec, ":")

# Convert to a named list
pars <- setNames(
    lapply(split_vec, function(x) as.numeric(x[2])),
    sapply(split_vec, function(x) x[1])
)

# Print the result
print(result_list)


A <- runif(35, 120, 199)
theta <- 7
B0 <- 40
k <- 0.2

curve(B0 + (A - B0) * (1 - exp(-k * x))^theta, from = 1, to = 35, n = 35)

curve((B0 * A * exp(k * x)) / ((A - B0) + B0 * exp(k * x)), from = 1, to = 35, n = 35)
