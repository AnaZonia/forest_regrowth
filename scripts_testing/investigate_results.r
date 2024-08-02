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
    "5y",
    "10y",
    "15y",
    "all"
)

datafiles <- lapply(intervals, function(file) {
    paste0("./data/", file, "_LULC_mat_dist.csv")
})

dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

head(dataframes[[1]])

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


results_optim <- read.csv("./data/countrywide_results_optim.csv")
results_lm <- read.csv("./data/countrywide_results_lm.csv")
results_rf <- read.csv("./data/countrywide_results_rf.csv")


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


tst <- read.csv("./data/amaz_15yr.csv")
colnames(tst)




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

total_string_count <- count_occurrences("nohup.out", "Time so far")

print(paste("Total occurrences of string:", total_string_count))

# should take 54 min to go through all optim iterations with 
224/80*19.4691248297691
