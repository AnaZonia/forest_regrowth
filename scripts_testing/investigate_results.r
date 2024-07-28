################################################################################
#                                                                              #
#                 Comparing and investigating data                             #
#                                                                              #
################################################################################

library(corrplot)

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



tst <- read.csv("./data/fit_results.csv")

# Order the data frame by the column rsq
tst_ordered <- tst %>%
    arrange(desc(rsq))

# View the first few rows of the ordered data frame
head(tst_ordered)
