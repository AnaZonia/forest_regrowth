
library(tidyverse)

setwd("/home/anazonia/Documents/forest_regrowth")


csv_files <- list.files(paste0("./0_data/UTM_10k"), pattern = "\\.csv$", full.names = TRUE)
df <- csv_files %>%
    map(read_csv) %>%
    bind_rows()

mean(df$biomass, na.rm = TRUE)


csv_files <- list.files(paste0("./0_data/grid_1k_amazon_secondary"), pattern = "\\.csv$", full.names = TRUE)
df <- csv_files %>%
    map(read_csv) %>%
    bind_rows()
