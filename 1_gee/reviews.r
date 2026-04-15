

library(tidyverse)

csv_files <- list.files(paste0("./0_data/grid_10k_amazon_secondary_allpixels"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()



head(df)

unique(df$ecoreg)

table(df$edge)

# df <- df %>%
#     filter(ecoreg == 518)

# nrow(df)

# edge_pixels <- df %>%
#     filter(edge == 1)

# non_edge_pixels <- df %>%
#     filter(edge == 0)

# mean(edge_pixels$biomass, na.rm = TRUE)
# mean(non_edge_pixels$biomass, na.rm = TRUE)
