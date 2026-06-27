# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#    Predictions for Regrowth by 2050 (pasture and secondary)
#
#                 Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

library(tidyverse)
library(terra)
library(scales) # for label formatting
library(foreach)
library(doParallel)

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- Future predictions - Bezerra map ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make future projections with the general model from error propagation
# get the projection for each pixel selected from the 10km grid cell
# get the area expected for that grid cell, then from that get the total

# get average carbon per area for those cells

# multiply that by the total area regrown to get the estimated

# ----------------------------------------------------

library(terra)

# veg <- forest vegetation
# gveg <- grassland vegetation
# mosc <- mosaic vegetation
# fores <- forestry

for (scenario in c("SSP1_RCP19", "SSP2_RCP45", "SSP3_RCP70")) {
    r <- rast(paste0("./0_data/LUCCMEBR_", scenario, "_land_cover_type_100km2_2015_2050.nc"))

    forest_2050 <- r$veg_8
    forest_2015 <- r$veg_1

    writeRaster(forest_2050, paste0("./0_data/forest_2050_", scenario, ".tif"), overwrite = TRUE)
    writeRaster(forest_2015, paste0("./0_data/forest_2015_", scenario, ".tif"), overwrite = TRUE)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------- Calculate future carbon sequestration ---------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

error_prop_results <- readRDS("./0_results/error_prop.rds")
pars <- colMeans(error_prop_results[[2]])

future <- read.csv("./0_data/future_scenarios_area.csv") %>%
        rename(asymptote = nearest_mature)

future <- subset(future, !is.na(cec)) # remove rows where cec is NA (deal with this later)
norm_future <- normalize_independently(future)$train_data

future <- future %>%
        select(
            growth_SSP1_RCP19, growth_SSP2_RCP45, growth_SSP3_RCP70
        ) %>%
        cbind(norm_future[, -c(1:3)]) %>%
        mutate(
            age = 35,
            across(c(growth_SSP1_RCP19, growth_SSP2_RCP45, growth_SSP3_RCP70), ~ .x * 100) # km2 → ha
        )

future <- dummy_cols(future,
    select_columns = "topography",
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)


# ── Compute SSP totals ────────────────────────────────────────────────────────

compute_ssp_totals <- function(future, growth_col, pars) {
    df <- subset(future, !is.na(future[[growth_col]]))
    df$pred <- growth_curve(pars, df, pars["lag"])
    df <- subset(df, !is.na(pred))

    # assume 25.8% of total biomass is belowground
    # convert biomass in Mg/ha to MgC/ha (assuming 50% C content)
    # total estimate in Teragram of Carbon (TgC / ha)
    # Biomass conversions: +belowground (25.8%) → MgC/ha (50% C) → TgC
    total_biomass_tgc <- sum(df$pred) * 1.258 * 0.5 / 1e6
    total_area <- sum(df[[growth_col]]) # ha

    list(
        total_TgC       = total_biomass_tgc,
        total_area_ha   = total_area,
        TgC_per_Mha     = total_biomass_tgc / (total_area / 1e6),
        data            = df
    )
}

ssp_results <- list(
    SSP1 = compute_ssp_totals(future, "growth_SSP1_RCP19", pars),
    SSP2 = compute_ssp_totals(future, "growth_SSP2_RCP45", pars),
    SSP3 = compute_ssp_totals(future, "growth_SSP3_RCP70", pars)
)


# ── Pastureland baseline ──────────────────────────────────────────────────────

pars_no_flood <- pars[names(pars) != "floodable_forests"]

data_1k <- import_data("grid_1k_amazon_pastureland", biome = 1, n_samples = "all", categorical = categorical)
coords <- data_1k$coords
data_1k <- apply_min_max_scaling(data_1k$df, train_stats)

# pastureland does not have an age column, so we create one
# assuming it starts regrowing at 2015
data_1k$age <- 35

pred <- growth_curve(pars_no_flood, data_1k)
area <- data_1k[[grep("area", names(data_1k), value = TRUE)]] / 10000 # m2 → ha

# Biomass conversions
pred_tgc_per_ha <- pred * 1.258 * 0.5 / 1e6 # TgC/ha per pixel

# Select top pixels by biomass until cumulative area = SSP1 area
order_idx <- order(pred_tgc_per_ha, decreasing = TRUE)
cum_area <- cumsum(area[order_idx])
n_needed <- which(cum_area >= ssp_results$SSP1$total_area_ha)[1]
sel <- order_idx[seq_len(n_needed)]

total_pasture_tgc <- sum(pred_tgc_per_ha[sel] * area[sel], na.rm = TRUE)
total_pasture_area <- sum(area[sel], na.rm = TRUE) # ha


# ── Build summary data frame ──────────────────────────────────────────────────

ssp_summary <- data.frame(
    category = factor(
        c("SSP1", "SSP2", "SSP3", "Top Priority\nPasture"),
        levels = c("SSP1", "SSP2", "SSP3", "Top Priority\nPasture")
    ),
    total_TgC = c(
        ssp_results$SSP1$total_TgC,
        ssp_results$SSP2$total_TgC,
        ssp_results$SSP3$total_TgC,
        total_pasture_tgc
    ),
    total_area_Mha = c(
        ssp_results$SSP1$total_area_ha,
        ssp_results$SSP2$total_area_ha,
        ssp_results$SSP3$total_area_ha,
        total_pasture_area
    ) / 1e6
) %>%
    mutate(TgC_per_Mha = total_TgC / total_area_Mha)


# ── Reusable plot function ────────────────────────────────────────────────────

plot_carbon_bars <- function(df, y_col, y_label) {
    ggplot(df, aes(x = category, y = .data[[y_col]])) +
        geom_bar(stat = "identity", width = 0.7, fill = "#043927") +
        scale_y_continuous(
            labels = scales::label_comma(),
            name   = y_label
        ) +
        labs(x = NULL) +
        theme_minimal(base_size = 16) +
        theme(
            axis.title.y = element_text(size = 30, color = "black", margin = margin(r = 15)),
            axis.text.x  = element_text(size = 25, color = "black", margin = margin(t = 15)),
            axis.text.y  = element_text(size = 25, color = "black"),
            panel.grid   = element_blank(),
            axis.line    = element_line(color = "black", linewidth = 0.8)
        )
}


# ── Figure 4d: Total carbon (TgC) ────────────────────────────────────────────

fig_4d_total <- plot_carbon_bars(
    ssp_summary,
    y_col   = "total_TgC",
    y_label = "Total carbon stored by 2050 (TgC)"
)

ggsave("0_results/figures/figure_4d_total_carbon.jpeg",
    plot = fig_4d_total, width = 10, height = 12, dpi = 300
)


# ── Figure 4e: Carbon per unit area (TgC / Mha) ──────────────────────────────

fig_4e_per_area <- plot_carbon_bars(
    ssp_summary,
    y_col   = "TgC_per_Mha",
    y_label = "Carbon stored by 2050 (TgC / Mha)"
)

ggsave("0_results/figures/figure_4e_carbon_per_area.jpeg",
    plot = fig_4e_per_area, width = 10, height = 12, dpi = 300
)