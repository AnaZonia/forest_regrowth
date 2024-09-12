    # ---------------------
    library(fastDummies)

    # Select only the dummy columns for PCA
    pca_cols <- grep("ecoreg_|soil_", names(data), value = TRUE)
    pca_cols <- c(
        grep("ecoreg_|soil_", names(data), value = TRUE)
        # "cwd", "mean_prec", "mean_si"
    )
    pca_data <- data[, pca_cols]
    pca_result <- prcomp(pca_data, scale. = TRUE)
    summary(pca_result)
    # scree plot
    plot(pca_result, type = "l")
    biplot(pca_result)
    eigenvalues <- (pca_result$sdev)^2
    print(eigenvalues)

    cumulative_var <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
    plot(cumulative_var, type = "b", ylab = "Cumulative Proportion of Variance Explained", xlab = "Number of Principal Components")

    n_components <- 4 # number of components you decided to keep
    pca_scores <- pca_result$x[, 1:n_components]
    head(pca_scores)
    cols_to_remove <- c(grep("ecoreg_|soil_", names(data), value = TRUE), "mean_prec", "num_fires_after_regrowth", "ts_fire_after_regrowth")
    cols_to_remove
    pars_names <- setdiff(pars_names, c("ecoreg", "soil", "mean_prec", "num_fires_after_regrowth", "ts_fire_after_regrowth"))
    pars_names
    data <- cbind(data[, setdiff(names(data), cols_to_remove)], pca_scores)
