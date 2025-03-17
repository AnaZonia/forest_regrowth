library(ggplot2)

tst <- function(k0){
    # Define parameters
    asymptote <- 200
    lag <- 30
    ages <- 0:100
    theta <- 5
    k <- 0.05

    # Compute values for full growth curve
    values <- asymptote * (1 - exp(-k * ages))^theta

    # Extract sampled values
    ages_sampled <- lag:50
    values_sampled <- values[lag:50]

    values_k0_added <- asymptote * (1 - exp(-(k0 + k * ages)))^theta
    # this changes the intercept and rate of growth
    # makes sense why lag shouldn't matter in this case - we're basically introducing an intercept term, so erasing it.
    # but how come it was working before????

    values_k0_multiplied <- asymptote * (1 - exp(-((k0 + k) * ages)))^theta
    # this only changes the rate of growth
    print(ages)

    extrapolated_ages <- (ages_sampled[1] - lag):ages_sampled[1] # 0:lag
    extrapolated_values <- asymptote * (1 - exp(-((k0 + k) * extrapolated_ages)))^theta
    
    print(extrapolated_ages)
    print(extrapolated_values[lag])
    print(values_k0_multiplied[100])
    
    # Plot using ggplot2
    q <- ggplot() +
        geom_line(aes(x = ages, y = values), color = "black", size = 1.2) + # Full curve
        geom_line(aes(x = ages_sampled, y = values_sampled), color = "red", size = 1.5) + # Sampled section
        geom_line(aes(x = ages, y = values_k0_added), color = "green", size = 1.5) + # Sampled section
        # geom_line(aes(x = ages, y = values_k0_multiplied), color = "purple", size = 1.5) + # Sampled section
        geom_line(aes(x = extrapolated_ages, y = extrapolated_values), color = "blue", size = 1.5) + # Extrapolated section
        labs(title = "Growth Curve with Delay", x = "Time", y = "Population Size") +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    print(q)

}

tst(0.3)


