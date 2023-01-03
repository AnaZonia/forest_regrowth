


agb_forest_age <- cbind(regrowth_cleaned, agbd = GEDI[match(regrowth_cleaned$xy,GEDI$xy),c("agbd")])
agb_forest_age <- agb_forest_age[complete.cases(agb_forest_age[, ncol(agb_forest_age)]), ]

plot(agb_forest_age$forest_age, agb_forest_age$agbd)

# saveRDS(agb_forest_age, 'agb_forest_age.rds')
# agb_forest_age = readRDS('agb_forest_age.rds')

sds <- aggregate(agbd ~ forest_age, agb_forest_age, sd)
means <- aggregate(agbd ~ forest_age, agb_forest_age, mean)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'mean', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = mean)) + 
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd)) +
  geom_point() + theme(text = element_text(size = 20))  

