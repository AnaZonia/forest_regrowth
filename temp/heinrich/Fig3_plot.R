library(pgirmess)
require(investr)
library(qpcR)
library(nlstools)


#Set the working directory will need to change for each user
setwd("D:\\Brazil_2020\\code_for_paper\\Data_jan2021\\Fig3_data_input\\")


########### Load the point data ####################
data1 =read.csv2("SW-Group.csv",sep=",")
data2 =read.csv2("SE-Group.csv",sep=",")
data3 =read.csv2("NW-Group.csv",sep=",")
data4 =read.csv2("NE-Group.csv",sep=",")

# class1: SW region 
colnames(data1)[1] <- "age"
neither_1<-data1[data1$threshold == 11,]
burnt_1<-data1[data1$threshold == 12,]
deforest_1<-data1[data1$threshold == 21,]
both_1<-data1[data1$threshold == 22,]

neither_1$majority_fiftyth = neither_1$majority_fiftyth/2
neither_1$Corrected_AGB = neither_1$Corrected_AGB/2

burnt_1$majority_fiftyth = burnt_1$majority_fiftyth/2
burnt_1$Corrected_AGB = burnt_1$Corrected_AGB/2

deforest_1$majority_fiftyth = deforest_1$majority_fiftyth/2
deforest_1$Corrected_AGB = deforest_1$Corrected_AGB/2

both_1$majority_fiftyth = both_1$majority_fiftyth/2
both_1$Corrected_AGB = both_1$Corrected_AGB/2

# class2: SE regions 

colnames(data2)[1] <- "age"
neither_2<-data2[data2$threshold == 11,]
burnt_2<-data2[data2$threshold == 12,]
deforest_2<-data2[data2$threshold == 21,]
both_2<-data2[data2$threshold == 22,]

neither_2$majority_fiftyth = neither_2$majority_fiftyth/2
neither_2$Corrected_AGB = neither_2$Corrected_AGB/2

burnt_2$majority_fiftyth = burnt_2$majority_fiftyth/2
burnt_2$Corrected_AGB = burnt_2$Corrected_AGB/2

deforest_2$majority_fiftyth = deforest_2$majority_fiftyth/2
deforest_2$Corrected_AGB = deforest_2$Corrected_AGB/2

both_2$majority_fiftyth = both_2$majority_fiftyth/2
both_2$Corrected_AGB = both_2$Corrected_AGB/2

# class3: NW region 

colnames(data3)[1] <- "age"
neither_3<-data3[data3$threshold == 11,]
burnt_3<-data3[data3$threshold == 12,]
deforest_3<-data3[data3$threshold == 21,]
both_3<-data3[data3$threshold == 22,]

neither_3$majority_fiftyth = neither_3$majority_fiftyth/2
neither_3$Corrected_AGB = neither_3$Corrected_AGB/2

burnt_3$majority_fiftyth = burnt_3$majority_fiftyth/2
burnt_3$Corrected_AGB = burnt_3$Corrected_AGB/2

deforest_3$majority_fiftyth = deforest_3$majority_fiftyth/2
deforest_3$Corrected_AGB = deforest_3$Corrected_AGB/2

both_3$majority_fiftyth = both_3$majority_fiftyth/2
both_3$Corrected_AGB = both_3$Corrected_AGB/2

# class4: NE region

colnames(data4)[1] <- "age"
neither_4<-data4[data4$threshold == 11,]
burnt_4<-data4[data4$threshold == 12,]
deforest_4<-data4[data4$threshold == 21,]
both_4<-data4[data4$threshold == 22,]

neither_4$majority_fiftyth = neither_4$majority_fiftyth/2
neither_4$Corrected_AGB = neither_4$Corrected_AGB/2

burnt_4$majority_fiftyth = burnt_2$majority_fiftyth/2
burnt_4$Corrected_AGB = burnt_4$Corrected_AGB/2

deforest_4$majority_fiftyth = deforest_4$majority_fiftyth/2
deforest_4$Corrected_AGB = deforest_4$Corrected_AGB/2

both_4$majority_fiftyth = both_4$majority_fiftyth/2
both_4$Corrected_AGB = both_4$Corrected_AGB/2

sum_maj_1 = sum(data1$majority_count)
sum_maj_2 = sum(data2$majority_count)
sum_maj_3 = sum(data3$majority_count)
sum_maj_4 = sum(data4$majority_count)
sum_all = sum(sum_maj_1, sum_maj_2, sum_maj_3, sum_maj_4)
######## Load the regrowth models ##########################

load("regrowth_model_SW_NoDist.RData")
load("regrowth_model_SW_Fire.RData")
load("regrowth_model_SW_Deforest.RData")
load("regrowth_model_SW_BothDist.RData")

load("regrowth_model_SE_NoDist.RData")
load("regrowth_model_SE_Fire.RData")
load("regrowth_model_SE_Deforest.RData")
load("regrowth_model_SE_BothDist.RData")

load("regrowth_model_NW_NoDist.RData")
load("regrowth_model_NW_Fire.RData")
load("regrowth_model_NW_Deforest.RData")
load("regrowth_model_NW_BothDist.RData")

load("regrowth_model_NE_NoDist.RData")
load("regrowth_model_NE_Fire.RData")
load("regrowth_model_NE_Deforest.RData")
load("regrowth_model_NE_BothDist.RData")


#########Plotting #############

# SOUTH-WEST
neither_legend_1 = paste0("Single deforestation, Not burnt (" ,round(sum(neither_1$majority_count)/sum_maj_1*100, 1), "%)")
burnt_legned_1 = paste0("Single deforestation, Burnt (" ,round(sum(burnt_1$majority_count)/sum_maj_1*100, 1), "%)")
deforest_legned_1 = paste0("Repeated deforestation, Not burnt (" ,round(sum(deforest_1$majority_count)/sum_maj_1*100, 1), "%)")
both_legned_1 = paste0("Repeated deforestation, Burnt (" ,round(sum(both_1$majority_count)/sum_maj_1*100, 1), "%)")


pdf("region_SW_byAGC_regrowth.pdf", width = 7.5, height = 6)
par( mar=c(5,6,4,1)+.05)
plotFit(CR_1_1,data=neither_1,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '',  xaxt="n",yaxt="n", ylab='')
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_1_2,data=burnt_1,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='',  xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_1_3,data=deforest_1,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='',  xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_1_4,data=both_1[1:29,],col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5),
        las =1,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.25,
        cex.sub = 1.5,
        shade=TRUE,  xlab='Age (years)',  ylab=expression('AGC (Mg C ha'^-1*')'), font.main=1, main = paste0("South-West & Central region (",round((sum_maj_1/sum_all)*100, 1), "% of SF)"))
# Add a legend
legend("topleft", inset = 0.01,legend=(c(neither_legend_1, burnt_legned_1, deforest_legned_1, both_legned_1)),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.25, pt.cex=1.25)
mtext("(Moderate MCWD - Low radiation - Moderate Precipitation - High SCC)", cex=1.1)
dev.off()

neither_legend_2 = paste0("Single deforestation, Not burnt (" ,round(sum(neither_2$majority_count)/sum_maj_2*100, 1), "%)")
burnt_legned_2 = paste0("Single deforestation, Burnt (" ,round(sum(burnt_2$majority_count)/sum_maj_2*100, 1), "%)")
deforest_legned_2 = paste0("Repeated deforestation, Not burnt (" ,round(sum(deforest_2$majority_count)/sum_maj_2*100, 1), "%)")
both_legned_2 = paste0("Repeated deforestation, Burnt (" ,round(sum(both_2$majority_count)/sum_maj_2*100, 1), "%)")



#jpeg('region_SE_byAGC_v2_old_growth_forest.jpg', res = 500, width = 7.5, height = 6, units = 'in')
pdf("region_SE_byAGC_regrowth.pdf", width = 7.5, height = 6)
par(new = TRUE, mar=c(5,6,4,1)+.1)
plotFit(CR_2_1,data=neither_2,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.1)
plotFit(CR_2_2,data=burnt_2,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.1)
plotFit(CR_2_3,data=deforest_2,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.1)
plotFit(CR_2_4,data=both_2[1:29,],col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5), 
        las =1,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.25,
        cex.sub = 1.5,
        shade=TRUE,  xlab='Age (years)',  ylab=expression('AGC (Mg C ha'^-1*')'), font.main=1, main = paste0("South-East & North region (",round((sum_maj_2/sum_all)*100, 1), "% of SF)"))
# Add a legend
legend("topleft", inset = 0.01,legend=c(neither_legend_2, burnt_legned_2, deforest_legned_2, both_legned_2),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, cex=1.25, pt.cex=1.25)
mtext("(Very high MCWD - Moderate radiation - Low Precipitation - Moderate SCC)" ,cex=1.1)
dev.off()



neither_legend_3 = paste0("Single deforestation, Not burnt (" ,round(sum(neither_3$majority_count)/sum_maj_3*100, 1), "%)")
burnt_legned_3 = paste0("Single deforestation, Burnt (" ,round(sum(burnt_3$majority_count)/sum_maj_3*100, 1), "%)")
deforest_legned_3 = paste0("Repeated deforestation, Not burnt (" ,round(sum(deforest_3$majority_count)/sum_maj_3*100, 1), "%)")
both_legned_3 = paste0("Repeated deforestation, Burnt (" ,round(sum(both_3$majority_count)/sum_maj_3*100, 1), "%)")




#jpeg('region_NW_byAGC_old_growth_forest.jpg', res = 500, width = 7.5, height = 6, units = 'in')
pdf("region_NW_byAGC_regrowth.pdf", width = 7.5, height = 6)
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_3_1,data=neither_3,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_3_2,data=burnt_3,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_3_3,data=deforest_3,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_3_4,data=both_3[1:29,],col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5), 
        las =1,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.25,
        cex.sub = 1.5,
        shade=TRUE,  xlab='Age (years)',  ylab=expression('AGC (Mg C ha'^-1*')'),font.main=1, main = paste0("North-West region (",round((sum_maj_3/sum_all)*100, 1), "% of SF)"))
# Add a legend
legend("topleft", inset = 0.01,legend=c(neither_legend_3, burnt_legned_3, deforest_legned_3, both_legned_3),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.25, pt.cex=1.25)
mtext("(Very Low MCWD - Low radiation - High Precipitation - Low SCC)",cex=1.1)
dev.off()





neither_legend_4 = paste0("Single deforestation, Not burnt (" ,round(sum(neither_4$majority_count)/sum_maj_4*100, 1), "%)")
burnt_legned_4 = paste0("Single deforestation, Burnt (" ,round(sum(burnt_4$majority_count)/sum_maj_4*100, 1), "%)")
deforest_legned_4 = paste0("Repeated deforestation, Not burnt (" ,round(sum(deforest_4$majority_count)/sum_maj_4*100, 1), "%)")
both_legned_4 = paste0("Repeated deforestation, Burnt (" ,round(sum(both_4$majority_count)/sum_maj_4*100, 1), "%)")

#jpeg('region_NE_byAGC_old_growth_forest.jpg', res = 500, width = 7.5, height = 6, units = 'in')
pdf("region_NE_byAGC_regrowth.pdf", width = 7.5, height = 6)
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_4_1,data=neither_4,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n", yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_4_2,data=burnt_4,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='',xaxt="n",yaxt="n")
par(new = TRUE, mar=c(5,6,4,1)+.05)
plotFit(CR_4_3,data=deforest_4,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='', xaxt="n",yaxt="n")
par(new = TRUE)
plotFit(CR_4_4,data=both_4[1:29,],col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5), 
        las =1,
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.25,
        cex.sub = 1.5,
        shade=TRUE,  xlab='Age (years)', ylab=expression('AGC (Mg C ha'^-1*')'),font.main=1, main = paste0("North-East & Central-North (",round((sum_maj_4/sum_all)*100, 1), "% of SF)"))
# Add a legend
legend("topleft", inset = 0.01,legend=c(neither_legend_4, burnt_legned_4, deforest_legned_4, both_legned_4),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, cex=1.25, pt.cex=1.25)
mtext("(Low MCWD - High radiation - High Precipitation - Low SCC)",cex=1.1)
dev.off()


