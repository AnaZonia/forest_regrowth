library(pgirmess)
require(investr)
library(qpcR)
library(nlstools)

#Set the working directory will need to change for each user
setwd("/home/aavila/forest_regrowth/heinrich_poorter/Fig1_data_input/Fig1_data_input/")

########### Load the point data ####################

# Creating a name of data and reading the csv table
dados2<-read.csv2("MCWD_assessment_may2020_CHIRPSv1.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2

colnames(dados2)[1] <- "age"

mcwd_0<-dados2[dados2$threshold == 0,]
mcwd_180<-dados2[dados2$threshold == -180,]
mcwd_277<-dados2[dados2$threshold == -277,]
mcwd_350<-dados2[dados2$threshold == -350,] 

# Creating a name of data and reading the csv table
dados2<-read.csv2("fire_assessment_v2.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2
#dados2<-read.csv2("combined_summary_stats_no_prim_for_byMCWD.csv",sep=",")
# List the specific data/column of datafile
colnames(dados2)[1] <- "age"

fire0<-dados2[dados2$threshold == 0,]
fire1<-dados2[dados2$threshold == 1,]
fire2<-dados2[dados2$threshold == 2,]


# Creating a name of data and reading the csv table
dados2<-read.csv2("soil_assessment_v2.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2
#dados2<-read.csv2("combined_summary_stats_no_prim_for_byMCWD.csv",sep=",")
# List the specific data/column of datafile
colnames(dados2)[1] <- "age"

#note values are still in log format - in the Figure and the Supplementary Table 8 they have been converted: e.g 10^-0.1 = 0.8
SF1<-dados2[dados2$threshold == -1,] #10^-1 = 0.1
SF03<-dados2[dados2$threshold == -0.5,] #10^-0.5 = 0.3
SF05<-dados2[dados2$threshold == -0.3,] # 10^-0.3 = 0.5
SF08 = dados2[dados2$threshold == -0.1,] #10^-0.1 = 0.8

# Creating a name of data and reading the csv table
dados2<-read.csv2("deforestation_cycle_assessment_v2.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2
#dados2<-read.csv2("combined_summary_stats_no_prim_for_byMCWD.csv",sep=",")
# List the specific data/column of datafile
colnames(dados2)[1] <- "age"

def1<-dados2[dados2$threshold == 1,]
def2<-dados2[dados2$threshold == 2,]
def3<-dados2[dados2$threshold == 3,]

# Creating a name of data and reading the csv table
dados2<-read.csv2("radiation_assessment_v2.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2

colnames(dados2)[1] <- "age"

rad_150<-dados2[dados2$threshold == 150,]
rad_170<-dados2[dados2$threshold == 170,]
rad_180<-dados2[dados2$threshold == 180,]
rad_187<-dados2[dados2$threshold == 187,]

# Creating a name of data and reading the csv table
dados2<-read.csv2("precipitation_assessment_may2020_CHIRPS_thirds.csv",sep=",")
dados2$Corrected_AGB = dados2$Corrected_AGB/2
dados2$majority_fiftyth = dados2$majority_fiftyth/2

colnames(dados2)[1] <- "age"

precip_0<-dados2[dados2$threshold == 0,]
precip_1920<-dados2[dados2$threshold == 1920,]
precip_2210<-dados2[dados2$threshold == 2210,]

######## Load the regrowth models ##########################

load("regrowth_model_MCWD_0.RData")
load("regrowth_model_MCWD_180.RData")
load("regrowth_model_MCWD_277.RData")
load("regrowth_model_MCWD_350.RData")

load("regrowth_model_fire_0.RData")
load("regrowth_model_fire_1.RData")
load("regrowth_model_fire_2.RData")

load("regrowth_model_SF_1.RData")
load("regrowth_model_SF_03.RData")
load("regrowth_model_SF_05.RData")
load("regrowth_model_SF_08.RData")


load("regrowth_model_def_1.RData")
load("regrowth_model_def_2.RData")
load("regrowth_model_def_3.RData")

load("regrowth_model_rad_150.RData")
load("regrowth_model_rad_170.RData")
load("regrowth_model_rad_180.RData")
load("regrowth_model_rad_187.RData")

load("regrowth_model_precip_0.RData")
load("regrowth_model_precip_1920.RData")
load("regrowth_model_precip_2210.RData")


load("regrowth_model_MCWD_0.RData")
model_CR_350

#########Plotting #############
#jpeg('all_variables_byAGC_aug2020_fire_corrected_v2.jpg', res = 500, width = 7.5, height = 6, units = 'in')
pdf("Fig1a_f.pdf", width = 7.5, height = 6)
par(mfrow = c(3,2), mar=c(0.7, 0.5, 1.0, 0.2), oma = c(4, 5, 0.2, 0.2),  xpd = NA) # (bottom, left, top, right) #orig: mar=c(0.5, 0.5, 0.2, 0.2)
par(family = "sans")

plotFit(model_CR_rad187,data=rad_187,col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5), shade=TRUE, ylab='', xlab = '', yaxt="n", xaxt="n",tck =-0.03)
axis(side = 1, labels=F, tck = -0.03)
par(new = TRUE)
plotFit(model_CR_rad180,data=rad_180,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='', yaxt="n",  xaxt="n", tck =-0.03)
par(new = TRUE)
plotFit(model_CR_rad170, data=rad_170,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', yaxt="n", xaxt="n", tck =-0.03)
par(new = TRUE)
plotFit(model_CR_rad150,data=rad_150,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'dodgerblue3',col.conf=adjustcolor('dodgerblue1', 0.5), shade=TRUE, ylab=expression('AGC (Mg C ha'^-1*')'), las =1,xlab = '', xaxt="n", tck =-0.03, cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5,
        cex.sub = 1.5)
#text(x=(38), y =(20), labels=c("test"))
axis(side = 1, labels=F, tck = -0.03)
# Add a legend
legend("topleft", inset = 0.01,legend=c("Very low radiation (<170)", "Low radiation (170 to 180)", "Moderate radiation (180 to 187)", "High radiation (>187)"),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.15, pt.cex=1.15)

plotFit(model_CR_0,data=mcwd_0,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n", tck =-0.03)
par(new = TRUE)
plotFit(model_CR_180,data=mcwd_180,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n", yaxt="n",tck =-0.03)
par(new = TRUE)
plotFit(model_CR_277,data=mcwd_277,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='', xaxt="n",yaxt="n", tck =-0.03)
par(new = TRUE)
plotFit(model_CR_350,data=mcwd_350,col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.3), shade=TRUE, xlab='',  ylab='', xaxt="n", yaxt="n", tck = -0.03)
axis(side = 1, labels=F, tck = -0.03)
axis(side = 2, labels=F,  tck = -0.03)
legend("topleft", inset = 0.01,legend=c("Very low MCWD (>-180)", "Low MCWD (-180 to -277)", "Moderate MCWD (-277 to -350)", "Very high MCWD (<-350)"),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.15, pt.cex=1.15)

plotFit(model_CR_precip0,data=precip_0,col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorchid4',col.conf=adjustcolor('darkorchid', 0.5), shade=TRUE, ylab='', xlab = '', yaxt="n", xaxt="n", tck =-0.03)
axis(side = 1, labels=F, tck = -0.03)
par(new = TRUE)
#plotFit(model_CR_precip1920,data=precip_1920,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE, xlab = '', ylab='',yaxt="n", xaxt="n", tck =-0.03)
#par(new = TRUE)
plotFit(model_CR_precip1920, data=precip_1920,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', yaxt="n", xaxt="n",tck =-0.03)
par(new = TRUE)
plotFit(model_CR_precip2210,data=precip_2210,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'dodgerblue3',col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE, xlab='',  ylab=expression('AGC (Mg C ha'^-1*')'),xaxt="n", las =1, tck = -0.03,  cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5,
        cex.sub = 1.5)
axis(side = 1, labels=F, tck = -0.03)
legend("topleft", inset = 0.01,legend=c("High precipitation (>2210)", "Moderate precipitation (1920 to 2210)", "Low precipitation (<1920)"),
       col=c("dodgerblue", "darkseagreen3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.12, pt.cex=1.15)

plotFit(model_CR_SF1,data=SF1,col = 'darkorchid3 ',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'darkorchid4', col.conf=adjustcolor('darkorchid ', 0.3), shade=TRUE,  xlab = '', ylab='',xaxt="n", yaxt="n",tck = -0.03)
par(new = TRUE)
plotFit(model_CR_SF03,data=SF03,col = 'darkorange3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkorange4',col.conf=adjustcolor('darkorange', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",yaxt="n", tck = -0.03)
par(new = TRUE)
plotFit(model_CR_SF05,data=SF05,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE, xlab='',  ylab='', xaxt="n",yaxt="n", tck = -0.03)
par(new = TRUE)
plotFit(model_CR_SF08,data=SF08,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'dodgerblue3',col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE, xlab='',  ylab='',xaxt="n", yaxt="n", tck = -0.03)
axis(side = 1, labels=F, tck = -0.03)
axis(side = 2, labels=F, tck = -0.03)
# Add a legend
legend("topleft", inset = 0.01,legend=c("High SCC (>0.8)", "Moderate SCC (0.5 to 0.8)", "Low SCC (0.3 to 0.5)", "Very low SCC (< 0.3)"),
       col=c("dodgerblue", "darkseagreen3", "darkorange3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.15, pt.cex=1.15)

plotFit(model_CR_fire0,data=fire0,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='',  xaxt="n", yaxt="n", tck = -0.03)
par(new = TRUE)
plotFit(model_CR_fire1,data=fire1,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n",  yaxt="n", tck = -0.03)
par(new = TRUE)
plotFit(model_CR_fire2,data=fire2,col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'darkorchid4', col.conf=adjustcolor('darkorchid', 0.3), shade=TRUE,  xlab = 'Age (years)', ylab=expression('AGC (Mg C ha'^-1*')'), las =1,tck =-0.03, cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5,
        cex.sub = 1.5)
axis(side = 2, labels=F, tck = -0.03)
# Add a legend
legend("topleft", inset = 0.01,legend=c("Fire occurance: 0", "Fire occurance: 1", "Fire occurance: 2+"),
       col=c("dodgerblue", "darkseagreen3", "darkorchid3"), pch = 19, box.lty = 0, bg="transparent", cex=1.15, pt.cex=1.15)


plotFit(model_CR_def1,data=def1,col = 'dodgerblue',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'dodgerblue3', col.conf=adjustcolor('dodgerblue1', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n", yaxt="n" ,tck = -0.03)
par(new = TRUE)
plotFit(model_CR_def2,data=def2,col = 'darkseagreen3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit= 'darkseagreen4',col.conf=adjustcolor('darkseagreen', 0.3), shade=TRUE,  xlab = '', ylab='', xaxt="n", yaxt="n", tck = -0.03)
par(new = TRUE)
plotFit(model_CR_def3,data=def3,col = 'darkorchid3',interval="confidence",xlim=c(0,40),ylim=c(0,110),pch=19,col.fit = 'darkorchid4', col.conf=adjustcolor('darkorchid', 0.3), shade=TRUE,  xlab = 'Age (years)', ylab='', yaxt="n", tck =-0.03, cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5,
        cex.sub = 1.5)
axis(side = 2, labels=F, tck = -0.03)
legend("topleft", inset = 0.01,legend=c("Deforested: 1", "Deforested: 2", "Deforested: 3+"),
       col=c("dodgerblue", "darkseagreen3", "darkorchid3"), pch = 19, box.lty = 0, cex=1.15, pt.cex=1.15)

dev.off()




