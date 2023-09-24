library(caret) 
library(caTools)
library(randomForest)
library(party)
library(foreign)
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(egg)

windowsFonts(Arial=windowsFont("TT Arial"))

set.seed(5)
random_iterations=floor(runif(30, min=0, max=1000))
random_iterations


conditional_importance = function(region_str, conditional_string){ 

  
  for (j in 1:length(random_iterations)){
  
  i = random_iterations[j]
  
  set.seed(i)
  

  if (region_str == "whole_Amazon"){
    folder = paste0(main_folder, region_str,"_sector")
    setwd(folder)
    filename = paste0(region_str , "_2500_cforest_importance_conditional",conditional_string,"_seed",i,".RDS") }
  else {
    folder = paste0(main_folder, "sample_size_1300\\",region_str,"_sector")
    setwd(folder)
    filename = paste0(region_str , "_cforest_importance_conditional",conditional_string,"_seed",i,".RDS")
    }
  
  print(paste(Sys.time(), "starting with"))
  print(paste0(filename))
  
  data = readRDS(filename)
  imp = data$importance
  imp$rank = rank(imp$Overall)
  
  imp = cbind(variables = rownames(imp), imp)
  
  datalist[[j]] = imp

}

variables = (
  c(
    "Forest age",
    "MCWD",
    "SCC",
    "Deforestation",
    "Fire",
    "Precipitation",
    "SW radiation"
  )
)
av_rank = rowMeans(sapply(datalist, "[[", "rank"))
standard_dev = rowSds(sapply(datalist, "[[", "rank"))
margin = qt(0.95, df=30)*standard_dev/sqrt(30)
stat_data = data.frame("Variable" = variables,"Average_rank" = av_rank,"CI" = margin)

row.names(stat_data) = variables

par()
variable_plot_order = rev(
  c(
    "Deforestation",
    "Fire",
    "SCC",
    "Precipitation",
    "MCWD",
    "SW radiation",
    "Forest age")
  )
#determine the color
if (region_str == "SW") {
  color_shade = "#06f001"
} else if (region_str == "SE") {
  color_shade = "#1704a0"
} else if (region_str == "NW") {
  color_shade = "#059a01"
} else if (region_str == "NE") {
  color_shade = "#05a9ff"
} else if (region_str =="whole_Amazon") {
  color_shade = "gray50"
}
  stat_data$variable_plot_order = factor(stat_data$Variable, levels =
                                           variable_plot_order)
  
  rank_graph = ggplot(stat_data, aes(x = variable_plot_order, y = Average_rank)) +
    #annotate("rect", xmin = -Inf, xmax = 5.5, ymin = -Inf, ymax = Inf,
     #        alpha = 0.2) + 
    #annotate("rect", xmin = 5.5, xmax = Inf, ymin = -Inf, ymax = Inf,
     #        alpha = 0.5) +
    #annotate("text", x = 3, y = 7.5, label = "Environmental \n drivers") +
    #annotate("text", x = 6.6, y = 7.5, label = "Disturbance (human) \n drivers") +
    geom_bar(stat = "identity", fill = color_shade, alpha = 0.7) +
    geom_errorbar(aes(ymin=Average_rank-CI, ymax=Average_rank+CI, width=0.15)) +
    ylim(0, 8) +
    scale_x_discrete(breaks = variable_plot_order) +
    geom_text(aes(
      label = (round(Average_rank, digits = 1)),
      hjust = 0.5
    ),
    size = 3, position = position_stack(vjust = 0.5)) + #vjust=0, hjust=-0.15
    labs(x = "") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        color = "black",
        size = 10,
        angle = 90, vjust = 0.5, hjust = 1
      ),
      axis.text.y = element_text(color = "black", size = 8),
      text = element_text(size = 8)
    ) +
    labs(y = expression("Least important " %->% " Most important"), theme(vjust =5,axis.title.y = element_text(size = 8))) + 
    geom_vline(xintercept=c(5.5), linetype= 2)

  return (rank_graph)
}

main_folder = "D:\\Brazil_2020\\code_for_paper\\Data_jan2021\\Fig1g_2b_e_variable_importance\\multiple_iterations\\"
setwd(main_folder)

datalist = list()
whole_Amazon = conditional_importance("whole_Amazon", "False")

whole_Amazon

Amz = ggarrange(whole_Amazon)
annotated_amz = annotate_figure(Amz,left = text_grob("Average Ranking", rot = 90, hjust=0.0, size=10),
                                bottom = text_grob("Drivers", size=10, vjust =-1.3) )

annotated_amz

ggsave(plot = annotated_amz , width = 5, height = 11, unit="cm", dpi = 1200, filename = "whole_amazon_var_importanceV2.pdf")

dev.off()

datalist = list()

setwd(main_folder)

NW_graph = conditional_importance("NW", "True")
NE_graph = conditional_importance("NE", "True")
SW_graph = conditional_importance("SW", "True")
SE_graph = conditional_importance("SE", "True")


figure = ggarrange(NW_graph +
                     theme(plot.margin = ggplot2::margin(r = 1)),
                   NE_graph + 
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           plot.margin = ggplot2::margin(r = 1, l = 1)),
                   SW_graph +
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           plot.margin = ggplot2::margin(r = 1, l = 1)),
                   SE_graph +
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           plot.margin = ggplot2::margin(r = 1, l = 1)),
                   ncol=4, nrow=1)

annotated=annotate_figure(figure,
                left = text_grob("Average Ranking", rot = 90, hjust=0.0, size=10),
                bottom = text_grob("Drivers", size=10, vjust =-1.3))


annotated

setwd(main_folder)

ggsave(plot = annotated ,width = 17, height = 9,unit="cm",dpi = 1200,filename = "regional_var_importanceV2.pdf")

dev.off()
