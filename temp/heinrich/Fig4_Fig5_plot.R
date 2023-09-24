library(raster)
library(ggpubr)
library(foreign)
library(ggplot2)
library(sf)
library(rgdal)
library(dplyr)
#library(ggforce)
#library(stringr)
#library(reshape)
#library(viridis)
#library(ggsci)
#library(erer)

#NOTE: FOR THIS SCRIPT TO WORK - NEED TO RUN "plot_byAllRegions_as_AGB.R.R" first!!! ###############
#This is to initiate the growth models. 
#script contains reference to AGB - but this is later converted to AGCarbon! 
#Where AGB_rate refers to the AGB at that particular age - rather than a growth rate. 

########## USER - please change following three lines to your directory ############
setwd("D:\\Brazil_2020\\code_for_paper\\Data\\Fig3_data_input\\")
#read in the source script - ensure to set working directroy above to following folder: "Fig3_data_input"
source("D:\\Brazil_2020\\code_for_paper\\Data\\Fig3_data_input\\load_AGB_regrowth_models.R")

setwd("D:\\Brazil_2020\\code_for_paper\\Data\\Fig4_5_carbon_sink_2017\\secondary_forest_by_region_and_disturbance\\")

pattern_name = paste0('sec_for_2016_.*ist.tif.vat.dbf$')
all_files = list.files(pattern = pattern_name, full.names = T)
all.data_2016 = str_remove(all_files, "AGB_")
all.data_2017 = str_remove(all_files, "2016_")
all.data_2017 = str_remove(all.data_2017, "AGB_")
all.data_2016 = unique(all.data_2016)
all.data_2017 = unique(all.data_2017)

###Only need to use this function once ######
###This function applies the correct growth model depending on region and disturbance type
#and will save this infomation as a column in the dataframe (dbf) of the file. 
apply_models = function(all.data){
  
  datalist = list()
  for (i in all.data){
    
    
    data_dbf = foreign::read.dbf((i))
    
    #if(grepl("2016", i) == FALSE){
     # data_dbf = data_dbf[-1,]
   # }
    
  
    data_dbf$area_hect = data_dbf$Count * 0.09
    data_dbf$AGB_f = fearnside[2:(nrow(data_dbf)+ 1), 3]
    #data_dbf$AGB_Fearn = data_dbf$AGB_f * data_dbf$area_hect
    data_dbf$AGB_Fearn = data_dbf$Value * 3.94 *data_dbf$area_hect 
    data_dbf$AGB_wang =  ((341.2*data_dbf$Value)/(35 + data_dbf$Value)) * data_dbf$area_hect #341.2
    
    if (grepl("SW", i) ==TRUE){
      print(paste0("this is SW region"))
      age = data.frame(age=seq(1, max(data_dbf$Value), 1))
      data_dbf$no_dist = predict(CR_1_1, age) * data_dbf$area_hect
      
      
      if (grepl("noDist", i) == TRUE){
        #to extract data as tiff later - need to conver to interger -
        #integers are less heavy
        #So I have multiplied by 10000 to retrain detail of growth - need to divide by 10,000 later on again.
        data_dbf$AGB_rate = as.integer(predict(CR_1_1, age) *10000)  
                                       
        data_dbf$AGB_total = predict(CR_1_1, age) * data_dbf$area_hect
        
        
      } else if (grepl("fire", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_1_2, age) *10000) 
        data_dbf$AGB_total = predict(CR_1_2, age) * data_dbf$area_hect
        

      }  else if (grepl("deforest", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_1_3, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_1_3, age) * data_dbf$area_hect

      } else if (grepl("both", i) == TRUE){
        data_dbf$AGB_rate = as.integer(predict(CR_1_4, age) *10000)
        
        data_dbf$AGB_total = predict(CR_1_4, age) * data_dbf$area_hect
        

      }
    } else if (grepl("SE", i) ==TRUE){
      print(paste0("this is SE region"))
      age = data.frame(age=seq(1, max(data_dbf$Value), 1))
      data_dbf$no_dist = predict(CR_2_1, age) * data_dbf$area_hect
      

      
      if (grepl("noDist", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_2_1, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_2_1, age) * data_dbf$area_hect
        

      } else if (grepl("fire", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_2_2, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_2_2, age) * data_dbf$area_hect
        

      }  else if (grepl("deforest", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_2_3, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_2_3, age) * data_dbf$area_hect
        

      } else if (grepl("both", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_2_4, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_2_4, age) * data_dbf$area_hect
        

      }
    } else if (grepl("NW", i) ==TRUE){
      print(paste0("this is NW region"))
      age = data.frame(age=seq(1, max(data_dbf$Value), 1))
      data_dbf$no_dist = predict(CR_3_1, age) * data_dbf$area_hect
      

      
      if (grepl("noDist", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_3_1, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_3_1, age) * data_dbf$area_hect
        

      } else if (grepl("fire", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_3_2, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_3_2, age) * data_dbf$area_hect
        

      }  else if (grepl("deforest", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_3_3, age) *10000)
        
        data_dbf$AGB_total = predict(CR_3_3, age) * data_dbf$area_hect
        

      } else if (grepl("both", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_3_4, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_3_4, age) * data_dbf$area_hect

      }
    } else if (grepl("NE", i) ==TRUE){
      print(paste0("this is NE region"))
      
      age = data.frame(age=seq(1, max(data_dbf$Value), 1))
      data_dbf$no_dist = predict(CR_4_1, age) * data_dbf$area_hect
      
      
      if (grepl("noDist", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_4_1, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_4_1, age) * data_dbf$area_hect
        
      } else if (grepl("fire", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_4_2, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_4_2, age) * data_dbf$area_hect
        
      }  else if (grepl("deforest", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_4_3, age) *10000) 
        
        data_dbf$AGB_total = predict(CR_4_3, age) * data_dbf$area_hect
        

      } else if (grepl("both", i) == TRUE){
        
        data_dbf$AGB_rate = as.integer(predict(CR_4_4, age) *10000)  
        data_dbf$AGB_total = predict(CR_4_4, age) * data_dbf$area_hect
        

      }
    }
    
    #foreign::write.dbf(data_dbf, i)
    datalist[[i]] = data_dbf
  }
  return(datalist)
}


predict_future = function(all.data){
  
  future_list = list()
  for (k in 1:16){
    j = all.data[k]
    
    data_dbf = foreign::read.dbf((j))
    
    print(paste0(j))
    
    age = data.frame(age=seq(1, max(data_dbf$Value), 1))
    
    output = matrix(ncol=14, nrow= max(data_dbf$Value))
    
    
    for (i in 0:13){
      
      print(paste0(i))
      #data_dbf = foreign::read.dbf((file))
      
      if (grepl("SW", j) ==TRUE){
        print(paste0("this is SW region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))

        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){

          data_dbf$AGB.j = predict(CR_1_1, ageX) * data_dbf$area_hect

        } else if (grepl("fire", j) == TRUE){
          data_dbf$AGB.j = predict(CR_1_2, ageX) * data_dbf$area_hect
        }  else if (grepl("deforest", j) == TRUE){
          data_dbf$AGB.j = predict(CR_1_3, ageX) * data_dbf$area_hect
        } else if (grepl("both", j) == TRUE){
          data_dbf$AGB.j = predict(CR_1_4, ageX) * data_dbf$area_hect
        }
      } else if (grepl("SE", j) ==TRUE){
        print(paste0("this is SE region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){
          data_dbf$AGB.j = predict(CR_2_1, ageX) * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          data_dbf$AGB.j = predict(CR_2_2, ageX) * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          data_dbf$AGB.j = predict(CR_2_3, ageX) * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          data_dbf$AGB.j = predict(CR_2_4, ageX) * data_dbf$area_hect
        }
      } else if (grepl("NW", j) ==TRUE){
        print(paste0("this is NW region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){
          data_dbf$AGB.j = predict(CR_3_1, ageX) * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          data_dbf$AGB.j = predict(CR_3_2, ageX) * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          data_dbf$AGB.j = predict(CR_3_3, ageX) * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          data_dbf$AGB.j = predict(CR_3_4, ageX) * data_dbf$area_hect
        }
      } else if (grepl("NE", j) ==TRUE){
        print(paste0("this is NE region"))
        
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        
        if (grepl("noDist", j) == TRUE){
          data_dbf$AGB.j = predict(CR_4_1, ageX) * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          data_dbf$AGB.j = predict(CR_4_2, ageX) * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          data_dbf$AGB.j = predict(CR_4_3, ageX) * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          data_dbf$AGB.j = predict(CR_4_4, ageX) * data_dbf$area_hect
          
        }
      }
      
      
      output[,i+1] = data_dbf$AGB.j
      
    }
  
    dataframe_model = data.frame("age"=age, output)
    future_list[[k]] = dataframe_model
  }
  
  return(future_list)
}


predict_future_confidence_interval = function(all.data){
  
  future_list_lower = list()
  future_list_upper = list()
  for (k in 1:16){
    j = all.data[k]
    
    data_dbf = foreign::read.dbf((j))
    
    print(paste0(j))
    
    age = data.frame(age=seq(1, max(data_dbf$Value), 1))
    
    output_lower = matrix(ncol=14, nrow= max(data_dbf$Value))
    output_upper = matrix(ncol=14, nrow= max(data_dbf$Value))
    
    for (i in 0:13){
      
      print(paste0(i))
      #data_dbf = foreign::read.dbf((file))
      
      if (grepl("SW", j) ==TRUE){
        print(paste0("this is SW region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        
        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){
          confidence_file = read.csv("SW_1_1_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          confidence_file = read.csv("SW_1_2_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          confidence_file = read.csv("SW_1_3_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          confidence_file = read.csv("SW_1_4_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }
      } else if (grepl("SE", j) ==TRUE){
        print(paste0("this is SE region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){
          confidence_file = read.csv("SE_2_1_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          confidence_file = read.csv("SE_2_2_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          confidence_file = read.csv("SE_2_3_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          confidence_file = read.csv("SE_2_4_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
        }
      } else if (grepl("NW", j) ==TRUE){
        print(paste0("this is NW region"))
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        if (grepl("noDist", j) == TRUE){
          confidence_file = read.csv("NW_3_1_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          confidence_file = read.csv("NW_3_2_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          confidence_file = read.csv("NW_3_3_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
          
        } else if (grepl("both", j) == TRUE){
          confidence_file = read.csv("NW_3_4_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          }
      } else if (grepl("NE", j) ==TRUE){
        print(paste0("this is NE region"))
        
        age = data.frame(age=seq(1, max(data_dbf$Value), 1))
        ageX = data.frame(age=age$age + i)
        
        
        if (grepl("noDist", j) == TRUE){
          confidence_file = read.csv("NE_4_1_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("fire", j) == TRUE){
          confidence_file = read.csv("NE_4_2_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }  else if (grepl("deforest", j) == TRUE){
          confidence_file = read.csv("NE_4_3_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        } else if (grepl("both", j) == TRUE){
          confidence_file = read.csv("NE_4_4_confidence_interval_as_AGC.csv")
          lower = confidence_file$lower*2
          upper = confidence_file$upper*2
          data_dbf$AGB.lower = lower[min(ageX):max(ageX)] * data_dbf$area_hect
          data_dbf$AGB.upper = upper[min(ageX):max(ageX)] * data_dbf$area_hect
          
        }
      }
      
      
      output_lower[,i+1] = data_dbf$AGB.lower
      output_upper[,i+1] = data_dbf$AGB.upper
      
    }
    
    dataframe_model_lower = data.frame("age"=age, output_lower)
    future_list_lower[[k]] = dataframe_model_lower
    
    dataframe_model_upper = data.frame("age"=age, output_upper)
    future_list_upper[[k]] = dataframe_model_upper
  }
  
  outputs = list("lower"=future_list_lower, "upper"= future_list_upper)
  return(outputs)
}

# datalist_2016 = apply_models(all.data_2016)
# all_data_2016 = do.call(rbind, datalist_2016)
# all_sums_2016=colSums(all_data_2016, na.rm = TRUE)
# 
# datalist_2017 = apply_models(all.data_2017)
# all_data_2017 = do.call(rbind, datalist_2017)
# all_sums_2017=colSums(all_data_2017, na.rm = TRUE)


#####work from here on now #######

#create empty lists to append to
changelist = list()
list_dbf_2016 = list()
list_dbf_2017 = list()
new_growth = list()


for (i in 1:length(all.data_2016)){
  
  file_2016 = all.data_2016[i]
  file_2017 =all.data_2017[i]
  dbf_2016 = foreign::read.dbf((file_2016))
  dbf_2017 = foreign::read.dbf((file_2017))
  dbf_2016$name = file_2016
  dbf_2017$name = file_2017
  
  pixel_loss = dbf_2017$Count[-1] - dbf_2016$Count
  #need to divide by 10,000 because I *10,000 earlier to keep data in interger format. Silly Me. 
  AGB_deforestation_loss = pixel_loss * (dbf_2016$AGB_rate/10000) * 0.09 
  net_change = dbf_2017$AGB_total[-1] - dbf_2016$AGB_total
  
  year_change = data.frame("pixel.loss" = pixel_loss, "AGB.deforestation.loss" = AGB_deforestation_loss, "net.change" = net_change)
  changelist[[i]] = year_change
  list_dbf_2016[[i]] = dbf_2016
  list_dbf_2017[[i]] = dbf_2017
  new_growth[[i]] = dbf_2017$AGB_total[1]
  }
  
change_data = do.call(rbind, changelist)
new_forest_data = do.call(rbind, new_growth)

all_data_2016 = do.call(rbind, list_dbf_2016)
all_sums_2016 = colSums(all_data_2016[2:6], na.rm = TRUE)
all_data_2017 = do.call(rbind, list_dbf_2017)
all_sums_2017=colSums(all_data_2017[2:6], na.rm = TRUE)

sw_data = filter(all_data_2017, grepl("SW", name))
se_data = filter(all_data_2017, grepl("SE", name))
nw_data = filter(all_data_2017, grepl("NW", name))
ne_data = filter(all_data_2017, grepl("NE", name))


sw_data_sum = (colSums(sw_data[2:6], na.rm=TRUE))/1000000/2

se_data_sum = (colSums(se_data[2:6], na.rm=TRUE))/1000000/2
nw_data_sum = (colSums(nw_data[2:6], na.rm=TRUE))/1000000/2
ne_data_sum = (colSums(ne_data[2:6], na.rm=TRUE))/1000000/2

total_area = all_sums_2017[2]
sw_area = colSums(sw_data[3])/total_area
se_area = colSums(se_data[3])/total_area
nw_area = colSums(nw_data[3])/total_area
ne_area = colSums(ne_data[3])/total_area


age_sw = aggregate(sw_data$Count, by =list(sw_data$Value), FUN=sum)
age_se = aggregate(se_data$Count, by =list(se_data$Value), FUN=sum)
age_nw = aggregate(nw_data$Count, by =list(nw_data$Value), FUN=sum)
age_ne = aggregate(ne_data$Count, by =list(ne_data$Value), FUN=sum)

by_regions_age = data.frame("age"=seq(1,32), "SW"= c(age_sw$x) , "SE"= c(age_se$x), "NW"= c(age_nw$x), "NE"= c(age_ne$x))
#write.csv(by_regions_age, "number_pixels_by_regions.csv")

by_regions = data.frame("region"=c("SW", "SE", "NW", "NE", "SW", "SE", "NW", "NE"), "AGB"= c(sw_data_sum[5],se_data_sum[5],nw_data_sum[5],ne_data_sum[5], (sw_data_sum[3]-sw_data_sum[5]) , (se_data_sum[3]-se_data_sum[5]),(nw_data_sum[3]-nw_data_sum[5]), (ne_data_sum[3]-ne_data_sum[5])),
                        "type"=c("Actual", "Actual","Actual","Actual","Potential gain","Potential gain","Potential gain","Potential gain"))

sum_data_2016 = aggregate(all_data_2016$AGB_total, by=list(all_data_2016$Value), FUN=sum)
sum_data_2017 = aggregate(all_data_2017$AGB_total, by=list(all_data_2017$Value), FUN=sum)

sum_data_2016_count = aggregate(all_data_2016$Count, by=list(all_data_2016$Value), FUN=sum)
sum_data_2017_count = aggregate(all_data_2017$Count, by=list(all_data_2017$Value), FUN=sum)


loss_deforestation = (sum(change_data$AGB.deforestation.loss))
net_change_2016_2017 = (sum(change_data$net.change))
new_forest_gain = (sum_data_2017$x[1])
existing_forest_growth = net_change_2016_2017 - loss_deforestation

total_change = loss_deforestation + existing_forest_growth + new_forest_gain 


predict_confidence = predict_future_confidence_interval(all.data_2017)
predict_conf_lower = do.call(rbind, predict_confidence$lower)
predict_conf_upper = do.call(rbind, predict_confidence$upper)
#convert all lower band values below 0 to 0 
#predict_conf_lower[predict_conf_lower <0] =0
predict_lower_age=aggregate(predict_conf_lower[,2:15], by=list(predict_conf_lower$age), FUN=sum)
predict_upper_age=aggregate(predict_conf_upper[,2:15], by=list(predict_conf_upper$age), FUN=sum)


predict_2017_2030 = predict_future(all.data_2017)
future_bind = do.call(rbind, predict_2017_2030)

predict_by_age=aggregate(future_bind[,2:15], by=list(future_bind$age), FUN=sum)

preserved_all = (all_sums_2017[5]/1000000/2) + cumsum(diff((colSums(predict_by_age[,2:15])/1000000/2)))
preserved_5plus = (all_sums_2017[5]/1000000/2) + cumsum(diff((colSums(predict_by_age[5:32,2:15])/1000000/2)))
preserved_10plus = (all_sums_2017[5]/1000000/2) + cumsum(diff((colSums(predict_by_age[10:32,2:15])/1000000/2)))
preserved_15plus = (all_sums_2017[5]/1000000/2) + cumsum(diff((colSums(predict_by_age[15:32,2:15])/1000000/2)))
preserved_20plus = (all_sums_2017[5]/1000000/2) + cumsum(diff((colSums(predict_by_age[20:32,2:15])/1000000/2)))

preserved_all = append(all_sums_2017[5]/1000000/2, preserved_all[1:13]) 
preserved_5plus = append(all_sums_2017[5]/1000000/2, preserved_5plus[1:13]) 
preserved_10plus = append(all_sums_2017[5]/1000000/2, preserved_10plus[1:13]) 
preserved_15plus = append(all_sums_2017[5]/1000000/2, preserved_15plus[1:13]) 
preserved_20plus = append(all_sums_2017[5]/1000000/2, preserved_20plus[1:13]) 


#do same for lower values
preserved_all_lower = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1/1000000/2))) +  (colSums(predict_lower_age[,2:15])/1000000/2)
preserved_5plus_lower = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[5:32]/1000000/2))) +  (colSums(predict_lower_age[5:32,2:15])/1000000/2)
preserved_10plus_lower = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[10:32]/1000000/2))) +  (colSums(predict_lower_age[10:32,2:15])/1000000/2)
preserved_15plus_lower =((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[15:32]/1000000/2))) +  (colSums(predict_lower_age[15:32,2:15])/1000000/2)
preserved_20plus_lower = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[20:32]/1000000/2))) +  (colSums(predict_lower_age[20:32,2:15])/1000000/2)

#preserved_all_lower = append((all_sums_2017[5]/1000000/2), preserved_all_lower[1:13])
#preserved_5plus_lower = append((all_sums_2017[5]/1000000/2), preserved_5plus_lower[1:13])
#preserved_10plus_lower = append((all_sums_2017[5]/1000000/2), preserved_10plus_lower[1:13])
#preserved_15plus_lower = append((all_sums_2017[5]/1000000/2), preserved_15plus_lower[1:13])
#preserved_20plus_lower = append((all_sums_2017[5]/1000000/2), preserved_20plus_lower[1:13])

#do same for upper values
preserved_all_upper =  ((sum(predict_by_age$X1/1000000/2))) - ((all_sums_2017[5]/1000000/2)) +  (colSums(predict_upper_age[,2:15])/1000000/2)
preserved_5plus_upper = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[5:32]/1000000/2))) +  (colSums(predict_upper_age[5:32,2:15])/1000000/2)
preserved_10plus_upper = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[10:32]/1000000/2))) +  (colSums(predict_upper_age[10:32,2:15])/1000000/2)
preserved_15plus_upper =((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[15:32]/1000000/2))) +  (colSums(predict_upper_age[15:32,2:15])/1000000/2)
preserved_20plus_upper = ((all_sums_2017[5]/1000000/2)) - ((sum(predict_by_age$X1[20:32]/1000000/2))) +  (colSums(predict_upper_age[20:32,2:15])/1000000/2)

#preserved_all_upper = append((all_sums_2017[5]/1000000/2), preserved_all_upper[1:13])
#preserved_5plus_upper = append((all_sums_2017[5]/1000000/2), preserved_5plus_upper[1:13])
#preserved_10plus_upper = append((all_sums_2017[5]/1000000/2), preserved_10plus_upper[1:13])
#preserved_15plus_upper = append((all_sums_2017[5]/1000000/2), preserved_15plus_upper[1:13])
#preserved_20plus_upper = append((all_sums_2017[5]/1000000/2), preserved_20plus_upper[1:13])


preserved_all = data.frame("Year"=seq(2017, 2030, 1), "values" = preserved_all)
preserved_all$lower =  preserved_all_lower
preserved_all$upper= preserved_all_upper
preserved_all$type = "A"

plot(preserved_all$values, ylim=c(240, 600))
par(new=T)
plot(preserved_all$lower, ylim=c(240, 600))
par(new=T)
plot(preserved_all$upper, ylim=c(240, 600))

preserved_5plus = data.frame("Year"=seq(2017, 2030, 1), "values" = preserved_5plus)
preserved_5plus$lower =  preserved_5plus_lower
preserved_5plus$upper= preserved_5plus_upper
preserved_5plus$type = "B"

preserved_10plus = data.frame("Year"=seq(2017, 2030, 1), "values" = preserved_10plus)
preserved_10plus$lower =  preserved_10plus_lower
preserved_10plus$upper= preserved_10plus_upper
preserved_10plus$type = "C"

preserved_15plus = data.frame("Year"=seq(2017, 2030, 1), "values" = preserved_15plus)
preserved_15plus$lower =  preserved_15plus_lower
preserved_15plus$upper= preserved_15plus_upper
preserved_15plus$type = "D"


preserved_20plus = data.frame("Year"=seq(2017, 2030, 1), "values" = preserved_20plus)
preserved_20plus$lower =  preserved_20plus_lower
preserved_20plus$upper= preserved_20plus_upper
preserved_20plus$type = "E"

future_cummulative = rbind(preserved_all, preserved_5plus,preserved_10plus,preserved_15plus, preserved_20plus)

preserved_all2 = (diff((colSums(predict_by_age[,2:15])/1000000/2)))
preserved_5plus2 = (diff((colSums(predict_by_age[5:32,2:15])/1000000/2)))
preserved_10plus2 = (diff((colSums(predict_by_age[10:32,2:15])/1000000/2)))
preserved_15plus2 = (diff((colSums(predict_by_age[15:32,2:15])/1000000/2)))
preserved_20plus2 = (diff((colSums(predict_by_age[20:32,2:15])/1000000/2)))

sum_data_2017_area = aggregate(all_data_2017$area_hect, by=list(all_data_2017$Value), FUN=sum)

area_preserved_all = sum(sum_data_2017_area$x)/1000000/2
area_preserved_5plus = sum(sum_data_2017_area$x[5:32])/1000000/2
area_preserved_10plus = sum(sum_data_2017_area$x[10:32])/1000000/2
area_preserved_15plus = sum(sum_data_2017_area$x[15:32])/1000000/2
area_preserved_20plus = sum(sum_data_2017_area$x[20:32])/1000000/2


preserved_all2 = data.frame("Year"=seq(2018, 2030, 1),"values"=(diff((colSums(predict_by_age[,2:15])/1000000/2))))
#preserved_all2$lower
#confidence intervals expressed as decimals
lower_perc= ((colSums(predict_by_age[,2:15])/1000000/2)-((colSums(predict_lower_age[,2:15])/1000000/2)))/(colSums(predict_by_age[,2:15])/1000000/2)
upper_perc= (((colSums(predict_upper_age[,2:15])/1000000/2))-(colSums(predict_by_age[,2:15])/1000000/2))/(colSums(predict_by_age[,2:15])/1000000/2)

preserved_all2$lower = preserved_all2$values * (1-lower_perc[2:14])
preserved_all2$upper = preserved_all2$values * (1 + upper_perc[2:14])
preserved_all2$type = "A"


preserved_5plus2 = data.frame("Year"=seq(2018, 2030, 1),"values"=(diff((colSums(predict_by_age[5:32,2:15])/1000000/2))))

lower_perc= ((colSums(predict_by_age[5:32,2:15])/1000000/2)-((colSums(predict_lower_age[5:32,2:15])/1000000/2)))/(colSums(predict_by_age[5:32,2:15])/1000000/2)
upper_perc= (((colSums(predict_upper_age[5:32,2:15])/1000000/2))-(colSums(predict_by_age[5:32,2:15])/1000000/2))/(colSums(predict_by_age[5:32,2:15])/1000000/2)
preserved_5plus2$lower = preserved_5plus2$values * (1-lower_perc[2:14])
preserved_5plus2$upper = preserved_5plus2$values * (1 + upper_perc[2:14])
preserved_5plus2$type = "B"

preserved_10plus2 = data.frame("Year"=seq(2018, 2030, 1),"values"=(diff((colSums(predict_by_age[10:32,2:15])/1000000/2))))

lower_perc= ((colSums(predict_by_age[10:32,2:15])/1000000/2)-((colSums(predict_lower_age[10:32,2:15])/1000000/2)))/(colSums(predict_by_age[10:32,2:15])/1000000/2)
upper_perc= (((colSums(predict_upper_age[10:32,2:15])/1000000/2))-(colSums(predict_by_age[10:32,2:15])/1000000/2))/(colSums(predict_by_age[10:32,2:15])/1000000/2)
preserved_10plus2$lower = preserved_10plus2$values * (1-lower_perc[2:14])
preserved_10plus2$upper = preserved_10plus2$values * (1 + upper_perc[2:14])
preserved_10plus2$type = "C"

preserved_15plus2 = data.frame("Year"=seq(2018, 2030, 1),"values"=(diff((colSums(predict_by_age[15:32,2:15])/1000000/2))))
lower_perc= ((colSums(predict_by_age[15:32,2:15])/1000000/2)-((colSums(predict_lower_age[15:32,2:15])/1000000/2)))/(colSums(predict_by_age[15:32,2:15])/1000000/2)
upper_perc= (((colSums(predict_upper_age[15:32,2:15])/1000000/2))-(colSums(predict_by_age[15:32,2:15])/1000000/2))/(colSums(predict_by_age[15:32,2:15])/1000000/2)
preserved_15plus2$lower = preserved_15plus2$values * (1-lower_perc[2:14])
preserved_15plus2$upper = preserved_15plus2$values * (1 + upper_perc[2:14])
preserved_15plus2$type = "D"

preserved_20plus2 = data.frame("Year"=seq(2018, 2030, 1),"values"=(diff((colSums(predict_by_age[20:32,2:15])/1000000/2))))
lower_perc= ((colSums(predict_by_age[20:32,2:15])/1000000/2)-((colSums(predict_lower_age[20:32,2:15])/1000000/2)))/(colSums(predict_by_age[20:32,2:15])/1000000/2)
upper_perc= (((colSums(predict_upper_age[20:32,2:15])/1000000/2))-(colSums(predict_by_age[20:32,2:15])/1000000/2))/(colSums(predict_by_age[20:32,2:15])/1000000/2)
preserved_20plus2$lower = preserved_20plus2$values * (1-lower_perc[2:14])
preserved_20plus2$upper = preserved_20plus2$values * (1 + upper_perc[2:14])
preserved_20plus2$type = "E"

future_nets = rbind(preserved_all2, preserved_5plus2, preserved_10plus2, preserved_15plus2, preserved_20plus2)

#future_nets = data.frame("Year"=seq(2018, 2030, 1), "All preserved"= preserved_all2, "5+ forests preserved"= preserved_5plus2, "10+ forests preserved"= preserved_10plus2,
 #                        "15+ forests preserved"= preserved_15plus2, "20+ forests preserved"= preserved_20plus2)


#all_data_2017$type = rownames(all_data_2017)
data_deforest = dplyr::filter(all_data_2017, grepl("deforest|both_dist", name))

#ggplot(sum_data, x = "Group.1", y = "")

fire_impact=dplyr::filter(all_data_2017, grepl("fire", name))
fire_dif = sum(fire_impact$no_dist - fire_impact$AGB_total)

def_impact=dplyr::filter(all_data_2017, grepl("deforest", name))
def_dif=sum(def_impact$no_dist - def_impact$AGB_total)

both_impact=dplyr::filter(all_data_2017, grepl("both_dist", name))
both_dif = sum(both_impact$no_dist - both_impact$AGB_total)



disturbance_only_2017 = data.frame(sinks = factor(c( "Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances")),
                              levels = c( "Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances"), 
                              AGB = c( fire_dif, def_dif, both_dif),
                              member = c("Potential sink","Potential sink","Potential sink"))

summary_data2_2016 = data.frame(sinks = factor(c("Wang et al.", "Fearnside and \nGuimarães (1996)","sink 2017","No \ndisturbance", "Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances")),
                          levels = c("Wang et al.", "Fearnside (1996)","No \ndisturbance", "2017 sink","Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances"), 
                          AGB = c(all_sums_2016[5],all_sums_2016[7],all_sums_2016[5],all_sums_2016[9], fire_dif, def_dif, both_dif),
                          member = c("Other1","Other2","sink 2017","Potential sink","sink 2017","sink 2017","sink 2017"),
                          disturbance = c("Other1","Other2","No", "both", "Yes", "Yes", "Yes"))


summary_net_changes = data.frame(sinks = factor(c("Ongoing growth", "New forest growth", "Deforestation", "Net_sink")),
                                 AGB = c(existing_forest_growth, new_forest_gain, loss_deforestation, total_change),
                                 member = c("Loss.gains", "Loss.gains", "Loss.gains", "net"),
                                 disturbance = c("Existing forest growth", "New forest growth","Deforestation", "Net change" ))

summary_data_2017_only = data.frame(sinks = factor(c("sink 2017" ,"Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances", "dummy")),
                               levels = c("2017 sink","Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances", "dummy"), 
                               AGB = c(all_sums_2017[5], fire_dif, def_dif, both_dif, 0),
                               member = c("sink 2017","sink 2017","sink 2017","sink 2017", "zummy"),
                               disturbance = c("No", "Yes", "Yes", "Yes", "Yes"))

summary_data_2017 = data.frame(sinks = factor(c("sink 2017","No \ndisturbance", "Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances", "netchange")),
                               levels = c( "2017 sink","No \ndisturbance","Fire \ndisturbance", "Deforestation \ndisturbance", "Both  \ndisturbances", "netchange"), 
                               AGB = c((all_sums_2017[5]-total_change),all_sums_2017[5], fire_dif, def_dif, both_dif, total_change),
                               member = c("A","B","B","B","B", "A"),
                               disturbance = c("No", "No", "Yes", "Yes", "Yes", "Yes"))


summary_data_2017$AGB= summary_data_2017$AGB/1000000/2

summary_data_2017_only$AGB= summary_data_2017_only$AGB/1000000/2


summary_net_changes$AGB = summary_net_changes$AGB/1000000/2

disturbance_only_2017$AGB= disturbance_only_2017$AGB/1000000/2


# plotting ----------------------------------------------------------------

p.main_new = ggplot(summary_data_2017, aes(x = member, y = AGB, fill=factor(disturbance, levels = c("Yes", "No")))) + 
  geom_bar(stat = "identity", color = "white") +
  coord_flip()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=12),axis.text.y = element_text(color="black", size=12), 
        text=element_text(size=12)) + 
  scale_fill_manual(values = c("#b8cafc",  "royalblue4" ,"gray50", "gray50")) + 
  scale_x_discrete(labels = c( "2017\nstock", "Potential\n stock", ""), name ="") + 
  scale_y_continuous(name = expression('AGC (Tg C)')) + 
  theme(legend.position="none") + 
  annotate("text", x=2, y = (summary_data_2017_only$AGB[1])/2, label=round(summary_data_2017_only$AGB[1], digits=1),size=4.7, color="white") + 
  annotate("text", x=1, y = (summary_data_2017_only$AGB[1])/2, label=round(all_sums_2016[5]/1000000/2, digits=1),size=4.7, color="white") + 
  annotate("text", x=1, y = 310, label="11.9",size=4.7, color="black")


p.main_new

ggsave("AGC_sink_aug2020v2_as_stock_long.pdf", p.main_new,width =5, height=2, units="in" , dpi=1200)

p.region=ggplot(by_regions, aes(x=region, y=AGB, fill=factor(type, levels = c("Potential gain", "Actual")))) +
  geom_bar(stat = "identity", color = "white") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=12),axis.text.y = element_text(color="black", size=12), 
        text=element_text(size=12)) + 
  scale_fill_manual(values = c("#b8cafc",  "royalblue4")) +
  scale_x_discrete(name="Regions") +
  scale_y_continuous(name = expression('AGC (Tg C)')) +
  geom_text(size=4.7, position = position_stack(vjust=0.5), label=round(by_regions$AGB, digits=1))+
  labs(fill="Carbon stocks")

p.region

ggsave("AGC_sink_by_region_aug2020.jpeg", p.region,width =5, height=4, units="in" , dpi=500)


p.zoom=ggplot(disturbance_only_2017, aes(x = member, y = AGB, fill=sinks)) + 
  geom_bar(stat = "identity", color = "white") +
  #coord_flip()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14), 
        text=element_text(size=14)) + 
  #scale_fill_manual(values = c("#DC267F","#648FFF", "#40B0A6")) + 
  scale_fill_manual(values = c("#DC267F","orangered2", "darkgoldenrod2")) + 
  scale_x_discrete(labels = c("Reduction to \npotential stock"), name ="") + 
  scale_y_continuous(name = expression('AGC (Tg C)')) + 
  theme(legend.title=element_blank(), legend.text = element_text(size=14))+
  geom_text(aes(label=round(AGB, digits=1)),position = position_stack(vjust = 0.5), size=4.7, color="white")+
  theme(legend.position="none") 

p.zoom
#ggsave("AGC_sink_reduction_2017_stock_sept2020.jpeg", p.zoom,width =2, height=3.5, units="in" , dpi=500)
ggsave("AGC_sink_reduction_2017_stock_sept2020.pdf", p.zoom,width =2, height=3.5, units="in" , dpi=500)

p.net = ggplot(summary_net_changes, aes(x = member, y=AGB, fill = factor(disturbance, levels = c("Existing forest growth", "New forest growth","Deforestation", "Net change")))) +
  geom_bar(stat="identity", color = "white") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=12),axis.text.y = element_text(color="black", size=12), 
        text=element_text(size=12)) + 
  scale_fill_manual(values = c("royalblue4", "royalblue1", "#FF5A45", "darkorchid3")) + 
  scale_x_discrete(labels = c("2017 \nGains & Losses", "2017 net sink"), name ="Sink fluxes") + 
  scale_y_continuous(name = expression('AGC (Tg C)')) + 
  geom_hline(yintercept=0, color="grey40", linetype = "dashed")+
  guides(fill=guide_legend(title="Fluxes"))+ 
  geom_text(aes(label=round(AGB, digits=1)),position = position_stack(vjust = 0.5), size=4.7, color="white")
  
p.net

ggsave("AGC_fluxes_sept_2020.jpeg", p.net,width =6, height=4, units="in" , dpi=500)


p.region_age = ggplot(by_regions_age, aes(x=age)) +
  geom_line(aes(y = SW, color = "SW&C"), size=2) +
  geom_line(aes(y = SE, color = "SE&N"), size=2) +
  geom_line(aes(y = NW, color = "NW"), size=2) +
  geom_line(aes(y = NE, color = "NE&NC"), size=2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14), 
        text=element_text(size=14)) + 
  scale_x_continuous(name="Age (Years)") +
  scale_y_continuous(name = "Number of pixels") +
  scale_color_manual(values = c('SW&C' = '#06f001', 'SE&N' = '#1704a0', 'NW'= '#059a01', 'NE&NC'='#05a9ff')) +
  labs(color = 'Regions') +
  theme(legend.position = c(0.8, 0.7), legend.title=element_text(size=14), legend.text=element_text(size=14))

p.region_age

ggsave("Age_by_region.jpeg", p.region_age,width =5, height=4, units="in" , dpi=500)

#melted_future = melt(future_nets, id="Year")

annual_acc = ggplot(data=future_nets, aes(x=Year, y =values, color=type, ymax=upper, ymin=lower)) +
  geom_line(size=1.1)+ 
  geom_ribbon(alpha=0.3, aes(fill=type), color=NA)+
  #scale_color_manual(values=c("#22118E", "#2B98C7", "#E8D40F", "#FE7C00", "firebrick"))+
  scale_color_manual(values=c("steelblue", "darkseagreen3", "darkorange2", "darkorchid2", "grey50"), 
                              labels =c("All ages preserved (13.8Mha)", "Age 5+ preserved (9.5Mha)", "Age 10+ preserved (5.6Mha)", "Age 15+ preserved (3.5Mha)", "Age 20+ preserved (2.2Mha)"))+
  scale_fill_manual(values=c("steelblue", "darkseagreen3", "darkorange2", "darkorchid2", "grey50"))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14), 
        text=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_blank())+
  scale_y_continuous(name =expression('Annual C accumulation (Tg C yr'^-1*')')) +
  #theme(legend.position = c(0.86, 0.80)) +
  #annotate(geom="point", x = 2025, y = 10.9, color = "red", size= 3) + 
  #annotate(geom="point", x = 2030, y = 9.81, color = "red", size= 3) +
  scale_x_continuous(breaks= c(2020, 2025, 2030)) + 
  theme(legend.position = "none")

annual_acc

ggsave("AGC_annual_future_with_uncertainty_aug2020.jpeg", annual_acc,width =6, height=4, units="in" , dpi=500)

getwd()


annual_increase = ggplot(data=future_cummulative, aes(x=Year, y =values)) +
  geom_line(size=1.1, aes(color=type))+ 
  geom_ribbon(alpha=0.3, linetype=0, aes(ymax=upper, ymin=lower,fill=type), show.legend=FALSE)+
  geom_ribbon(alpha=0.3, linetype=0, aes(ymax=values, ymin=values,fill=type))+
  scale_color_manual(values=c("steelblue", "darkseagreen3", "darkorange2", "darkorchid2", "grey50"), 
                     labels =c("All ages preserved (13.8Mha)", "Age 5+ preserved (9.5Mha)", "Age 10+ preserved (5.6Mha)", "Age 15+ preserved (3.5Mha)", "Age 20+ preserved (2.2Mha)"))+
  scale_fill_manual(values=c("steelblue", "darkseagreen3", "darkorange2", "darkorchid2", "grey50"),
                    labels =c("All ages preserved (13.8Mha)", "Age 5+ preserved (9.5Mha)", "Age 10+ preserved (5.6Mha)", "Age 15+ preserved (3.5Mha)", "Age 20+ preserved (2.2Mha)"))+#)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14), 
        text=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.background = element_blank())+
  scale_y_continuous(name =expression('Preserved stock (Tg C)')) +
  theme(legend.position = c(0.25, 0.83)) +
  scale_x_continuous(breaks= c(2020, 2025, 2030))

annual_increase

ggsave("AGC_sink_accumulation_with_uncertainty_aug2020.jpeg", annual_increase,width =6, height=4, units="in" , dpi=500)

future_figure = ggarrange(annual_increase,annual_acc,ncol=2, nrow=1)
future_figure
ggsave("AGC_future_Fig5.pdf", future_figure, width = 12, height=4, units = "in")

getwd()
