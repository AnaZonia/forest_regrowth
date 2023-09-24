
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4479398.svg)](https://doi.org/10.5281/zenodo.4479398)

# Secondary Forest regrowth in Amazonia, Brazil
This repository contains the code used to produce the figures in the main paper of Heinrich et al.

## Scope
The code in this repositroy can be used to reproduce the main Figures in Heinrich et al. The maps of secondary forest extent and ages were developed based on the Land-Use Land-cover dataset Mapbiomas 3.1. This step is described in much greater detail in Silva et al., 2020 and so is not shown here. The code and maps of secondary forest age and extent can be downloaded here: https://github.com/celsohlsj/gee_brazil_sv. \
In this repository the scripts used to make the main figures in the paper are shown. Please note that the code has not been amended for wider use and still contains set working directories for use with University of Bristol systems, you will need to change these for the scripts to run. This repository contains:

1. "Fig1a_f_plot.R" - this file will read .csv files containing infomation on the median-bias corrected Aboveground Biomass (AGB) values as a function of age, disaggregated by limits of a given driving variable (this is the raw data seen as dots in the figure). The script will then read in the regrowth models under different driving conditions expressed in terms of Aboveground Carbon (AGC) produced in this reseach contained in RData format. Set working directory - at or near lines 8 a(ctrl +F: setwd to check). The following packages will need to be installed (and their dependencies):
  * nlstools
  * investr
  * qpcr
  * pgirmess

2. "Fig3_plot.R" - similar to (1), this script will read in .csv files containing the infomation of the median-bias corrected AGB values as a function of age. This time the data is disaggregated by the 4 regions identified according to similarities in the climatological variables (Short-wave radiation, precipitation, and the Maximum Cumulative Water Deficit index - MCWD). In each region, the regrowth models are disaggregated by disturbance (fires and/or repeated deforestations, or neither). The same packages as (1) are required. Ultimately, 16 models are read in - from 4 regions, with 4 submodels. Set working directory - at or near lines 8 (ctrl +F: setwd to check). 

3. "Fig1g_2b_e_plot.R" - this script will plot the average ranked importance of each driving variable for the whole Amazon and the regions shown in Figure 1g and 2b-e of the main paper. Set working directory - at or near line 136 (ctrl +F: setwd to check). The following packages will need to be installed (and their dependencies):
  * caret
  * caTools
  * randomForest
  * party
  * foreign
  * tidyverse
  * matrixStats
  * ggplot2
  * ggpubr
  * egg

4. "Fig4_Fig5_plot.R" - this script uses the .dbf file of the .tiff files for each secondary forest pixel in each region shown in Figure 2a, and further disaggregated by the kind of disturbance experienced. Note, this script uses another script as source code, this is found in the same directory as the data needed to run this script. Set working directory - at or near lines 21 - 25 (ctrl +F: setwd to check). The following packages will need to be installed (and their dependencies):
  * raster
  * ggpubr
  * foreign
  * ggplot2
  * sf
  * rgdal
  * dplyr

## Usage
The data produced in this project were produced using a combination of programming languages due to differences in the author's preferences and expertise. Overall the initial data analysis was carried out in Python2.7 (Arcpy) and Python.3.7. Most of the post-processing of the inital data was then carried out in R (v3.6) - code for which is shown here.  

## Dataset access and credit
The final data can be accessed from the Zenodo repository:https://zenodo.org/deposit/4290633 

All code shown here is based on the work by Heinrich et al. and forms the basis my (Viola Heinrich) PhD research. If you make use of any of the code, methodology or data, please cite the main paper Heinrich et al. 2021 and any associated Data and Code DOIs. Thank you.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.


