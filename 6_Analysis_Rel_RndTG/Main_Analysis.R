## Last Modification 6/6/2021
## SYS-Mut Real-TG Running the model

graphics.off() 
rm(list=ls())
start_time <- Sys.time()

#--------------- Setup parallel back end to use many processors
library(iterators)
library(parallel)
library(foreach)
library(doParallel)

#----------------
library(stringi)
library(stringr)
library(RTCGAToolbox)
library(lattice)
library(gridExtra)
library(gapminder)
library(grid)
library(fpc)
library(DescTools)
library(data.table)
library(dplyr)

################# Library
library(fitdistrplus)
library(survival)
library(npsurv)
library(lsei)
library(readxl)
library(writexl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(DescTools)
library(MASS)


#---------------- Set-path
CurrentDir=getwd()
setwd(CurrentDir)


source("../1_PreProcess/Init_Parameters.R")
InitParams<<-Init_Parameters()
CANCER_TYPE<<-InitParams$CANCER

#--------------- Parameters
source("Parameters.R")
Params<<-Parameters()

#--------------- Setup parallel back-end to use many processors

source("Average_RnG_WI_NewFilter_logRank.R")
Average_RnG_WI_NewFilter_logRank()

source("Result_Average_RnG_WI_NewFilter.R")
Result_Average_RnG_WI_NewFilter()

source("Charts_Rnd_Real_WI_NewFilter.R")
Charts_Rnd_Real_WI_NewFilter()

source("Finding_Outliers_Using_Lognorm_WI_NewFilter.R")
Finding_Outliers_Using_Lognorm_WI_NewFilter()

source("Finding_Outliers_Using_Lognorm_WI_NewFilter_BasedonrandomRuns.R")
Finding_Outliers_Using_Lognorm_WI_NewFilter_BasedonrandomRuns()

# source("All_GI_Imp_ECDF_10RUN_EstimVal.R")
# All_GI_Imp_ECDF_10RUN_EstimVal()
# 
# source("All_GI_Imp_ECDF_10RUN_TrueVal.R")
# All_GI_Imp_ECDF_10RUN_TrueVal()

source("All_GI_Imp_ECDF_RelTG_TrueVal.R")
All_GI_Imp_ECDF_RelTG_TrueVal()

# source("All_GI_Imp_ECDF_EstimRelTG_Param.R")
# All_GI_Imp_ECDF_EstimRelTG_Param()