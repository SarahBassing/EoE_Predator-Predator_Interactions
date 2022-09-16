  #'  ---------------------------------
  #'  Cleaning raw camera data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean deployment and deteciton data, make sure annual datasets 
  #'  are consistent, filter and merge.
  #'  ---------------------------------
  
  #'  Load libraries
  library(tidyverse)
  
  #'  Camera deployment data
  load("./Data/IDFG camera data/eoe2020s_cameras.RData")
  load("./Data/IDFG camera data/eoe2020w_cameras.RData")
  load("./Data/IDFG camera data/eoe2021s_cameras.RData")
  
  load("./Data/IDFG camera data/wolf2019s_cameras.RData")
  load("./Data/IDFG camera data/wolf2020s_cameras.RData")
  load("./Data/IDFG camera data/wolf2021s_cameras.RData")
  
  #'  Match columns
  names(cams_s20_eoe)
  names(cams_w20_eoe)
  names(cams_s21_eoe)
  names(cams_s19_wolf)
  names(cams_s20_wolf)
  names(cams_s21_wolf)
  
  #'  Which cams were deployed >1 year on EoE project
  same_cams <- cams_s21_eoe$CamID %in% cams_s20_eoe$CamID
  two_yr_cams <- cams_s21_eoe
  two_yr_cams$cams_s20 <- same_cams
  two_yr_cams <- dplyr::select(two_yr_cams, c(Region, Gmu, Setup, LocationID, Lat, Long, CamID, cams_s20))
  