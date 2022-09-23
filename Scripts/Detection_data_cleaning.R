  #'  ---------------------------------
  #'  Cleaning raw detection data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean camea deteciton data, make sure annual datasets are consistent, 
  #'  filter into manageable data sets, and connect to camera locations.
  #'  ---------------------------------
  
  #'  Load libraries
  library(sf)
  library(raster)
  library(ggplot2)
  library(tidyverse)
  
  #'  Camera locations
  cams_eoe <- read.csv("./Data/IDFG camera data/cams_eoe_skinny.csv") %>%
    dplyr::select(-X)
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(-X)
  # cams_wolf <- read.csv("./Data/IDFG camera data/cams_wolf_skinny.csv") %>%
  #   dplyr::select(-X)
  
  #'  Detection data
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM.RData")
  
  
  
  
  
  
  
  
  