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
  #'  Motion triggered images
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM.RData")
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/wolf20s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/wolf21s_allM.RData")
  
  #'  Time triggered images
  load("./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  #'  Filter down to "Keepers"
  keepers <- function(dets) {
    keep_dets <- filter(dets, SetLocation == "KeepSet")
    return(keep_dets)
  }
  eoe_motion <- list(eoe20s_allM, eoe20w_allM, eoe21s_allM)
  eoe_motion_skinny <- lapply(eoe_motion, keepers)
  # save(eoe_motion_skinny, file = "./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  eoe_time <- list(eoe20s_allT, eoe20w_allT, eoe21s_allT)
  eoe_time_skinny <- lapply(eoe_time, keepers)
  # save(eoe_time_skinny, file = "./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  
  wolf_motion <- list(wolf19s_allM, wolf20s_allM, wolf21s_allM)
  wolf_motion_skinny <- lapply(wolf_motion, keepers)
  # save(wolf_motion_skinny, file = "./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  
  wolf_time <- list(wolf19s_allT, wolf20s_allT, wolf21s_allT)
  wolf_time_skinny <- lapply(wolf_time, keepers)
  # save(wolf_time_skinny, file = "./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  #'  Load in keeper data
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  
  
  
  