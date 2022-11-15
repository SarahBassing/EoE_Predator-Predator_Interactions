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
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(-X)
  # cams_eoe <- read.csv("./Data/IDFG camera data/cams_eoe_skinny.csv") %>%
  #   dplyr::select(-X)
  # cams_wolf <- read.csv("./Data/IDFG camera data/cams_wolf_skinny.csv") %>%
  #   dplyr::select(-X)
  
  #' #'  Full detection data sets
  #' #'  Motion triggered images
  #' load("./Data/IDFG camera data/Split datasets/eoe20s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe20w_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe21s_allM.RData")
  #' 
  #' load("./Data/IDFG camera data/Split datasets/wolf19s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf20s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf21s_allM.RData")
  #' 
  #' #'  Time triggered images
  #' load("./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  #' 
  #' load("./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  #'  Keeper data sets
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  #'  Add camera setup info (P/U/O/A) and consistent naming structure to match 
  #'  camera location data
  add_deploy_info <- function(dets, season, pred) { 
    sub_cams <- cams_eoe_long %>%
      filter(Season == season) %>%
      filter(Setup == pred)
    pcams <- dets$CamID %in% sub_cams$CamID
    dets$pcams <- pcams
    dets <- dets %>%
      mutate(Setup = ifelse(pcams == "TRUE", "P", "U"),
             NewLocationID = paste0("GMU", Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-pcams)
    return(dets)
  }
  #'  EoE cameras
  eoe_seasons <- list("Smr20", "Wtr20", "Smr21")
  eoe_motion_list <- mapply(add_deploy_info, season = eoe_seasons, dets = eoe_motion_skinny, pred = "predator", SIMPLIFY = FALSE) 
  eoe_time_list <- mapply(add_deploy_info, season = eoe_seasons, dets = eoe_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  #'  Wolf cameras
  wolf_seasons <- list("Smr20", "Wtr20", "Smr21")
  wolf_motion_list <- mapply(add_deploy_info, season = wolf_seasons, dets = wolf_motion_skinny, pred = "predator", SIMPLIFY = FALSE)  
  wolf_time_list <- mapply(add_deploy_info, season = wolf_seasons, dets = wolf_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  
  
  
    
  eoe20s_noon <- eoe_time_list[[1]]
  
  