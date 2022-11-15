  #'  ---------------------------------
  #'  Cleaning raw detection data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean camea deteciton data, make sure annual datasets are consistent, 
  #'  filter into manageable data sets, and connect to camera locations.
  #'  ---------------------------------
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(sf)
  library(raster)
  library(ggplot2)
  library(tidyverse)
  
  #'  Camera locations
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(-X)
  cams_wolf_long <- read.csv("./Data/IDFG camera data/cams_wolf_long.csv") %>%
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
  
  #'  Keeper data sets (all animal/human/vehicle images & noon timelapse images)
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  #'  Add camera setup info (P/U/O/A) and consistent naming structure to match 
  #'  camera location data
  #'  EoE cameras
  eoe_deploy_info <- function(dets, season, pred) { 
    #'  Filter to specific season and predator setup cameras
    sub_cams <- cams_eoe_long %>%
      filter(Season == season) %>%
      filter(Setup == pred)
    #'  Identify which observations came from predator setup cameras (TRUE)
    pcams <- dets$CamID %in% sub_cams$CamID
    #'  Add column to detection data
    dets$pcams <- pcams
    #'  Create new camera location name that matches camera deployment data
    dets <- dets %>%
      mutate(Setup = ifelse(pcams == "TRUE", "P", "U"),
             NewLocationID = paste0("GMU", Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-pcams) %>%
      relocate(NewLocationID, .before = File) %>%
      relocate(Setup, .after = NewLocationID)
    return(dets)
  }
  eoe_seasons <- list("Smr20", "Wtr20", "Smr21")
  eoe_motion_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_motion_skinny, pred = "predator", SIMPLIFY = FALSE) 
  eoe_noon_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  #'  Double check it worked
  eoe20s_noon <- eoe_noon_list[[1]]
  
  #'  Wolf cameras
  wolf_deploy_info <- function(dets, season, abund, abund_occu) { 
    #'  Filter to specific season and Abundance cameras
    sub_abund_cams <- cams_wolf_long %>%
      filter(Season == season) %>%
      filter(Target == abund)
    #'  Identify which observations came from Abundance cameras (TRUE)
    acams <- dets$CamID %in% sub_abund_cams$CamID
    #'  Add column to detection data
    dets$acams <- acams
    #'  Filter to specific season and Abund_Occu cameras
    sub_abundoccu_cams <- cams_wolf_long %>%
      filter(Season == season) %>%
      filter(Target == abund_occu)
    #'  Identify which observations came from Abund_Occu cameras (TRUE)
    aocams <- dets$CamID %in% sub_abundoccu_cams$CamID
    #'  Add column to detection data
    dets$aocams <- aocams
    #'  Create new camera location name that matches camera deployment data
    dets <- dets %>%
      mutate(Setup = ifelse(acams == "TRUE", "A", "O"),
             Setup = ifelse(aocams == "TRUE", "B", Setup),
             NewLocationID = paste0("GMU", Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-c(acams, aocams)) %>%
      relocate(NewLocationID, .before = File) %>%
      relocate(Setup, .after = NewLocationID)
    return(dets)
  }
  wolf_seasons <- list("Smr19", "Smr20", "Smr21")
  wolf_motion_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_motion_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)  
  wolf_noon_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_time_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)  
  #'  Double check it worked
  wolf20s_noon <- wolf_noon_list[[2]]
  
  