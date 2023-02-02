  #'  ----------------------------------
  #'  Spatial variation in detection
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  January 2023
  #'  ----------------------------------
  #'  Explore extent to which detections correlate with camera placement
  #'  ----------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(tidyverse)
  library(ggplot2)

  #'  Load detection histories
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20w_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")

  #'  Load camera station, problem cams, and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Target", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) #' Remove duplicate camera where height changed slightly
  
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  format_cam_station <- function(cams, season) {
    cams <- dplyr::select(cams, c("NewLocationID", "Lat", "Long")) %>%
      mutate(GMU = str_extract(NewLocationID, "[^_]+"),
             Season = season) %>%
      left_join(cams_eoe_long[cams_eoe_long$Season == season,], by = c("NewLocationID", "Season")) %>%
      #'  Drop handful of duplicated rows (not sure why this happens with left_join)
      unique() %>%
      #'  Reduce categories to random, road, trail
      mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | 
                                     CameraFacing == "decommision" | CameraFacing == "decommission", "road", CameraFacing),
             CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing),
             CameraFacing = ifelse(is.na(CameraFacing), "trail", CameraFacing)) %>% # predator cameras so either trail or road
      arrange(NewLocationID)
    return(cams)
  }
  cams_eoe20s <- format_cam_station(eoe_probcams_20s, season = "Smr20")
  cams_eoe20w <- format_cam_station(eoe_probcams_20w, season = "Wtr20")
  cams_eoe21s <- format_cam_station(eoe_probcams_21s, season = "Smr21") 
  
  # load("./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  # load("./Data/Covariates_extracted/Covariates_EoE_Wtr20.RData")
  # load("./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  
  
  #'  ----------------------------------------------------
  ####  Summary of detection events by camera deployment  ####
  #'  ----------------------------------------------------
  #'  Function to tally detections across primary period based on camera deployment
  #'  and plot distribution of detection/non-detections based on deployment type/location
  detection_by_placement <- function(dets, cams, sppseason) {
    #'  Sum total detections per camera and indicate if 1 or more detection occurred
    dets <- as.data.frame(dets) %>%
      rownames_to_column(var = "NewLocationID") %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(ndets = rowSums(.[2:11], na.rm = T),
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels =c("random", "road", "trail"))) %>%
      rowwise() %>%
      mutate(det = max(c_across(2:11), na.rm = T)) %>%
      dplyr::select(c(Setup, CameraFacing, det, ndets)) %>%
      filter(det != "-Inf")
    
    #'  Table data by camera set up and trail type
    print(sppseason)
    det_by_setup <- dplyr::select(dets, c(Setup, det))
    print(table(det_by_setup))
    print(ggplot(dets, aes(x = det, fill = Setup)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      labs(fill="") +
      ggtitle(paste(sppseason, "detections by camera deployment")))
    
    det_bytrailtype <- dplyr::select(dets, c(CameraFacing, det))
    print(table(det_bytrailtype))
    print(ggplot(dets, aes(x = det, fill = CameraFacing)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
      labs(fill="") +
      ggtitle(paste(sppseason, "detections by camera deployment")))
    
    ndet_by_setup <- dplyr::select(dets, c(Setup, ndets))
    print(table(ndet_by_setup))
    print(ggplot(dets, aes(x = ndets, fill = Setup)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      labs(fill="") +
      ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    ndet_bytrailtype <- dplyr::select(dets, c(CameraFacing, ndets))
    print(table(ndet_bytrailtype))
    print(ggplot(dets, aes(x = ndets, fill = CameraFacing)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
      labs(fill="") +
      ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    return(dets)
  }
  wolf20s <- detection_by_placement(DH_eoe20s_predators[[5]][[1]], cams = cams_eoe20s, sppseason = "Wolf Smr20")
  wolf21s <- detection_by_placement(DH_eoe21s_predators[[5]][[1]], cams = cams_eoe21s, sppseason = "Wolf Smr21")
  coy20s <- detection_by_placement(DH_eoe20s_predators[[3]][[1]], cams = cams_eoe20s, sppseason = "Coyote Smr20")
  coy20s <- detection_by_placement(DH_eoe21s_predators[[3]][[1]], cams = cams_eoe21s, sppseason = "Coyote Smr21")
  lion20s <- detection_by_placement(DH_eoe20s_predators[[4]][[1]], cams = cams_eoe20s, sppseason = "Lion Smr20")
  lion20s <- detection_by_placement(DH_eoe21s_predators[[4]][[1]], cams = cams_eoe21s, sppseason = "Lion Smr21")
  bear20s <- detection_by_placement(DH_eoe20s_predators[[1]][[1]], cams = cams_eoe20s, sppseason = "Bear Smr20")
  bear20s <- detection_by_placement(DH_eoe21s_predators[[1]][[1]], cams = cams_eoe21s, sppseason = "Bear Smr21")
  bob20s <- detection_by_placement(DH_eoe20s_predators[[2]][[1]], cams = cams_eoe20s, sppseason = "Bobcat Smr20")
  bob20s <- detection_by_placement(DH_eoe21s_predators[[2]][[1]], cams = cams_eoe21s, sppseason = "Bobcat Smr21")
  
  
  #### table human detections, then make plots  ###
  

  
  
  
  
  
  