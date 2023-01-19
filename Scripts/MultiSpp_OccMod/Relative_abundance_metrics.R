  #'  -------------------------------
  #'  Relative abundance index
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2023
  #'  -------------------------------
  #'  Script to generate relative abundance indices for different species detected
  #'  at each camera site and compare them using Pearson's correlation coefficient
  #'  to determine if there are notable differences in how they represent 
  #'  relative abundance at the camera site. 
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(tidyverse)
  # library(timetk)
  
  #'  ----------------
  ####  Read in data  ####
  #'  ----------------
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe20w_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  
  
  #'  ------------------------
  ###  Filter detection data  ####
  #'  ------------------------
  #'  1) Remove sequential problem images from larger image set
  thin_detections <- function(dets, seqprobs) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    return(skinny_dets)
  }
  eoe20s_allM_skinny <- thin_detections(eoe20s_allM, eoe_seqprob_20s)
  eoe20w_allM_skinny <- thin_detections(eoe20w_allM, eoe_seqprob_20w)
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s)
  
  #'  2) Filter detection data to time period and species of interest
  #'  Species of interest are ungulates, humans/vehicles, and domestic animals
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      filter(Species == "elk" | Species == "moose" | Species == "muledeer" | Species == "whitetaileddeer" |
               Species == "rabbit_hare" | Species == "human" | Species == "dog_domestic" | 
               Species == "horse" |  Species == "cattle_cow" | Species == "cat_domestic" |
               Vehicle == "TRUE") %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time),
        Species = ifelse(Vehicle == "TRUE", "human_motorized", Species)
      ) %>%
      dplyr::select(-Vehicle) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      #'  Remove maintenance images so not included in human/motorized image sets
      filter(OpState != "maintenance") %>%
      arrange(NewLocationID)
    return(dets)
  }
  eoe20s_dets <- detections(eoe20s_allM_skinny, start_date = "2020-07-01", end_date = "2020-09-15")
  eoe20w_dets <- detections(eoe20w_allM_skinny, start_date = "2020-12-01", end_date = "2021-02-01")
  eoe21s_dets <- detections(eoe21s_allM_skinny, start_date = "2021-07-01", end_date = "2021-09-15")
  
  
  #'  ------------------------------------------
  ####  Summarize independent detection events  ####
  #'  ------------------------------------------
  #'  Generate "independent" detection events based on defined amount of time 
  #'  elapsing between sequential images of same species at a camera, then count
  #'  the number of unique detection events of a given species at each camera.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    #'  Generate unique detection events
    det_events <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Flag images of same species at same camera as being a different detection event
      #'  when time since last image of that species is greater than defined time interval
      group_by(Species, NewLocationID) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds!
      ungroup() %>%
      #'  Retain only the first image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup()
      
    #' Count number of detection events per species per camera
    n_dets <- det_events %>%
      group_by(NewLocationID, Species) %>%
      summarise(n_dets = n()) %>%
      ungroup()
    return(n_dets)
  }
  eoe20s_5min_ndets <- unique_detections(eoe20s_dets, elapsed_time = 300) # (5*60 = 300 seconds)
  eoe20w_5min_ndets <- unique_detections(eoe20w_dets, elapsed_time = 300)
  eoe21s_5min_ndets <- unique_detections(eoe21s_dets, elapsed_time = 300)
  #'  List 5min-elapsed detection events
  eoe_5min_list <- list(eoe20s_5min_ndets, eoe20w_5min_ndets, eoe21s_5min_ndets)
  
  eoe20s_30min_ndets <- unique_detections(eoe20s_dets, elapsed_time = 1800) #(30*60 = 1800 seconds)
  eoe20w_30min_ndets <- unique_detections(eoe20w_dets, elapsed_time = 1800)
  eoe21s_30min_ndets <- unique_detections(eoe21s_dets, elapsed_time = 1800)
  #'  List 30min-elapsed detection events
  eoe_30min_list <- list(eoe20s_30min_ndets, eoe20w_30min_ndets, eoe21s_30min_ndets)

  # mismatches <- setdiff(eoe20s_5min_ndets[,1:3], eoe20s_30min_ndets[,1:3])
  
  #'  ---------------------------------
  ####  Hours when detection occurred  ####
  #'  ---------------------------------
  #'  Reduce image data to the hour when it occurred and thin to a single detection
  #'  per hour.
  #'  ---------------------------------
  detection_hour <- function(dets) {
    det_hr <- dets %>%
      #'  Remove empty images
      filter(Species != "none") %>%
      #'  Floor time to the hour of each image
      mutate(floor_DT = floor_date(posix_date_time, unit = "hour")) %>%
      #'  Thin to just the first image per hour for a given species per camera
      group_by(NewLocationID, Species, floor_DT) %>%
      slice(1L) %>%
      ungroup() %>%
      group_by(NewLocationID, Species) %>%
      #'  Count total number of hours each species was detected per camera
      summarize(n_dets = n()) %>%
      ungroup()
    return(det_hr)
  }
  eoe20s_dethr <- detection_hour(eoe20s_dets)
  eoe20w_dethr <- detection_hour(eoe20w_dets)
  eoe21s_dethr <- detection_hour(eoe21s_dets)
  #'  List hour of detection data sets
  eoe_dethr_list <- list(eoe20s_dethr, eoe20w_dethr, eoe21s_dethr)
  
  #'  --------------------
  ####  Correlation test  ####
  #'  --------------------
  #'  Test for correlation between different metrics of relative abundance
  compare_relative_abund <- function(ndets_5min, ndets_30min, n_dethrs, spp) {
    #'  Filter to species of interest
    ndets_5min <- filter(ndets_5min, Species == spp) %>% dplyr::select(n_dets) %>%
      rename(ndets_5min = n_dets)
    ndets_30min <- filter(ndets_30min, Species == spp) %>% dplyr::select(n_dets) %>%
      rename(ndets_30min = n_dets)
    n_dethrs <- filter(n_dethrs, Species == spp) %>% dplyr::select(n_dets) %>%
      rename(dethr = n_dets)
    #'  Combine all datasets
    ra_metrics <- cbind(ndets_5min, ndets_30min, n_dethrs) 
      
    
    #'  Run correlation test using Pearson's correlation
    corr_all <- cor(ra_metrics)
    corr_all <- as.data.frame(round(corr_all, 2))
    
    print(corr_all)
    return(corr_all)
  }
  elk_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "elk")
  moose_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "moose")
  md_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "muledeer")
  wtd_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "whitetaileddeer")
  bunny_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "rabbit_hare")
  
  human_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "human")
  motorized_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "human_motorized")
  cattle_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "cattle_cow")
  
  #'  Saving different metrics of relative abundance index
  #'  Detection hr highly correlated with other metrics for all species & years; r = 0.93 - 1.0
  #'  (5min and 30min metrics less correlated for moose in some years; r = 0.78)
  save(eoe_5min_list, file = "./Data/Relative abundance data/EoE_RelativeN_5minElapsed.RData")
  save(eoe_30min_list, file = "./Data/Relative abundance data/EoE_RelativeN_30minElapsed.RData")
  save(eoe_dethr_list, file = "./Data/Relative abundance data/EoE_RelativeN_HrOfDetection.RData")

  
  #'  ----------------------------------------------------
  ####  Relative abundance indices per sampling occasion  ####
  #'  ----------------------------------------------------
  #'  Number of unique detection events per sampling occasion, based on 5-min vs 
  #'  30-min elapsing between sequential images to define an independent event 
  unique_dets_per_occasion <- function(dets, elapsed_time, startDate, occlength) {
    #'  Generate unique detection events
    det_events30 <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Indicate which 1 week sampling occasion each image belongs to 
      mutate(SamplingOccasion = 1 + as.numeric(Date - as.Date(startDate)) %/% occlength) %>%
      #'  Flag images of same species at same camera as being a different detection event
      #'  when time since last image of that species is greater than defined time interval
      group_by(NewLocationID, Species) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds! 
      ungroup() %>%
      #'  Retain only the first image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup()

    #' Count total number of detection events per species per camera
    n_dets <- det_events %>%
      group_by(NewLocationID, Species, SamplingOccasion) %>%
      summarise(n_dets = n()) %>%
      ungroup()
    return(n_dets)
  }
  eoe20s_5min_ndets_occ <- unique_dets_per_occasion(eoe20s_dets, elapsed_time = 300, startDate = "2020-07-01", occlength = 7) # 7 day sampling occasion
  eoe20w_5min_ndets_occ <- unique_dets_per_occasion(eoe20w_dets, elapsed_time = 300, startDate = "2020-12-01", occlength = 7)
  eoe21s_5min_ndets_occ <- unique_dets_per_occasion(eoe21s_dets, elapsed_time = 300, startDate = "2021-07-01", occlength = 7)
  #'  List 5min-elapsed detection events
  eoe_5min_sampocc_list <- list(eoe20s_5min_ndets_occ, eoe20w_5min_ndets_occ, eoe21s_5min_ndets_occ)
  
  eoe20s_30min_ndets_occ <- unique_dets_per_occasion(eoe20s_dets, elapsed_time = 1800, startDate = "2020-07-01", occlength = 7) # 7 day sampling occasion
  eoe20w_30min_ndets_occ <- unique_dets_per_occasion(eoe20w_dets, elapsed_time = 1800, startDate = "2020-12-01", occlength = 7)
  eoe21s_30min_ndets_occ <- unique_dets_per_occasion(eoe21s_dets, elapsed_time = 1800, startDate = "2021-07-01", occlength = 7)
  #'  List 30min-elapsed detection events
  eoe_30min_sampocc_list <- list(eoe20s_30min_ndets_occ, eoe20w_30min_ndets_occ, eoe21s_30min_ndets_occ)
  
  
  #'  Number of hours when at least one detection occurred per sampling occasion
  det_hour_per_occasion <- function(dets, startDate, occlength) {
    det_hr <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Remove empty images
      filter(Species != "none") %>%
      #'  Indicate which 1 week sampling occasion each image belongs to
      mutate(SamplingOccasion = 1 + as.numeric(Date - as.Date(startDate)) %/% occlength) %>%
      #'  Floor time to the hour of each image
      mutate(floor_DT = floor_date(posix_date_time, unit = "hour")) %>%
      #'  Thin to just the first image per hour for a given species and camera
      group_by(NewLocationID, Species, floor_DT) %>%
      slice(1L) %>%
      ungroup() %>%
      #'  Count total number of hours each species was detected per sampling occasion
      group_by(NewLocationID, Species, SamplingOccasion) %>%
      summarize(n_dets = n()) %>%
      ungroup()
    return(det_hr)
  }
  eoe20s_dethr_occ <- det_hour_per_occasion(eoe20s_dets, startDate = "2020-07-01", occlength = 7)
  eoe20w_dethr_occ <- det_hour_per_occasion(eoe20w_dets, startDate = "2020-12-01", occlength = 7)
  eoe21s_dethr_occ <- det_hour_per_occasion(eoe21s_dets, startDate = "2021-07-01", occlength = 7)
  #'  List hour of detection data sets
  eoe_dethr_sampocc_list <- list(eoe20s_dethr_occ, eoe20w_dethr_occ, eoe21s_dethr_occ)
  
  
  ####  CURRENTLY USING THIS TO FLAG CAMS WHERE CUMSUM FUNCTION ISN'T WORKING RIGHT  ####
  
  #'  Test for correlation between different metrics of relative abundance
  compare_relative_abund <- function(ndets_5min, ndets_30min, n_dethrs, spp) {
    #'  Filter to species of interest
    ndets_5min <- filter(ndets_5min, Species == spp) %>% mutate(datastream = "5min")#%>% dplyr::select(n_dets) %>%
      #rename(ndets_5min = n_dets)
    ndets_30min <- filter(ndets_30min, Species == spp) %>% mutate(datastream = "30min")#%>% dplyr::select(n_dets) %>%
      #rename(ndets_30min = n_dets)
    n_dethrs <- filter(n_dethrs, Species == spp) %>% mutate(datastream = "dethr")#%>% dplyr::select(n_dets) %>%
      #rename(dethr = n_dets)
    #'  Combine all datasets
    ra_metrics <- full_join(ndets_5min, ndets_30min, by = c("NewLocationID", "Species", "SamplingOccasion")) %>%
      full_join(n_dethrs, by = c("NewLocationID", "Species", "SamplingOccasion"))
    return(ra_metrics)
    
    #' #'  Run correlation test using Pearson's correlation
    #' corr_all <- cor(ra_metrics)
    #' corr_all <- as.data.frame(round(corr_all, 2))
    #' 
    #' print(corr_all)
    #' return(corr_all)
  }
  elk_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "elk")
  moose_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "moose")
  md_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "muledeer")
  wtd_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "whitetaileddeer")
  bunny_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "rabbit_hare")
  
  human_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "human")
  vehicl_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "human_vehicle")
  cattle_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_sampocc_list, ndets_30min = eoe_30min_sampocc_list, n_dethrs = eoe_dethr_sampocc_list, spp = "cattle_cow")
  
  
  #'  Saving total hours detected as relative abundance index
  save(eoe_5min_sampocc_list, file = "./Data/Relative abundance data/EoE_RelativeN_5minElapsed_SamplingOcc.RData")
  save(eoe_30min_sampocc_list, file = "./Data/Relative abundance data/EoE_RelativeN_30minElapsed_SamplingOcc.RData")
  save(eoe_dethr_sampocc_list, file = "./Data/Relative abundance data/EoE_RelativeN_HrOfDetection_SamplingOcc.RData")
  
  
  
  
  