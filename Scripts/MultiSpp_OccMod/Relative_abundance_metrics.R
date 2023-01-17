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
               Species == "human" | Species == "dog_domestic" | #Species == "deer_speciesunknown" | 
               Species == "horse" |  Species == "cattle_cow" | Species == "cat_domestic" |
               Vehicle == "TRUE") %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time),
        Species = ifelse(Vehicle == "TRUE", "human_vehicle", Species)
      ) %>%
      dplyr::select(-Vehicle) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
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
      summarize(n_dets = n())
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
    n_dethrs <- filter(n_dethrs, Species == spp) %>% dplyr::select(c(NewLocationID, n_dets)) %>%
      rename(dethr = n_dets)
    #'  Combine all datasets
    ra_metrics <- cbind(ndets_5min, ndets_30min, n_dethrs) %>% 
      dplyr::select(-NewLocationID)
      
    
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
  
  human_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "human")
  vehicl_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "human_vehicle")
  cattle_corr <- mapply(compare_relative_abund, ndets_5min = eoe_5min_list, ndets_30min = eoe_30min_list, n_dethrs = eoe_dethr_list, spp = "cattle_cow")
  
  #'  Save relative abundance index data
  save(eoe_5min_list, file = "./Data/Relative abundance data/EoE_RelativeN_5minElapsed.RData")
  save(eoe_30min_list, file = "./Data/Relative abundance data/EoE_RelativeN_30minElapsed.RData")
  save(eoe_dethr_list, file = "./Data/Relative abundance data/EoE_RelativeN_HrOfDetection.RData")

  
  
  
  
  