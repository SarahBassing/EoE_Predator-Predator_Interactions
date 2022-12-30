  #'  ---------------------------------
  #'  Defining unique detection events
  #'  Sarah Bassing
  #'  ID CRU - Predator Interactions
  #'  Dec. 2020
  #'  ---------------------------------
  #'  Explor the duratin of time used to define a unique detection event. How much
  #'  does it change the number of detections for a given species?
  #'  ---------------------------------
    
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
  
  
  #'  ------------------------
  ###  Filter detection data  ####
  #'  ------------------------
  #'  Filter detection data to focal species and time period of interest
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Count") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time)
      ) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID))
    return(dets)
  }
  eoe20s_dets <- detections(eoe20s_allM, start_date = "2020-07-01", end_date = "2020-09-15")
  eoe20w_dets <- detections(eoe20w_allM, start_date = "2020-12-01", end_date = "2021-02-01")
  eoe21s_dets <- detections(eoe21s_allM, start_date = "2021-07-01", end_date = "2021-09-15")
  
  
  #'  -----------------------------------------
  ####  Generate independent detection events  ####
  #'  -----------------------------------------
  #'  Create a column grouping images into "independent" detection events based 
  #'  on defined amount of time that must elapse between images of the same species 
  #'  at the same camera site. Filter to the first image of each event and then
  #'  split into species-specific detection events.
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    det_events <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Flag images of same species at same camera as being a different detection event
      #'  when time since last image of that species is greater than defined time interval
      group_by(NewLocationID, Species) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds! 
      ungroup() %>%
      #'  Retain only the first image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup()
    
    return(det_events)
  }
  #'  2 minutes must elapse between images of a species
  eoe20s_dets_2min <- unique_detections(eoe20s_dets, elapsed_time = 120)
  eoe20w_dets_2min <- unique_detections(eoe20w_dets, elapsed_time = 120)
  eoe21s_dets_2min <- unique_detections(eoe21s_dets, elapsed_time = 120)
  #'  5 minutes must elapse between images of a species
  eoe20s_dets_5min <- unique_detections(eoe20s_dets, elapsed_time = 300)
  eoe20w_dets_5min <- unique_detections(eoe20w_dets, elapsed_time = 300)
  eoe21s_dets_5min <- unique_detections(eoe21s_dets, elapsed_time = 300)
  #'  10 minutes must elapse between images of a species
  eoe20s_dets_10min <- unique_detections(eoe20s_dets, elapsed_time = 600)
  eoe20w_dets_10min <- unique_detections(eoe20w_dets, elapsed_time = 600)
  eoe21s_dets_10min <- unique_detections(eoe21s_dets, elapsed_time = 600)
  #'  30 minutes must elapse between images of a species
  eoe20s_dets_30min <- unique_detections(eoe20s_dets, elapsed_time = 1800)
  eoe20w_dets_30min <- unique_detections(eoe20w_dets, elapsed_time = 1800)
  eoe21s_dets_30min <- unique_detections(eoe21s_dets, elapsed_time = 1800)
  #'  60 minutes must elapse between images of a species
  eoe20s_dets_60min <- unique_detections(eoe20s_dets, elapsed_time = 3600)
  eoe20w_dets_60min <- unique_detections(eoe20w_dets, elapsed_time = 3600)
  eoe21s_dets_60min <- unique_detections(eoe21s_dets, elapsed_time = 3600)  
  
  #'  List different definitions of independence together by season and filter
  smr20_dets <- list(eoe20s_dets_2min, eoe20s_dets_5min, eoe20s_dets_10min, eoe20s_dets_30min, eoe20s_dets_60min)
  wtr20_dets <- list(eoe20w_dets_2min, eoe20w_dets_5min, eoe20w_dets_10min, eoe20w_dets_30min, eoe20w_dets_60min)
  smr21_dets <- list(eoe21s_dets_2min, eoe21s_dets_5min, eoe21s_dets_10min, eoe21s_dets_30min, eoe21s_dets_60min)
  
  
  #'  Filter detection data by species and count number of unique detection events
  predators_only <- function(dat, indi_time) {
    n_predator_dets <- filter(dat, Species == "bear_black" | Species == "bobcat" | 
                          Species == "coyote" | Species == "mountain_lion" | 
                          Species == "wolf") %>%
      group_by(Species) %>%
      summarize(n = n()) %>%
      ungroup() 
    return(n_predator_dets)
  }
  
  #'  Create table summing unique detections per species and independence definition
  smr20_preds <- lapply(smr20_dets, predators_only) %>%
    do.call(rbind.data.frame, .) %>%
    #'  Add column indicating which definition of independence was used to generate detections
    mutate(Ind_event_time = rep(c("2_min_elapsed", "5_min_elapsed", "10_min_elapsed", 
                                  "30_min_elapsed", "60_min_elapsed"), each = 5)) %>%
    #'  Spread into a wide table
    spread(Ind_event_time, n) %>% 
    relocate(`5_min_elapsed`, .after = `2_min_elapsed`) %>%
    relocate(`10_min_elapsed`, .after = `5_min_elapsed`) %>%
    #'  Calculate average number of detections across different definitions of independence
    mutate(mean_n_dets = rowMeans(subset(., select = c(`2_min_elapsed`, `60_min_elapsed`))),
           #'  Calculate number of detections lost/gained depending on definition of independence
           diff_n_dets = `2_min_elapsed` - `60_min_elapsed`,
           Season = "Summer 2020")
  wtr20_preds <- lapply(wtr20_dets, predators_only) %>%
    do.call(rbind.data.frame, .) %>%
    #'  Add column indicating which definition of independence was used to generate detections
    mutate(Ind_event_time = rep(c("2_min_elapsed", "5_min_elapsed", "10_min_elapsed", 
                                  "30_min_elapsed", "60_min_elapsed"), each = 4)) %>%
    #'  Spread into a wide table
    spread(Ind_event_time, n) %>% 
    relocate(`5_min_elapsed`, .after = `2_min_elapsed`) %>%
    relocate(`10_min_elapsed`, .after = `5_min_elapsed`) %>%
    #'  Calculate average number of detections across different definitions of independence
    mutate(mean_n_dets = rowMeans(subset(., select = c(`2_min_elapsed`, `60_min_elapsed`))),
           #'  Calculate number of detections lost/gained depending on definition of independence
           diff_n_dets = `2_min_elapsed` - `60_min_elapsed`, 
           Season = "Winter 2020-2021")
  smr21_preds <- lapply(smr21_dets, predators_only) %>%
    do.call(rbind.data.frame, .) %>%
    #'  Add column indicating which definition of independence was used to generate detections
    mutate(Ind_event_time = rep(c("2_min_elapsed", "5_min_elapsed", "10_min_elapsed", 
                                  "30_min_elapsed", "60_min_elapsed"), each = 5)) %>%
    #'  Spread into a wide table
    spread(Ind_event_time, n) %>% 
    relocate(`5_min_elapsed`, .after = `2_min_elapsed`) %>%
    relocate(`10_min_elapsed`, .after = `5_min_elapsed`) %>%
    #'  Calculate average number of detections across different definitions of independence
    mutate(mean_n_dets = rowMeans(subset(., select = c(`2_min_elapsed`, `60_min_elapsed`))),
           #'  Calculate number of detections lost/gained depending on definition of independence
           diff_n_dets = `2_min_elapsed` - `60_min_elapsed`, 
           Season = "Summer 2021")
  
  #'  Merge together
  indi_definition_tbl <- rbind(smr20_preds, wtr20_preds, smr21_preds) %>%
    relocate(Season, .after = "Species") %>%
    arrange(Species, Season)
  
  #'  Save tables
  write.csv(indi_definition_tbl, "./Outputs/Tables/Independence_definition_ndets.csv")
  
  
  
  
  