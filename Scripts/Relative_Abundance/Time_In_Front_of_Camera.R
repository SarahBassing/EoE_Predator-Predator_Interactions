  #'  --------------------------------
  #'  Time In Front of Camera (TIFC)
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  August 2023
  #'  --------------------------------
  #'  Calculate relative abundance at each site using the Time In Front of Camera
  #'  (TIFC) approach (Huggard et al. 2018, Warbington & Boyce 2020, Becker et al. 2022).
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R (Data_Formatting folder)
  #'  Sampling effort data generated in Detection_histories_for_occmod.R (MultiSpp_OccMod folder)
  #'  -------------------------------
  
  #'  Load libraries
  library(data.table)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  --------------------------------
  ####  Read & format detection data  ####
  #'  --------------------------------
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_2023-08-09.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_2023-08-09.RData")
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  
  #'  ------------------------
  ####  Filter detection data  ####
  #'  ------------------------
  #'  1) Remove sequential problem images from larger image set
  thin_detections <- function(dets, seqprobs) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    return(skinny_dets)
  }
  eoe20s_allM_skinny <- thin_detections(eoe20s_allM, eoe_seqprob_20s) 
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s) 
  
  #'  2) Filter  data to time period and species of interest
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "File", "Location_Relative_Project", 
                    "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      filter(Species != "none") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time)) %>%
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "elk" | Species == "human" | Species == "moose" |
               Species == "mountain_lion" | Species == "muledeer" | Species == "rabbit_hare" |
               Species == "whitetaileddeer" | Species == "wolf" | Species == "cattle_cow") %>%
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
  eoe21s_dets <- detections(eoe21s_allM_skinny, start_date = "2021-07-01", end_date = "2021-09-15")
  
  #'  ---------------------------------------
  ####  Define series of consecutive images  ####
  #'  ---------------------------------------
  #'  Identify series of consecutive images of the same species following Becker et al. 2022
  #'  "Series" defined as consecutive images with intervals <120 seconds between
  #'  two consecutive images. Intermediate image without the species ends the series.
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    #'  Generate unique detection events
    dat <- arrange(dets, NewLocationID, posix_date_time)
    det_events <- c()
    det_events[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$NewLocationID[i-1] != dat$NewLocationID[i]) det_events[i] = i
      else (if (dat$Species[i-1] != dat$Species[i]) det_events[i] = i
            else (if (difftime(dat$posix_date_time[i], dat$posix_date_time[i-1], units = c("secs")) > elapsed_time) det_events[i] = i
                  else det_events[i] = det_events[i-1]))
    }
    
    det_events <- as.factor(det_events)
    
    #'  Add new column to larger data set
    det_events <- cbind(as.data.frame(dat), det_events)
    
    return(det_events)
  }
  eoe20s_120s_gaps <- unique_detections(eoe20s_dets, elapsed_time = 120) 
  eoe21s_120s_gaps <- unique_detections(eoe21s_dets, elapsed_time = 120)
  #'  List 5min-elapsed detection events
  eoe_5min_list <- list(eoe20s_120s_gaps, eoe21s_120s_gaps)
  
  
  
  
  