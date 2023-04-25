  #'  ---------------------------------------
  #'  Calculate time between detections
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------------
  #'  Calculate time-between-detections of one predator species following the 
  #'  detection of another predator species at the same camera site.
  #'  ---------------------------------------
  

  #'  Load libraries
  library(data.table)
  library(lubridate)
  library(chron)
  # library(sp)
  # library(raster)
  library(tidyverse)
  
  #'  --------------------------------
  ####  Read & format detection data  ####
  #'  --------------------------------
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe20w_sequential_probimgs.RData")
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
  eoe20w_allM_skinny <- thin_detections(eoe20w_allM, eoe_seqprob_20w)
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s)
  
  #'  2) Filter detection data to time period and species of interest
  #'  Species of interest are ungulates, humans/vehicles, and domestic animals
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      filter(Species != "none") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time),
        Species = ifelse(Vehicle == "TRUE", "human_motorized", Species),
        Category = ifelse(Species == "bobcat" | Species == "bear_black" | 
                                   Species == "coyote" | Species == "mountain_lion" | 
                                   Species == "wolf", "Predator", "Other")) %>%
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
      ungroup()
    return(det_events)
  }
  eoe20s_5min_dets <- unique_detections(eoe20s_dets, elapsed_time = 300) # (5*60 = 300 seconds)
  eoe20w_5min_dets <- unique_detections(eoe20w_dets, elapsed_time = 300)
  eoe21s_5min_dets <- unique_detections(eoe21s_dets, elapsed_time = 300)
  #'  List 5min-elapsed detection events
  eoe_5min_list <- list(eoe20s_5min_dets, eoe20w_5min_dets, eoe21s_5min_dets)
  
  #'  Filter data to the first/last image from each unique detection event
  first_last_image <- function(dets) {
    #'  First image of each detection event
    firstimg <- dets[dets$Category == "Predator",] %>%
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Det_type = "first") %>%
      arrange(NewLocationID, posix_date_time)
    #'  Last image of each detection event
    lastimg <- dets[dets$Category == "Predator",] %>%
      group_by(NewLocationID, Species, det_events) %>%
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last") %>%
      arrange(NewLocationID, posix_date_time)
    #'  Last image of all OTHER detection events
    lastother <- dets[dets$Category == "Other",] %>%
      group_by(NewLocationID, Species, det_events) %>%
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last") %>%
      arrange(NewLocationID, posix_date_time)
    #'  Merge last image of each predator/other and first image of each predator
    firstlast_img <- rbind(firstimg, lastimg, lastother) %>%
      arrange(NewLocationID, posix_date_time)
    return(firstlast_img)
  }
  firstlast_img <- lapply(eoe_5min_list, first_last_image)
  
  #'  Detections with only one image get duplicated b/c its the first & last image
  #'  Identify duplicate images (ignoring last column where they differ) and filter
  #'  to just one image per detection event
  drop_duplicates <- function(dets) {
    dups <- dets %>%
      dplyr::select("NewLocationID", "Date", "Time", "posix_date_time",
                    "TriggerMode", "Species", "Category", "Det_type") %>%
      group_by_at(vars(-Det_type)) %>%
      filter(n() > 1) %>%
      filter(Det_type == "last")
    #'  Remove the duplicate images from the larger data set & arrange by date/time/location
    no_dups <- anti_join(dets, dups) %>%
      arrange(NewLocationID, posix_date_time)
    return(no_dups)
  }
  firstlast_img <- lapply(firstlast_img, drop_duplicates)
  
  #'  -------------------------
  ####  Filter detection data  ####
  #'  -------------------------
  #'  Group multiple detection events of same category (but of different species) 
  #'  when they occur sequentially, then reduce groups of "other" category to a 
  #'  single observation (e.g., we only care that a predator sequence was broken
  #'  up but don't need all "other" detections).
  thin_dat <- function(dets) {
    dat <- arrange(dets, NewLocationID, posix_date_time) %>%
      dplyr::select(-det_events)
    caps_new <- c()
    caps_new[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$NewLocationID[i-1] != dat$NewLocationID[i]) caps_new[i] = i
      else(if (dat$Category[i-1] != dat$Category[i]) caps_new[i] = i
           else caps_new[i] = caps_new[i-1])
    }
    
    caps_new <- as.factor(caps_new)
    
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), caps_new)
    
    #'  Remove all extra detections when multiple detections of same category occur in a row
    predspp <- capdata[capdata$Category == "Predator",] #%>%
      # group_by(caps_new) %>%
      # slice(1L) %>%
      # ungroup()
    lasteverythingelse <- capdata[capdata$Category != "Predator",] %>%
      group_by(caps_new) %>% 
      slice_tail() %>%
      ungroup()
    #'  Combine into final data set
    dets <- rbind(predspp, lasteverythingelse) %>%
      arrange(NewLocationID, posix_date_time)
    return(dets)
  }
  full_predator_sequences <- lapply(firstlast_img, thin_dat) #firstlast_img_skinny
  
  #'  Flag sequential detections of two different predators
  flag_sequential_predators <- function(dets) {
    #'  Assign same ID to all detection events from the same camera
    dat <- arrange(dets, NewLocationID, posix_date_time)
    cam <- c()
    cam[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$NewLocationID[i-1] != dat$NewLocationID[i]) cam[i] = i
      else cam[i] = cam[i-1]
    }
    #'  Identify images where predator spp2 is detected right after predator spp1
    second_pred <- c()
    second_pred[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] != dat$Category[i]) second_pred[i] = "N"
      else(if (dat$Species[i-1] != dat$Species[i]) second_pred[i] = "Y"
           else(second_pred[i] = "N"))
      
    }
    #'  Identify images where predator spp1 is detected right before predator spp2
    first_pred <- c()
    first_pred[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] != dat$Category[i]) first_pred[i-1] = "N"
      else(if (dat$Species[i-1] != dat$Species[i]) first_pred[i-1] = "Y"
           else(first_pred[i-1] = "N"))
      
    }
    #'  Add "N" to very end of vector so the length matches number of rows in dets
    first_pred[nrow(dat)] <- "N"
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), cam, second_pred, first_pred) %>%
      #'  Make sure no "Other" category observations get labeled "Y"
      mutate(second_pred = ifelse(second_pred == "Y" & Category == "Other", "N", second_pred),
             first_pred = ifelse(first_pred == "Y" & Category == "Other", "N", first_pred),
             #'  Create column flagging detections of interest
             pred_pair = ifelse(second_pred == "Y", "Y", "N"),
             pred_pair = ifelse(first_pred == "Y", "Y", pred_pair)) %>%
      #'  Drop extra columns
      dplyr::select(-c(cam, second_pred, first_pred)) %>%
      #'  Retain only observations of two different predators detected in a row
      filter(pred_pair == "Y")
    return(capdata)
  }
  predator_pairs <- lapply(full_predator_sequences, flag_sequential_predators) 
  
  
  #'  ---------------------------------------------
  ####  Calculate times between detection events   ####
  #'  ---------------------------------------------
  #'  Function to calculate time between detection events of two focal species
  #'  Data structured so only last image of spp1 and first image of spp2 per
  #'  detection event are included in data frame.
  tbd <- function(detection_data, spp1, unittime) {
    #'  Create empty vector to be filled
    detection_data$TimeSinceLastDet <- c()
    #'  Fill first element of the vector to get it started
    detection_data$TimeSinceLastDet[1] <- 0
    #'  Loop through each row to calculate elapsed time since previous detection
    for (i in 2:nrow(detection_data)){
      #'  If previous detection was spp1, set time to 0
      if (detection_data$Category[i-1] == spp1) detection_data$TimeSinceLastDet[i] = 0
      #'  If current detection is spp2 and follows detection of spp1, calculate
      #'  the difference in time from previous detection to current detection
      if (detection_data$Category[i] != spp1) detection_data$TimeSinceLastDet[i] = difftime(detection_data$DateTime[i], detection_data$DateTime[i-1], units = unittime)
    }
    #'  Retain only prey observations (don't need the actual predator detections)
    # detection_data <- filter(detection_data, Category == "Prey")
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_pred.prey_smr <- tbd(predator_pairs, spp1 = "Predator", unittime = "min")
  
  
  
  
  
  
  
  