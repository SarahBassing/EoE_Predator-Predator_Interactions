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
  
  #'  -----------------------------------------
  ####  Generate independent detection events  ####
  #'  -----------------------------------------
  #'  Generate "independent" detection events based on defined amount of time 
  #'  elapsing between sequential images of same species at a camera, then count
  #'  the number of unique detection events of a given species at each camera.
  #'  Currently using 5 min as interval (5*60 = 300 seconds).
  #'  Then filter to the first and last image of each predator detection and 
  #'  just the last image of all other species.
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
  
  #'  --------------------------------------
  ####  Filter to specific pairs of images  ####
  #'  --------------------------------------
  #'  1. Thin image set to only single "other" image within sequential group of 
  #'     other images.
  #'  2. Flag pairs of different predator species detected back-to-back within 
  #'     the same set of sequential detection events.
  #'  3. Identify instances where order of last pred1 - first pred2 is incorrect 
  #'     owing to duplicated observations for same image.
  #'  4. Remove these images from larger last/first data set.
  #'  5. Calculate time-between-detections of different species
  #'  6. Thin image set to just detections of different predator species
  #'  --------------------------------------
  
  #' #'  Detections with only one image get duplicated b/c its the first & last image
  #' #'  Identify duplicate images (ignoring last column where they differ) and filter
  #' #'  to just one image per detection event
  #' drop_duplicates <- function(dets) {
  #'   dups <- dets %>%
  #'     dplyr::select("NewLocationID", "Date", "Time", "posix_date_time",
  #'                   "TriggerMode", "Species", "Category", "Det_type") %>%
  #'     group_by_at(vars(-Det_type)) %>%
  #'     filter(n() > 1) %>%
  #'     filter(Det_type == "last")
  #'   #'  Remove the duplicate images from the larger data set & arrange by date/time/location
  #'   no_dups <- anti_join(dets, dups) %>%
  #'     arrange(NewLocationID, posix_date_time)
  #'   return(no_dups)
  #' }
  #' firstlast_img <- lapply(firstlast_img, drop_duplicates)
  
  #'  Group multiple detection events of same category (but of different species) 
  #'  when they occur sequentially, then reduce groups of "other" category to a 
  #'  single observation (e.g., we only care that a predator sequence was broken
  #'  up but don't need all "other" detections).
  thin_dat_by_category <- function(dets) {
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
  full_predator_sequences <- lapply(firstlast_img, thin_dat_by_category) #firstlast_img_skinny
  
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
      #'  Rearrange observations so that duplicate images with different species
      #'  are ordered by sequence of when each species first arrived and left
      #'  i.e., arrange by predator that was detected first and is leaving (last)
      #'  followed by predator that was detected second and just showed up (first)
      group_by(NewLocationID) %>%
      arrange(posix_date_time) %>% #, desc(Det_type)
      ungroup() %>%
      arrange(NewLocationID) %>%
      #'  Make sure no "Other" category observations get labeled "Y"
      mutate(second_pred = ifelse(second_pred == "Y" & Category == "Other", "N", second_pred),
             first_pred = ifelse(first_pred == "Y" & Category == "Other", "N", first_pred),
             #'  Create column flagging detections of interest
             pred_pair = ifelse(second_pred == "Y", "Y", "N"),
             pred_pair = ifelse(first_pred == "Y", "Y", pred_pair),
             same_time = ifelse(lag(posix_date_time) == posix_date_time, "same", "diff"),
             same_time = ifelse(is.na(same_time), "diff", same_time),
             #'  Create a unique ID so I can more easily remove specific observations
             uniqueID = paste0(NewLocationID, "_", posix_date_time, "_", Species, "_", Det_type))
    return(capdata)
  }
  predator_pairs <- lapply(full_predator_sequences, flag_sequential_predators) 
  
  #'  Pull out predator sequences that contain observations of different species
  #'  detected at the same exact time. Identify which images need to be removed
  #'  b/c order of last/first species detected gets messy, especially when only
  #'  one image of both species so first/last image is same for both species. UGH.
  id_duplicate_times <- function(capdata) {
    double_obs <-  capdata %>%
      group_by(caps_new) %>%
      #'  Filter to any image within group of where... 
      #'  at least one image has same time as another
      filter(any(same_time == "same")) %>%
      #'  and at least one is part of a predator pairing
      filter(any(pred_pair == "Y")) %>%
      ungroup() %>%
      group_by(NewLocationID) %>%
      arrange(posix_date_time, desc(Det_type)) %>%
      ungroup() %>%
      arrange(NewLocationID)
    return(double_obs)
  } 
  double_dets <- lapply(predator_pairs, id_duplicate_times)
  
  #'  List all uniqueIDs of observations that need to be removed in a giant vector
  rm_20s <- c("GMU10A_P_104_2020-09-03 14:52:14_bobcat_first", "GMU10A_P_104_2020-09-03 14:52:15_coyote_last", "GMU10A_P_15_2020-08-07 00:27:15_bobcat_first", "GMU10A_P_23_2020-07-22 01:35:58_bobcat_first",
    "GMU10A_P_40_2020-08-03 00:16:27_coyote_last", "GMU10A_P_40_2020-08-03 00:16:27_wolf_first", "GMU10A_P_41_2020-07-16 22:01:43_coyote_last", "GMU10A_P_5_2020-09-15 07:27:58_bobcat_first",
    "GMU10A_P_59_2020-09-14 06:53:24_coyote_last", "GMU10A_P_86_2020-07-27 01:47:41_wolf_first", "GMU10A_P_86_2020-07-27 01:47:42_coyote_last", "GMU10A_P_86_2020-08-11 02:31:28_bobcat_first",
    "GMU10A_P_86_2020-08-11 02:31:44_coyote_last", "GMU10A_P_86_2020-08-11 02:40:03_coyote_last", "GMU10A_P_86_2020-08-11 02:40:03_bobcat_first", "GMU6_P_17_2020-07-01 22:36:24_coyote_last", 
    "GMU6_P_17_2020-07-01 22:36:25_mountain_lion_last", "GMU6_P_18_2020-08-15 22:20:59_bobcat_first", "GMU6_P_18_2020-08-15 22:21:07_coyote_last", "GMU6_P_37_2020-07-28 00:11:29_coyote_first",
    "GMU6_P_37_2020-07-28 00:11:30_wolf_last", "GMU6_P_38_2020-08-28 21:03:42_wolf_last", "GMU6_P_38_2020-08-28 21:03:42_wolf_first", "GMU6_P_45_2020-08-01 03:56:05_coyote_last", 
    "GMU6_P_56_2020-08-23 23:24:19_coyote_last", "GMU6_P_56_2020-08-28 00:58:18_bobcat_first", "GMU6_P_56_2020-08-28 00:58:21_coyote_last", "GMU6_P_58_2020-07-16 02:08:36_coyote_first", 
    "GMU6_P_58_2020-07-16 02:08:47_wolf_last", "GMU6_P_66_2020-07-31 22:16:15_bobcat_first", "GMU6_P_66_2020-07-31 22:16:18_coyote_last", "GMU6_P_94_2020-07-10 05:24:20_bobcat_last", 
    "GMU6_P_94_2020-09-09 01:52:26_mountain_lion_first", "GMU6_P_94_2020-09-09 01:53:41_bobcat_last", "GMU6_U_130_2020-08-21 23:57:41_wolf_first", "GMU6_U_130_2020-08-21 23:58:48_coyote_last",
    "GMU6_U_90_2020-09-04 13:27:34_mountain_lion_first", "GMU6_U_90_2020-09-04 13:27:50_bear_black_last")
  rm_21s <- c()
  
  
  
  

  #'  Function to remove specific observations that screw up the last-first detection order
  remove_obs <- function(pred_pairs, rm_obs) {
    #'  Grab all rows that meet criteria based on uniqueID
    obs_to_remove <- pred_pairs[pred_pairs$uniqueID %in% rm_obs,]
    #'  Remove these rows from the larger df
    reduced_predator_pairs <- anti_join(pred_pairs, obs_to_remove)
    return(reduced_predator_pairs)
  }
  predator_pairs_20s <- remove_obs(predator_pairs[[1]], rm_20s)
  predator_pairs_20s <- remove_obs(predator_pairs[[3]], rm_21s)
  #'  Being lazy and not doing this for wtr20 data right now
  
  #'  Re-list cleaned predator-pair data sets
  predator_pairs <- list(predator_pairs_20s, predator_pairs[[2]], predator_pairs_21s)
  
  ####### need to fix last/first column so they are every other - double check they are correct

  

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
  
  
  
  
  ### eventually thin to just predator data, no other observations in between
  
  
  