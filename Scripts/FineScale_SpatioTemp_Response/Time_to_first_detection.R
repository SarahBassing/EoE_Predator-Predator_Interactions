  #'  ---------------------------------------
  #'  Calculate time to detections
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------------
  #'  Calculate time from a random starting point to the first detection of each
  #'  species. To be used as a comparison to the predator-predator TBD analysis.
  #'  ---------------------------------------
  
  #'  Load libraries
  library(data.table)
  library(lubridate)
  library(chron)
  library(purrr)
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
  
  #'  Load corrected species ID and update larger datasets
  newSppID <- read.csv("./Data/IDFG camera data/questionable_images_doublecheck_SBBupdated.csv") %>%
    mutate(posix_date_time = as.POSIXct(posix_date_time, format="%m/%d/%Y %H:%M", tz="UTC"), 
           Year = as.numeric(format(posix_date_time, "%Y"))) %>%
    distinct() %>%
    dplyr::select(c("CamID", "NewLocationID", "File", "Location_Relative_Project", "posix_date_time", "Species", "NewSpecies", "Year"))
  #'  Filter by year and only observations with misidentified species
  newSppID_20s <- filter(newSppID, Year == "2020") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  newSppID_21s <- filter(newSppID, Year == "2021") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  #'  Snag unique ID for each image
  change_sppID_20s <- as.vector(newSppID_20s$Location_Relative_Project)
  change_sppID_21s <- as.vector(newSppID_21s$Location_Relative_Project)
  
  #'  Function to correc species misclassifications
  remove_obs <- function(dets, prob_images, correctSppID) {
    #'  Grab all rows that have a misID's species
    obs_to_remove <- dets[dets$Location_Relative_Project %in% prob_images,]
    #'  Switch original species ID with correct one
    newspp <- left_join(obs_to_remove, correctSppID, by = c("CamID", "NewLocationID", "File", "Species")) %>%
      dplyr::select(-c("Location_Relative_Project.y", "posix_date_time.y")) %>%
      relocate(NewSpecies, .after = Species) %>%
      dplyr::select(-Species) %>%
      rename(Species = "NewSpecies", Location_Relative_Project = "Location_Relative_Project.x", posix_date_time = "posix_date_time.x")
    #'  Remove these rows from the larger df
    dets_reduced <- anti_join(dets, obs_to_remove)
    #'  Replace with observation containing correct species ID
    dets_updated <- bind_rows(dets_reduced, newspp) %>%
      #'  Arrange everything back in order
      arrange(NewLocationID, posix_date_time, File)
    return(dets_updated)
  }
  eoe20s_allM_new <- remove_obs(eoe20s_allM, prob_images = change_sppID_20s, correctSppID = newSppID_20s)
  eoe21s_allM_new <- remove_obs(eoe21s_allM, prob_images = change_sppID_21s, correctSppID = newSppID_21s)
  
  #'  Re-assign
  eoe20s_allM <- eoe20s_allM_new
  eoe21s_allM <- eoe21s_allM_new
  
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
      dplyr::select("NewLocationID", "CamID", "File", "Location_Relative_Project", "Date", "Time", "posix_date_time", "TriggerMode",
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
  
  #'  Function to generate random sample of times within study period
  random_start_times <- function(start_date, end_date, nobs, firstlast) {
    #'  Create sequence of date-times across focal date range, incrementally changing by second
    datetime_sequence <- seq(as.POSIXct(start_date), as.POSIXct(end_date), by = "1 secs")
    print(head(datetime_sequence))
    
    #'  Sample random date-times from sequence
    datetime_random <- sample(datetime_sequence, size = nobs, replace = TRUE)
    print(head(datetime_random))
    
    #'  Make sure datetime is formatted correctly
    datetime_random <- as.POSIXct(datetime_random, format = "%Y-%m-%d %H:%M:%S")
    
    #'  Generate additional columns to add to datetime
    date_random <- as.Date(datetime_random)
    datetime_random <- as.data.frame(datetime_random)
    category <- rep("random", nobs)
    spp <- rep("NA", nobs)
    
    #'  Bind together into a single dataframe
    random_posixct <- cbind(datetime_random, date_random, category, spp)
    colnames(random_posixct) <- c("posix_date_time", "Date", "Category", "Species")
    print(head(random_posixct))
    
    return(random_posixct)
  }
  # firstlast_random_20s <- random_start_times(start_date = "2020-07-01", end_date = "2020-09-15", nobs = 1000, firstlast = firstlast_img[[1]])
  # firstlast_random_21s <- random_start_times(start_date = "2021-07-01", end_date = "2021-09-15", nobs = 1000, firstlast = firstlast_img[[3]])
  # firstlast_random <- list(firstlast_random_20s, firstlast_random_21s)
  # 
  # NewLocationID_20s <- unique(firstlast_img[[1]]$NewLocationID)
  # NewLocationID_21s <- unique(firstlast_img[[3]]$NewLocationID)
  # random_times_list_20s <- random_times_list_21s <- list()
  # for(i in 1:length(NewLocationID)){
  #   random_times_list_20s[[i]] <- random_start_times(start_date = "2020-07-01", end_date = "2020-09-15", nobs = 1000, firstlast = firstlast_img[[1]])
  #   random_times_list_21s[[i]] <- random_start_times(start_date = "2021-07-01", end_date = "2021-09-15", nobs = 1000, firstlast = firstlast_img[[3]])
  # }
  # names(random_times_list) <- NewLocationID
  # firstlast_random_20s <- map_df(random_times_list, ~as.data.frame(.x), .id="id")
  
  #'  Function to resample random times for each NewLocationID and add to larger first/last dataset
  add_random_times <- function(firstlast, start_date, end_date, nobs) {
    #'  Snag NewLocationID for each camera site
    NewLocationID <- unique(firstlast$NewLocationID)
    
    #'  Create empty list to fill with randomly sampled datetimes
    random_times_list <- list()
    
    #'  Loop through each camera and sample a collection of random datetimes
    for(i in 1:length(NewLocationID)){
      random_times_list[[i]] <- random_start_times(start_date = start_date, end_date = end_date, nobs = nobs, firstlast = firstlast)
    }
    #'  Name the list based on NewLocationID
    names(random_times_list) <- NewLocationID
    #'  Reduce to a single large data frame with additional column mapping back to NewLocationID
    random_df <- map_df(random_times_list, ~as.data.frame(.x), .id="id") %>%
      rename(NewLocationID = "id")
    
    #'  Join with first/last detection data
    firstlast_random <- full_join(firstlast, random_df, by = c("NewLocationID", "Date", "posix_date_time", "Species", "Category")) %>%
      arrange(NewLocationID, posix_date_time)
    
    return(firstlast_random)
  }
  firstlast_img_random_20s <- add_random_times(start_date = "2020-07-01", end_date = "2020-09-15", nobs = 1000, firstlast = firstlast_img[[1]])
  firstlast_img_random_21s <- add_random_times(start_date = "2021-07-01", end_date = "2021-09-15", nobs = 1000, firstlast = firstlast_img[[3]])
  
    
    
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
  full_predator_sequences <- lapply(firstlast_random, thin_dat_by_category) 
  
  
  
  
  
    
  
  
  
  
  