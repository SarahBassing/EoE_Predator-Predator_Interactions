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
  
  #'  Time between detections data
  load("./Data/Time_btwn_Detections/TBD_all_predator_pairs_2023-05-22.RData")
  
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
    # skinny_dets <- skinny_dets %>%
    #   filter(NewLocationID != "GMU6_U_29")
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
  firstlast_img_random <- list(firstlast_img_random_20s, firstlast_img_random_21s)
    
    
  #'  --------------------------------------
  ####  Filter to specific seets of images  ####
  #'  --------------------------------------
  #'  1. Thin image set to only single image within sequential group of each
  #'     category. Retain first image from predator and random categories, last
  #'     image from other category.
  #'  2. Flag & retain sequences when a random time is followed by a predator img.
  #'  3. Double check problem date ranges don't fall within time from random start
  #'     to first image of a predator detection.
  #'  4. Reduce to only sites where we have back-to-back predator detections as 
  #'     well to keep spatial extent of analyses consistent.
  #'  5. Calculate time-between-detections of different species
  #'  --------------------------------------

  #'  1) Group multiple detection events of same category (but of different species)
  #'  when they occur sequentially, then reduce groups to a single observation 
  #'  (e.g., we only care that a predator sequence was broken up but don't need 
  #'  all "other" detections).
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
    #'  Save the first of a group of predator detections
    predspp <- capdata[capdata$Category == "Predator",] %>%
    group_by(caps_new) %>%
    slice(1L) %>%
    ungroup()
    #'  Save the first of a group of random detections
    randtime <- capdata[capdata$Category == "random",] %>%
      group_by(caps_new) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Det_type = "random")
    #'  Save the last of a group other detections (non-predator or non-random)
    lasteverythingelse <- capdata[capdata$Category == "Other",] %>%
      group_by(caps_new) %>%
      slice_tail() %>%
      ungroup()
    #'  Combine into final data set
    dets <- rbind(predspp, randtime, lasteverythingelse) %>%
      arrange(NewLocationID, posix_date_time)
    return(dets)
  }
  thinned_sequences <- lapply(firstlast_img_random, thin_dat_by_category) 
  
  tst <- thinned_sequences[[1]]

  #'  2) Filter to just the random time followed by a predator time
  random_to_predator <- function(dat) {
    rnd_pred <- dat %>%
      group_by(NewLocationID) %>%
      #'  Flag pairings of random time followed by a predator deteciton
      mutate(keep = ifelse(Category == "random" & lead(Category == "Predator"), "Y", "N"),
             keep = ifelse(Category == "Predator" & lag(Category == "random"), "Y", keep)) %>%
      ungroup() %>%
      #'  Filter to just those pairings
      filter(keep == "Y")
    return(rnd_pred)
  } 
  random_predators <- lapply(thinned_sequences, random_to_predator)
  
  #'  3) Review cameras with problem time periods and the dates of random - predator 
  #'  detections to make sure elapsed time to detections do not encompass
  #'  problematic time periods (don't want inoperable cameras to bias T2D data)
  problem_dates_and_t2b_dets <- function(rand_preds, seqprobs, start_date, end_date) {
    #'  Format and thin problem camera data
    prob_date_range <- seqprobs %>%
      mutate(Date = as.Date(Date, format = "%d-%b-%Y")) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      group_by(NewLocationID) %>%
      #'  Group problem time periods separately if camera is operable for > 1 day 
      #'  within larger problem date range
      mutate(New_problem = cumsum(c(1, diff(Date) > 1))) %>%
      #'  Filter to just start and end of each problem period
      group_by(NewLocationID, New_problem) %>%
      filter(row_number()==1 | row_number()==n()) %>%
      ungroup()
    
    #'  Identify cameras that were temporarily inoperable and also had b2b predator detections
    reduced_prob_dates <- prob_date_range[prob_date_range$NewLocationID %in% rand_preds$NewLocationID,]
    #'  Identify cameras that had b2b predator detections but were also inoperable for some period of time
    t2b_prob_cams <- rand_preds[rand_preds$NewLocationID %in% prob_date_range$NewLocationID,]
    
    #'  Create data frame with start and end of problem dates
    prob_dates <- reduced_prob_dates %>%
      group_by(NewLocationID) %>%
      mutate(Prob_start = posix_date_time,
             Prob_end = lead(posix_date_time)) %>%
      ungroup() %>%
      #'  Retain correct date ranges in few instances where there are multiple 
      #'  different problem time periods at a camera site
      group_by(NewLocationID, New_problem) %>%
      slice(1L) %>%
      ungroup() %>%
      dplyr::select(c("NewLocationID", "Prob_start", "Prob_end")) 
    
    #'  Create data frame with start and end of t2d data that potentially overlap problem dates
    t2d <- t2b_prob_cams %>%
      group_by(NewLocationID) %>%
      mutate(Random_time = posix_date_time,
             Predator_detection = lead(posix_date_time),
             Detected_predator = lead(Species)) %>%
      ungroup() %>%
      filter(Category == "random") %>%
      dplyr::select(c("NewLocationID", "Species", "Category", "Random_time", "Predator_detection", "Detected_predator"))
    
    #'  Flag which observations fall within problematic date range. These snuck in 
    #'  and should be removed from data set.
    t2d_prob_overlap <- t2d %>%
      left_join(prob_dates, by = "NewLocationID") %>%
      mutate(sneaky_mf = ifelse(Random_time >= Prob_start & Random_time <= Prob_end, 1, NA),
             sneaky_mf = ifelse(Predator_detection >= Prob_start & Predator_detection <= Prob_end, 1, sneaky_mf)) %>%
      filter(!is.na(sneaky_mf))
    print(t2d_prob_overlap)
    
    #'  Create data frame of just images to remove based on those listed in t2d_prob_overlap
    #'  (need to remove predator obs AND random starting times assocaited w/ them)
    problem_cam <- rand_preds[rand_preds$NewLocationID %in% t2d_prob_overlap$NewLocationID,]
    problem_random_obs <- problem_cam[problem_cam$posix_date_time %in% t2d_prob_overlap$Random_time,]
    problem_pred_obs <- problem_cam[problem_cam$posix_date_time %in% t2d_prob_overlap$Predator_detection,]
    problem_obs <- rbind(problem_random_obs, problem_pred_obs) %>% arrange(posix_date_time)
    #'  Double check correct images are being removed (posix_date_time should match
    #'  Random_time and Predator_detection columns from t2d_prob_overlap)
    print(problem_obs)
    
    #'  Remove observations from larger time-to-detection data frame if they fall
    #'  within a problematic time period at a camera site
    thinned_t2d_dat <- anti_join(rand_preds, problem_obs)
      
    return(thinned_t2d_dat)
  }
  rand_predators_20s <- problem_dates_and_t2b_dets(random_predators[[1]], eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15")
  rand_predators_21s <- problem_dates_and_t2b_dets(random_predators[[2]], eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15")

  #'  4) Remove camera data from any cameras not included in the Time_between_detections
  #'  analysis to keep any site-specific heterogeneity (e.g., local density) consistent
  #'  across analyses
  tbd_cams_20s <- unique(tbd_pred_pairs_all$NewLocationID[tbd_pred_pairs_all$Year == "Smr20"])
  tbd_cams_21s <- unique(tbd_pred_pairs_all$NewLocationID[tbd_pred_pairs_all$Year == "Smr21"])
  rand_predators_20s_skinny <- rand_predators_20s[rand_predators_20s$NewLocationID %in% tbd_cams_20s,]
  rand_predators_21s_skinny <- rand_predators_21s[rand_predators_21s$NewLocationID %in% tbd_cams_21s,]
  
  random_predators_skinny <- list(rand_predators_20s_skinny, rand_predators_21s_skinny)
  
  
  #'  ----------------------------------------
  ####  Calculate times to detection events   ####
  #'  ----------------------------------------
  #'  Function to calculate time to detection event from a random starting point
  #'  Data structured so only first random time followed by first image of focal 
  #'  species are included in data frame.
  t2d <- function(detection_data, unittime) {
    detection_data <- detection_data %>%
      group_by(NewLocationID) %>%
      mutate(TimeSinceLastDet = difftime(posix_date_time, lag(posix_date_time), units = unittime),
             TimeSinceLastDet = ifelse(Det_type == "random", 0, TimeSinceLastDet),
             MinutesSinceLastDet = round(TimeSinceLastDet, 4),
             HoursSinceLastDet = round((TimeSinceLastDet/60), 2),
             DaysSinceLastDet = round((TimeSinceLastDet/1440), 2)) %>%
      ungroup() %>%
      filter(Category == "Predator") %>%
      dplyr::select(-c(TriggerMode, OpState, Category, Count, caps_new, Det_type, keep)) 
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  t2d_preds <- lapply(random_predators_skinny, t2d, unittime = "min")
  
  #'  Add column for year
  t2d_preds[[1]]$Year <- "Smr20"
  t2d_preds[[2]]$Year <- "Smr21"
  t2d_preds_all <- rbind(t2d_preds[[1]], t2d_preds[[2]])

  save(t2d_preds_all, file = paste0("./Data/Time_btwn_Detections/T2D_all_predators_", Sys.Date(), ".RData"))  
  
    
  #'  ----------------------------------------
  ####  Summary stats and data visualization  ####
  #'  ----------------------------------------
  #'  Total average
  mean(t2d_preds_all$MinutesSinceLastDet, na.rm = TRUE); sd(t2d_preds_all$MinutesSinceLastDet, na.rm = TRUE)
  mean(t2d_preds_all$HoursSinceLastDet, na.rm = TRUE); sd(t2d_preds_all$HoursSinceLastDet, na.rm = TRUE)
  
  #'  Summary stats for average time to first detection per predator species
  avg_t2d <- t2d_preds_all %>%
    group_by(Species) %>%
    summarise(total_obs = n(),
              mean_tbd_hr = round(mean(HoursSinceLastDet), 2),
              sd = round(sd(HoursSinceLastDet), 2),
              se = round(sd(HoursSinceLastDet)/sqrt(total_obs), 2)) %>%
    ungroup() %>%
    arrange(desc(total_obs))
  
  # write.csv(avg_tbd_focal_pairs, "./Outputs/Tables/Summary_stats_TBD.csv")
  
  hist(t2d_preds_all$HoursSinceLastDet, breaks = 50, main = "Elapsed time to first detection of predators", xlab = "Hours between detection events")
  hist(t2d_preds_all$DaysSinceLastDet, breaks = 50, main = "Elapsed time to first detection of predators", xlab = "Days between detection events")

  #'  Visualize data
  bear_his <- t2d_preds_all[t2d_preds_all$Species == "bear_black",] %>%
    ggplot(aes(x=HoursSinceLastDet, fill = Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10) +
    scale_fill_manual(values="#69b3a2")
  bob_his <- t2d_preds_all[t2d_preds_all$Species == "bobcat",] %>%
    ggplot(aes(x=HoursSinceLastDet, fill = Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10) +
    scale_fill_manual(values="#404080")
  coy_his <- t2d_preds_all[t2d_preds_all$Species == "coyote",] %>%
    ggplot(aes(x=HoursSinceLastDet, fill = Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10) +
    scale_fill_manual(values="#69b3a2")
  lion_his <- t2d_preds_all[t2d_preds_all$Species == "mountain_lion",] %>%
    ggplot(aes(x=HoursSinceLastDet, fill = Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10) +
    scale_fill_manual(values="#69b3a2")
  wolf_his <- t2d_preds_all[t2d_preds_all$Species == "wolf",] %>%
    ggplot(aes(x=HoursSinceLastDet, fill=Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10) +
    scale_fill_manual(values="#404080")
  plot(bear_his)
  plot(bob_his)
  plot(coy_his)
  plot(lion_his)
  plot(wolf_his)
  
  #'  Boxplots of all data, organized by sample size
  pred_box <- ggplot(t2d_preds_all, aes(x = Species, y=HoursSinceLastDet, fill=Species)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(data = avg_t2d,
              aes(Species, Inf, label = paste("n =",total_obs)), vjust = 8) +
    ggtitle("Boxplots summarizing time to first detection of predators") +
    xlab("Predator speceis") + ylab("Hours to first detection") +
    labs(fill = "Species")
  plot(pred_box)
  
  
  
  