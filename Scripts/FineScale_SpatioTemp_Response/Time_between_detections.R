  #'  ---------------------------------------
  #'  Calculate time between detections
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------------
  #'  Filter detection data down to just back-to-back pairs for focal species.
  #'  Calculate time-between-detections of one predator species following the 
  #'  detection of another predator species at the same camera site.
  #'  Summarize and visualize elapsed time between detections.
  #'  ---------------------------------------
  
  #'  Load libraries
  library(data.table)
  library(lubridate)
  library(chron)
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
      arrange(posix_date_time, File) %>% #, desc(Det_type)
      ungroup() %>%
      arrange(NewLocationID) %>%
      #'  Make sure no "Other" category observations get labeled "Y"
      mutate(second_pred = ifelse(second_pred == "Y" & Category == "Other", "N", second_pred),
             first_pred = ifelse(first_pred == "Y" & Category == "Other", "N", first_pred),
             #'  Change very last/very first predator image at a site to "N" in the few 
             #'  instances where it looks like sequential detections of different predator 
             #'  species but it's b/c it's a new camera location
             second_pred = ifelse(lag(NewLocationID) != NewLocationID & second_pred == "Y", "N", second_pred),
             first_pred = ifelse(NewLocationID != lead(NewLocationID) & first_pred == "Y", "N", first_pred),
             #'  Create column flagging detections of interest
             pred_pair = ifelse(second_pred == "Y", "Y", "N"),
             pred_pair = ifelse(first_pred == "Y", "Y", pred_pair),
             #'  Flag images with identical date/times
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
      arrange(posix_date_time, File) %>% #, desc(Det_type)
      ungroup() %>%
      arrange(NewLocationID)
    return(double_obs)
  } 
  double_dets <- lapply(predator_pairs, id_duplicate_times)
  
  #'  List all uniqueIDs of observations that need to be removed in a giant vector
  #'  CURRENT LIST is just to make my life easier for now - pretty sure these are incorrect
  #'  species classifications that need to be corrected. Once corrected I can revert to rm_20s <- c(NA)
  rm_20s <- c(#"GMU10A_P_59_2020-09-14 06:53:24_coyote_last", # this one actually needs to be removed - duplicated single detection of coyote but next row is a new camera so need to remove this "last" image
              "GMU6_P_14_2020-08-28 08:56:54_coyote_first", "GMU6_P_14_2020-08-28 08:56:55_bear_black_first",
              "GMU6_P_14_2020-08-28 08:58:35_coyote_last", "GMU6_P_14_2020-08-28 09:00:02_bear_black_last",
              "GMU10A_P_89_2020-08-22 00:32:19_coyote_last", "GMU10A_P_89_2020-08-22 00:32:20_bobcat_last",
              "GMU10A_P_89_2020-08-22 00:32:20_bobcat_first", "GMU10A_P_89_2020-08-22 10:04:19_coyote_first")
              #"GMU10A_P_104_2020-09-03 14:52:14_bobcat_first", "GMU10A_P_104_2020-09-03 14:52:15_coyote_last", "GMU10A_P_15_2020-08-07 00:27:16_mountain_lion_last",
              #"GMU10A_P_23_2020-07-22 01:35:58_bobcat_last", "GMU10A_P_23_2020-07-22 01:36:01_coyote_last", "GMU10A_P_40_2020-08-03 00:16:27_coyote_first", "GMU10A_P_40_2020-08-03 00:16:27_wolf_last",
              #"GMU10A_P_86_2020-07-27 01:47:41_wolf_last",  #"GMU10A_P_5_2020-09-15 07:27:58_bobcat_last",
              #"GMU10A_P_86_2020-07-27 01:47:42_coyote_last", "GMU10A_P_86_2020-08-11 02:31:28_bobcat_last", "GMU10A_P_86_2020-08-11 02:31:44_coyote_last",
              #"GMU10A_P_86_2020-08-11 02:40:03_bobcat_last", "GMU10A_P_86_2020-08-11 02:40:03_coyote_last", "GMU6_P_17_2020-07-01 22:36:24_coyote_last", "GMU6_P_17_2020-07-01 22:36:25_mountain_lion_last",
              #"GMU6_P_18_2020-08-15 22:20:59_bobcat_last", "GMU6_P_18_2020-08-15 22:21:07_coyote_last",
              #"GMU6_P_37_2020-07-28 00:11:29_coyote_last", "GMU6_P_37_2020-07-28 00:11:30_wolf_last", "GMU6_P_38_2020-08-28 21:03:42_wolf_last", "GMU6_P_38_2020-08-28 21:03:51_coyote_last",
              #"GMU6_P_56_2020-08-23 23:24:19_coyote_last", "GMU6_P_56_2020-08-23 23:24:15_bobcat_last",
              #"GMU6_P_56_2020-08-28 00:58:18_bobcat_last", "GMU6_P_56_2020-08-28 00:58:21_coyote_last", #"GMU6_P_58_2020-07-16 02:08:36_coyote_last", "GMU6_P_58_2020-07-16 02:08:47_wolf_last",
              #"GMU6_P_66_2020-07-31 22:16:15_bobcat_last", "GMU6_P_66_2020-07-31 22:16:18_coyote_last", # GMU6_P_66_2020-08-01 23:53:30_coyote_last & GMU6_P_66_2020-08-01 23:53:31_bobcat_first should be removed IF the unk deer classification in the middle of them is real
              #"GMU6_P_94_2020-07-10 05:24:20_bobcat_last", "GMU6_P_94_2020-07-10 05:24:20_coyote_last", "GMU6_P_94_2020-09-09 01:52:26_mountain_lion_last",
              #"GMU6_P_94_2020-09-09 01:53:41_bobcat_last", "GMU6_U_130_2020-08-21 23:57:41_wolf_last", "GMU6_U_130_2020-08-21 23:58:48_coyote_last", #
              #"GMU6_U_90_2020-09-04 13:27:34_mountain_lion_last", "GMU6_U_90_2020-09-04 13:27:50_bear_black_last"
              #"GMU10A_U_169_2020-08-07 22:40:17_wolf_last", "GMU10A_U_169_2020-08-07 22:40:19_coyote_last")
  rm_21s <- c(NA)#"GMU6_P_9_2021-09-02 21:57:37_bear_black_last", "GMU6_P_9_2021-09-02 21:57:38_wolf_last")

  #'  List all uniqueIDs of observations that need to switch the Det_type from
  #'  "first" to "last". These are truly the first image from each detection BUT
  #'  the first predator is still present when the second predator is detected so
  #'  need to artificially force these images to be the "last" image before the
  #'  second predator shows up
  change_to_last_20s <- c(NA)
  #                         "GMU10A_P_15_2020-08-07 00:27:15_bobcat_first", "GMU10A_P_23_2020-07-22 01:35:58_coyote_first",
  #                         "GMU10A_P_5_2020-09-15 07:27:50_coyote_first", "GMU10A_P_86_2020-07-27 01:47:12_coyote_first",
  #                         "GMU10A_P_86_2020-07-29 00:55:50_coyote_first", "GMU10A_P_86_2020-08-11 02:31:26_coyote_first",
  #                         "GMU10A_P_86_2020-08-11 02:39:57_coyote_first", "GMU10A_U_169_2020-08-07 22:36:03_coyote_first",
  #                         "GMU6_P_17_2020-07-01 22:35:41_mountain_lion_first", "GMU6_P_37_2020-07-28 00:11:27_wolf_first",
  #                         "GMU6_P_38_2020-08-28 21:03:00_coyote_first", "GMU6_P_56_2020-08-28 00:58:17_coyote_first",
  #                         "GMU6_P_58_2020-07-16 02:08:28_wolf_first", "GMU6_P_66_2020-07-31 22:16:13_coyote_first",
  #                         "GMU6_U_130_2020-08-21 23:54:00_coyote_first", "GMU6_U_90_2020-09-04 13:27:28_bear_black_first")
  change_to_last_21s <- c(NA) #"GMU6_P_9_2021-09-02 21:57:09_wolf_first")

  #'  Function to remove specific observations that screw up the last-first detection order
  remove_obs <- function(pred_pairs, rm_obs, new_last) {
    #'  Grab all rows that meet criteria based on uniqueID
    obs_to_remove <- pred_pairs[pred_pairs$uniqueID %in% rm_obs,]
    #'  Remove these rows from the larger df
    reduced_predator_pairs <- anti_join(pred_pairs, obs_to_remove)
    #'  Grab all rows that need Det_type changed and change them
    change_det_type <- reduced_predator_pairs[reduced_predator_pairs$uniqueID %in% new_last,]
    last_det_type <- change_det_type %>% mutate(Det_type = "last")
    #'  Remove observations with the original Det_type and replace with new
    reduced_predator_pairs <- anti_join(reduced_predator_pairs, change_det_type) %>%
      bind_rows(last_det_type) %>%
      #'  Arrange everything back in order
      arrange(NewLocationID, posix_date_time, File) #, desc(Det_type)
    return(reduced_predator_pairs)
  }
  predator_pairs_20s <- remove_obs(predator_pairs[[1]], rm_obs = rm_20s, new_last = change_to_last_20s)
  predator_pairs_21s <- remove_obs(predator_pairs[[3]], rm_obs = rm_21s, new_last = change_to_last_21s)
  #'  Being lazy and not doing this for wtr20 data right now

  #'  Re-list cleaned predator-pair data sets
  predator_pairs_thinned <- list(predator_pairs_20s, predator_pairs[[2]], predator_pairs_21s)


  #' #'  --------------------------------
  #' #####  Suspicious series of images  ####
  #' #'  --------------------------------
  #' #'  Pull out series of images to review based on the capture identifier
  #' #'  In these cases, 2 species were detected within ~1 min of each other, often
  #' #'  with spp2 showing up between images of spp1, so need to make sure these 
  #' #'  are not topological errors or miss-classifications
  #' caps_20s <- c(539, 1930, 2353, 4049, 5108, 8653, 8684, 8703, 8753, 18180, 18499, 
  #'               22218, 22494, 22533, 24840, 25000, 26015, 34178, 34234, 35783, 40143, 
  #'               11119, 13811)
  #' caps_21s <- c(44018)
  #' obs_to_check_20s <- predator_pairs[[1]][predator_pairs[[1]]$caps_new %in% caps_20s,]
  #' obs_to_check_21s <- predator_pairs[[3]][predator_pairs[[3]]$caps_new %in% caps_21s,]
  #' obs_to_check <- rbind(obs_to_check_20s, obs_to_check_21s) %>%
  #'   dplyr::select(-c(Category, Det_type, caps_new, cam, second_pred, first_pred, pred_pair, same_time, uniqueID))
  #' # write.csv(obs_to_check, "./Outputs/Tables/Images_to_double_check.csv")
  
  
  #'  ------------------------------------------
  ####  Final filtering to just focal pairings  ####
  #'  ------------------------------------------
  #'  Reduce to sequential detections of different predator species
  back_to_back_predators <- function(pred_pairs) {
    b2b_pred <- filter(pred_pairs, pred_pair == "Y") %>%
      arrange(NewLocationID, posix_date_time, Species)
    return(b2b_pred)
  }
  b2b_predators <- lapply(predator_pairs_thinned, back_to_back_predators) #predator_pairs
  
  #' #'  Repeat 2 coyote & 2 bobcat observations so sequence of detections has 
  #' #'  correct number of last/first observations
  #' GMU6_P_18 <- b2b_predators[[1]] %>% filter(NewLocationID == "GMU6_P_18")
  #' extra_obs1 <- GMU6_P_18 %>% filter(uniqueID == "GMU6_P_18_2020-08-15 22:19:16_coyote_first") %>%
  #'   mutate(Det_type = "last")
  #' GMU6_P_18_new <- bind_rows(GMU6_P_18, extra_obs1) %>%
  #'   arrange(posix_date_time, File)
  #' 
  #' GMU6_P_94 <- b2b_predators[[1]] %>% filter(NewLocationID == "GMU6_P_94")
  #' extra_obs2a <- GMU6_P_94 %>% filter(uniqueID == "GMU6_P_94_2020-07-10 05:24:19_coyote_first") %>%
  #'   mutate(Det_type = "last")
  #' extra_obs2b <- GMU6_P_94 %>% filter(uniqueID == "GMU6_P_94_2020-09-09 01:52:23_bobcat_first") %>%
  #'   mutate(Det_type = "last")
  #' GMU6_P_94_new <- bind_rows(GMU6_P_94, extra_obs2a, extra_obs2b) %>%
  #'   arrange(posix_date_time, File)
  #' 
  #' GMU6_P_56 <- b2b_predators[[1]] %>% filter(NewLocationID == "GMU6_P_56")
  #' extra_obs3 <- GMU6_P_56 %>% filter(uniqueID == "GMU6_P_56_2020-08-23 23:24:14_bobcat_first") %>%
  #'   mutate(Det_type = "last")
  #' GMU6_P_56_new <- bind_rows(GMU6_P_56, extra_obs3) %>%
  #'   arrange(posix_date_time, File)
  #' 
  #' #'  Remove these rows from the larger df and replace with updated dataframe
  #' b2b_predators[[1]] <- anti_join(b2b_predators[[1]], GMU6_P_94) %>%
  #'   bind_rows(GMU6_P_94_new) %>%
  #'   anti_join(GMU6_P_56) %>%
  #'   bind_rows(GMU6_P_56_new) %>%
  #'   anti_join(GMU6_P_18) %>%
  #'   bind_rows(GMU6_P_18_new) %>%
  #'   #'  Arrange everything back in order
  #'   arrange(NewLocationID, posix_date_time, File) #desc(Det_type)
  
  #'  Flag sequences that are out of order (e.g., last - last - first - first instead of last - first - last - first)
  #'  Usually occurs when next row is switching to a new camera location (ignore, 
  #'  not an issue for following code) OR there's only one image of a species & it's 
  #'    1) duplicated so that detection has a "first" and "last" image (ignore, 
  #'    not an issue for following code);
  #'    2) it's sandwiched between images of another species, likely due to species 
  #'    miss-classification, which MUST be corrected or can bias future TBD analyses
  flag_out_of_order_sequences <- function(dets) {
    whoops1 <- filter(dets, Det_type == lag(Det_type))
    whoops2 <- filter(dets, Det_type == lead(Det_type))
    whoops <- rbind(whoops1, whoops2) %>%
      distinct() %>%
      arrange(NewLocationID, posix_date_time, File)
    return(whoops)
  }
  double_check_order_20s <- flag_out_of_order_sequences(b2b_predators[[1]])
  double_check_order_21s <- flag_out_of_order_sequences(b2b_predators[[3]])
  
  #'  Double check everything looks alright
  tst <- b2b_predators[[1]]
  tst2 <- b2b_predators[[3]]
  
  #'  Flag back-to-back predator pairings where the gap between detections may
  #'  fall within problematic time period (when camera was inoperable)
  flag_artificial_b2b_pairs <- function(pred_pairs, seqprobs, start_date, end_date) {
    prob_date_range <- seqprobs %>%
      mutate(Date = as.Date(Date, format = "%d-%b-%Y")) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Filter to start and end of problem period
      group_by(NewLocationID) %>%
      filter(row_number()==1 | row_number()==n()) %>%
      ungroup()
    #'  Filter to just cameras where b2b predator observations occurred
    reduced_prob_dates <- prob_date_range[prob_date_range$NewLocationID %in% pred_pairs$NewLocationID,]
    b2b_prob_cams <- pred_pairs[pred_pairs$NewLocationID %in% prob_date_range$NewLocationID,]
    #'  Print how many cameras with problem date ranges are also in the b2b predator data set?
    print(nrow(reduced_prob_dates))
    print(nrow(b2b_prob_cams))
    
    #'  Do dates in reduced_prob_dates fall between dates in b2b_prob_cams for each predator pairing at a given site???
    
    return(reduced_prob_dates)
  }
  b2b_predators_20s <- flag_artificial_b2b_pairs(b2b_predators[[1]], eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15")
  b2b_predators_21s <- flag_artificial_b2b_pairs(b2b_predators[[3]], eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15")
  
  
  #'  ---------------------------------------------
  ####  Calculate times between detection events   ####
  #'  ---------------------------------------------
  #'  Function to calculate time between detection events of two focal species
  #'  Data structured so only last image of spp1 and first image of spp2 per
  #'  detection event are included in data frame.
  tbd <- function(detection_data, det_type, unittime) {
    detection_data <- detection_data %>%
      group_by(NewLocationID) %>%
      mutate(TimeSinceLastDet = difftime(posix_date_time, lag(posix_date_time), units = unittime),
             TimeSinceLastDet = ifelse(Det_type == "last", 0, TimeSinceLastDet),
             MinutesSinceLastDet = round(TimeSinceLastDet, 4),
             HoursSinceLastDet = round((TimeSinceLastDet/60), 2),
             DaysSinceLastDet = round((TimeSinceLastDet/1440), 2),
             Previous_Spp = lag(Species),
             Predator_pair = paste0(Previous_Spp, "-", Species)) %>%
      ungroup() %>%
      # relocate(Previous_Spp, .before = Species) %>%
      # relocate(Predator_pair, .after = Species) %>%
      filter(Det_type == "first") %>%
      dplyr::select(-c(TriggerMode, OpState, Category, Count, caps_new, cam, same_time, pred_pair, uniqueID)) #Det_type, second_pred, first_pred, 
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_pred_pairs <- lapply(b2b_predators, tbd, det_type = "last", unittime = "min")
  
  #'  Add a column for GMU
  gmu_col <- function(tbd) {
    tbd <- mutate(tbd, GMU = sub("_.*", "", NewLocationID))
    return(tbd)
  }
  tbd_pred_pairs <- lapply(tbd_pred_pairs, gmu_col)
  
  #'  Add column for year
  tbd_pred_pairs[[1]]$Year <- "Smr20"
  tbd_pred_pairs[[3]]$Year <- "Smr21"
  tbd_pred_pairs_all <- rbind(tbd_pred_pairs[[1]], tbd_pred_pairs[[3]])
  
  save(tbd_pred_pairs_all, file = paste0("./Data/Time_btwn_Detections/TBD_all_predator_pairs_", Sys.Date(), ".RData"))
  
  
  #'  ----------------------------------------
  ####  Summary stats and data visualization  ####
  #'  ----------------------------------------
  #'  Total average
  mean(tbd_pred_pairs_all$MinutesSinceLastDet, na.rm = TRUE); sd(tbd_pred_pairs_all$MinutesSinceLastDet, na.rm = TRUE)
  mean(tbd_pred_pairs_all$HoursSinceLastDet, na.rm = TRUE); sd(tbd_pred_pairs_all$HoursSinceLastDet, na.rm = TRUE)
  
  #'  Average by species pairing
  avg_tbd <- tbd_pred_pairs_all %>%
    group_by(Predator_pair) %>%
    summarise(total_obs = n(),
              mean_tbd = mean(HoursSinceLastDet),
              sd_tbd = sd(HoursSinceLastDet)) %>%
    ungroup()
  
  hist(tbd_pred_pairs_all$HoursSinceLastDet, breaks = 50, main = "Elapsed time between sequential detections of predators", xlab = "Hours between detection events")
  hist(tbd_pred_pairs_all$DaysSinceLastDet, breaks = 50, main = "Elapsed time between sequential detections of predators", xlab = "Days between detection events")
  
  #'  Split by species-pairs
  wolfbear <- filter(tbd_pred_pairs_all, Predator_pair == "wolf-bear_black" | Predator_pair == "bear_black-wolf")
  wolflion <- filter(tbd_pred_pairs_all, Predator_pair == "wolf-mountain_lion" | Predator_pair == "mountain_lion-wolf")
  wolfcoy <- filter(tbd_pred_pairs_all, Predator_pair == "wolf-coyote" | Predator_pair == "coyote-wolf")
  lionbear <- filter(tbd_pred_pairs_all, Predator_pair == "mountain_lion-bear_black" | Predator_pair == "bear_black-mountain_lion")
  lionbob <- filter(tbd_pred_pairs_all, Predator_pair == "mountain_lion-bobcat" | Predator_pair == "bobcat-mountain_lion")
  coybob <- filter(tbd_pred_pairs_all, Predator_pair == "coyote-bobcat" | Predator_pair == "bobcat-coyote")
  
  #'  Focus in on specific pairings
  focal_pairs <- rbind(wolfbear, wolfcoy, wolflion, lionbear, lionbob, coybob) %>%
    mutate(Predator_pair = gsub("_", " ", Predator_pair))
  
  #'  Summary stats of just focal species pairings
  avg_tbd_focal_pairs <- focal_pairs %>%
    group_by(Predator_pair) %>%
    summarise(total_obs = n(),
              mean_tbd_hr = round(mean(HoursSinceLastDet), 2),
              sd = round(sd(HoursSinceLastDet), 2),
              se = round(sd(HoursSinceLastDet)/sqrt(total_obs), 2)) %>%
    ungroup() %>%
    arrange(desc(total_obs)) %>%
    mutate(Predator_pair = gsub("_", " ", Predator_pair))
  # write.csv(avg_tbd_focal_pairs, "./Outputs/Tables/Summary_stats_TBD.csv")
  
  #'  Visualize tbd data for data sets with >50 observations
  coy_bob_his <- ggplot(coybob, aes(x=HoursSinceLastDet, fill=Predator_pair)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 5) +
    scale_fill_manual(values=c("#69b3a2", "#404080")) 
  wolf_coy_his <- ggplot(wolfcoy, aes(x=HoursSinceLastDet, fill=Predator_pair)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 5) +
    scale_fill_manual(values=c("#69b3a2", "#404080"))
  plot(coy_bob_his)
  plot(wolf_coy_his)
  
  #'  Add summary data to each observation
  focal_pairs <- full_join(focal_pairs, avg_tbd_focal_pairs, by = "Predator_pair") 
  
  #'  Boxplots of all data, organized by sample size
  pred_pair_box <- ggplot(focal_pairs, aes(x = reorder(Predator_pair, total_obs), y=HoursSinceLastDet, fill=Predator_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(data = avg_tbd_focal_pairs,
              aes(Predator_pair, Inf, label = paste("n =",total_obs)), vjust = 8) +
    ggtitle("Boxplots summarizing elapsed time between detections of predators") +
    xlab("Predator pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First - Second Predator")
  plot(pred_pair_box)
  
  
  