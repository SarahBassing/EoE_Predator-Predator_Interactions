  #'  ---------------------------------------
  #'  Calculate time between detections 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ---------------------------------------
  #'  Filter detection data down to back-to-back pairs for focal species as well 
  #'  as non-target species followed by a focal species.
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
  newnewSppID <- read.csv("./Data/IDFG camera data/More questionable IDs_ST_SBB.csv") %>% 
    mutate(posix_date_time = gsub("<a0>", " ", posix_date_time),
           posix_date_time = as.POSIXct(posix_date_time, format="%Y-%m-%d %H:%M", tz="UTC"),
           Year = as.numeric(format(posix_date_time, "%Y"))) %>%
    dplyr::select(c("CamID", "NewLocationID", "File", "Location_Relative_Project", "posix_date_time", "Species", "NewSpecies", "Year"))
  newSppID <- rbind(newSppID, newnewSppID) %>% distinct() 
  #'  Filter by year and only observations with misidentified species
  newSppID_20s <- filter(newSppID, Year == "2020") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  newSppID_21s <- filter(newSppID, Year == "2021") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  #'  Snag unique ID for each image
  change_sppID_20s <- as.vector(newSppID_20s$Location_Relative_Project)
  change_sppID_21s <- as.vector(newSppID_21s$Location_Relative_Project)
  
  #'  Function to correct species misclassifications
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
  #'  5. Double check inoperable date periods don't overlap back-to-back detections.
  #'  6. Calculate time-between-detections of different species
  #'  7. Thin image set to just detections of different predator species
  #'  --------------------------------------
  
  #'  Group multiple detection events of same category (but of different species)
  #'  when they occur sequentially, then reduce groups of "other" category to a
  #'  single observation and reduce predator detections to first and last image.
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
    capdata <- cbind(as.data.frame(dat), caps_new) %>%
      arrange(NewLocationID, posix_date_time)
    
    #'  Remove all extra detections when multiple detections of same category occur in a row
    predspp <- capdata[capdata$Category == "Predator",] %>%
      group_by(Species, caps_new) %>%
      slice(1, n()) %>%
      ungroup()
    lastother <- capdata[capdata$Category == "Other" & lead(capdata$Category == "Predator"),]
    # lasteverythingelse <- capdata[capdata$Category != "Predator",] %>%
    #   group_by(caps_new) %>%
    #   slice_tail() %>%
    #   ungroup()
    #'  Combine into final data set
    dets <- rbind(predspp, lastother) %>% #lasteverythingelse
      arrange(NewLocationID, posix_date_time) %>%
      #'  Create a unique ID so I can more easily remove specific observations
      mutate(uniqueID = paste0(NewLocationID, "_", posix_date_time, "_", Species, "_", Det_type))
    return(dets)
  }
  full_predator_sequence <- lapply(firstlast_img, thin_dat_by_category) 
  
  
  #'  ------------------------------------------
  ####  Remove problem obs, dates, and cameras  ####
  #'  ------------------------------------------
  #####  1) Remove instances of known species miss-classification  #####
  #'  List all uniqueIDs of observations that need to be removed in a giant vector
  #'  CURRENT LIST is just to make my life easier for now - these are incorrect
  #'  species classifications that need to be corrected. Once corrected I can revert to rm_20s <- c(NA)
  # rm_20s <- c("GMU6_P_14_2020-08-28 08:58:35_coyote_last", "GMU6_P_14_2020-08-28 09:00:02_bear_black_last")
  rm_20s <- c(NA)
  rm_21s <- c(NA)
  
  #'  Function to remove specific observations that screw up the last-first detection order
  remove_obs <- function(spp_pairs, rm_obs) {
    #'  Grab all rows that meet criteria based on uniqueID
    obs_to_remove <- spp_pairs[spp_pairs$uniqueID %in% rm_obs,]
    #'  Remove these rows from the larger df
    reduced_spp_pairs <- anti_join(spp_pairs, obs_to_remove) %>%
      #'  Arrange everything back in order
      arrange(NewLocationID, posix_date_time, File) 
    return(reduced_spp_pairs)
  }
  spp_pairs_20s <- remove_obs(full_predator_sequence[[1]], rm_obs = rm_20s)
  spp_pairs_21s <- remove_obs(full_predator_sequence[[3]], rm_obs = rm_21s)
  #'  Being lazy and not doing this for wtr20 data right now
  
  #'  Re-list cleaned species-pair data sets
  predator_pairs_thinned <- list(spp_pairs_20s, spp_pairs_21s)
  
  #####  2) Remove observations that fall within inoperable time periods  ####
  #'  Review cameras with problem time periods and the dates of species to predator 
  #'  detections to make sure elapsed time between detections does not encompass
  #'  problematic time periods (don't want inoperable cameras to bias TBD data)
  problem_dates_and_b2b_dets <- function(spp_pairs, seqprobs, start_date, end_date) {
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
    reduced_prob_dates <- prob_date_range[prob_date_range$NewLocationID %in% spp_pairs$NewLocationID,]
    #'  Identify cameras that had b2b predator detections but were also inoperable for some period of time
    b2b_prob_cams <- spp_pairs[spp_pairs$NewLocationID %in% prob_date_range$NewLocationID,]
    
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
    
    #'  Create data frame with start and end of b2b data that potentially overlap problem dates
    b2b <- b2b_prob_cams %>%
      group_by(NewLocationID) %>%
      mutate(Spp1_detection = posix_date_time,
             Spp2_detection = lead(posix_date_time),
             Previous_predator = lead(Species)) %>%
      ungroup() %>%
      filter(!is.na(Spp2_detection)) %>%
      dplyr::select(c("NewLocationID", "Species", "Category", "Spp1_detection", "Spp2_detection", "Previous_predator"))
    
    #'  Flag which observations fall within problematic date range. These snuck in 
    #'  and should be removed from data set.
    b2b_prob_overlap <- b2b %>%
      left_join(prob_dates, by = "NewLocationID") %>%
      mutate(sneaky_mf = ifelse(Spp1_detection >= Prob_start & Spp1_detection <= Prob_end, 1, NA),
             sneaky_mf = ifelse(Spp2_detection >= Prob_start & Spp2_detection <= Prob_end, 1, sneaky_mf)) %>%
      filter(!is.na(sneaky_mf))
    print(b2b_prob_overlap)
    
    #'  Create data frame of just images to remove based on those listed in tbd_prob_overlap
    problem_cam <- spp_pairs[spp_pairs$NewLocationID %in% b2b_prob_overlap$NewLocationID,]
    problem_spp1_obs <- problem_cam[problem_cam$posix_date_time %in% b2b_prob_overlap$Spp1_detection,]
    problem_spp2_obs <- problem_cam[problem_cam$posix_date_time %in% b2b_prob_overlap$Spp2_detection,]
    problem_obs <- rbind(problem_spp1_obs, problem_spp2_obs) %>% arrange(posix_date_time)
    #'  Double check correct images are being removed (posix_date_time should match
    #'  Spp1_detection and Spp2_detection columns from tbd_prob_overlap)
    print(problem_obs)
    
    #'  Remove observations from larger time-to-detection data frame if they fall
    #'  within a problematic time period at a camera site
    thinned_b2b_dat <- anti_join(spp_pairs, problem_obs)
    
    return(thinned_b2b_dat)
  }
  b2b_dets_20s <- problem_dates_and_b2b_dets(predator_pairs_thinned[[1]], eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15") 
  b2b_dets_21s <- problem_dates_and_b2b_dets(predator_pairs_thinned[[2]], eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15") 
  
  b2b_dets <- list(b2b_dets_20s, b2b_dets_21s)
  #'  Double check everything looks alright
  tst <- b2b_dets[[1]]
  tst2 <- b2b_dets[[2]]
  
  
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
             Focal_predator = Species,
             Previous_Spp = lag(Species),
             Species_pair = paste0(Previous_Spp, "-", Focal_predator)) %>%
      ungroup() %>%
      filter(Det_type == "first") %>%
      filter(!is.na(TimeSinceLastDet)) %>%
      dplyr::select(-c(TriggerMode, OpState, Species, Category, Count, Det_type, caps_new, uniqueID)) %>%
      relocate(Focal_predator, .after = posix_date_time) %>%
      relocate(Previous_Spp, .after = Focal_predator) %>%
      relocate(Species_pair, .after = Previous_Spp)
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  #'  Note: there should be NO negative values! If there are negative values this
  #'  means the script is calculating times between detections across camera sites
  tbd_spp_pairs <- lapply(b2b_dets, tbd, det_type = "last", unittime = "min")
  #'  Double check everything looks alright
  tst3 <- tbd_spp_pairs[[1]]
  tst4 <- tbd_spp_pairs[[2]]
  
  #'  Add a column for GMU
  gmu_col <- function(tbd) {
    tbd <- mutate(tbd, GMU = sub("_.*", "", NewLocationID))
    return(tbd)
  }
  tbd_spp_pairs <- lapply(tbd_spp_pairs, gmu_col)
  
  #'  Add column for year
  tbd_spp_pairs[[1]]$Year <- "Smr20"
  tbd_spp_pairs[[2]]$Year <- "Smr21"
  tbd_spp_pairs_all <- rbind(tbd_spp_pairs[[1]], tbd_spp_pairs[[2]])
  
  save(tbd_spp_pairs_all, file = paste0("./Data/Time_btwn_Detections/TBD_all_spp-pred_pairs_", Sys.Date(), ".RData"))
  
  
  #'  ----------------------------------------
  ####  Summary stats and data visualization  ####
  #'  ----------------------------------------
  #'  Total average
  mean(tbd_spp_pairs_all$MinutesSinceLastDet, na.rm = TRUE); sd(tbd_spp_pairs_all$MinutesSinceLastDet, na.rm = TRUE)
  mean(tbd_spp_pairs_all$HoursSinceLastDet, na.rm = TRUE); sd(tbd_spp_pairs_all$HoursSinceLastDet, na.rm = TRUE)
  
  #'  Average by species pairing
  avg_tbd <- tbd_spp_pairs_all %>%
    group_by(Species_pair) %>%
    summarise(total_obs = n(),
              mean_tbd = mean(HoursSinceLastDet),
              sd_tbd = sd(HoursSinceLastDet)) %>%
    ungroup()
  
  hist(tbd_spp_pairs_all$HoursSinceLastDet, breaks = 50, main = "Elapsed time between sequential detections of Spp1 and a predator", xlab = "Hours between detection events")
  hist(tbd_spp_pairs_all$DaysSinceLastDet, breaks = 50, main = "Elapsed time between sequential detections of Spp1 and a predator", xlab = "Days between detection events")
  
  #'  Split by species-pairs
  wolfbear <- filter(tbd_spp_pairs_all, Species_pair == "wolf-bear_black" | Species_pair == "bear_black-wolf")
  wolflion <- filter(tbd_spp_pairs_all, Species_pair == "wolf-mountain_lion" | Species_pair == "mountain_lion-wolf")
  wolfcoy <- filter(tbd_spp_pairs_all, Species_pair == "wolf-coyote" | Species_pair == "coyote-wolf")
  lionbear <- filter(tbd_spp_pairs_all, Species_pair == "mountain_lion-bear_black" | Species_pair == "bear_black-mountain_lion")
  lionbob <- filter(tbd_spp_pairs_all, Species_pair == "mountain_lion-bobcat" | Species_pair == "bobcat-mountain_lion")
  coybob <- filter(tbd_spp_pairs_all, Species_pair == "coyote-bobcat" | Species_pair == "bobcat-coyote")
  
  #'  Focus in on specific pairings
  focal_pairs <- rbind(wolfbear, wolfcoy, wolflion, lionbear, lionbob, coybob) %>%
    mutate(Species_pair = gsub("_", " ", Species_pair))
  
  #'  Summary stats of just focal species pairings
  avg_tbd_focal_pairs <- focal_pairs %>%
    group_by(Species_pair) %>%
    summarise(total_obs = n(),
              mean_tbd_hr = round(mean(HoursSinceLastDet), 2),
              sd = round(sd(HoursSinceLastDet), 2),
              se = round(sd(HoursSinceLastDet)/sqrt(total_obs), 2)) %>%
    ungroup() %>%
    arrange(desc(total_obs)) %>%
    mutate(Species_pair = gsub("_", " ", Species_pair))
  # write.csv(avg_tbd_focal_pairs, "./Outputs/Tables/Summary_stats_TBD.csv")
  
  #'  Visualize tbd data for data sets with >50 observations
  coy_bob_his <- ggplot(coybob, aes(x=HoursSinceLastDet, fill=Species_pair)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 5) +
    scale_fill_manual(values=c("#69b3a2", "#404080")) 
  wolf_coy_his <- ggplot(wolfcoy, aes(x=HoursSinceLastDet, fill=Species_pair)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 5) +
    scale_fill_manual(values=c("#69b3a2", "#404080"))
  plot(coy_bob_his)
  plot(wolf_coy_his)
  
  #'  Add summary data to each observation
  focal_pairs <- full_join(focal_pairs, avg_tbd_focal_pairs, by = "Species_pair") 
  
  #'  Boxplots of all data, organized by sample size
  pred_pair_box <- ggplot(focal_pairs, aes(x = reorder(Species_pair, total_obs), y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(data = avg_tbd_focal_pairs,
              aes(Species_pair, Inf, label = paste("n =",total_obs)), vjust = 8) +
    ggtitle("Boxplots summarizing elapsed time between detections of predators") +
    xlab("Predator pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First - Second Predator")
  plot(pred_pair_box)
  
  #'  Focal species boxplots
  bear_response <- tbd_spp_pairs_all[tbd_spp_pairs_all$Focal_predator == "bear_black",] %>%
    ggplot(aes(x = Species_pair, y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Boxplots summarizing elapsed time between detections of any species and black bear") +
    xlab("Species pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First Species - Bear")
  bob_response <- tbd_spp_pairs_all[tbd_spp_pairs_all$Focal_predator == "bobcat",] %>%
    ggplot(aes(x = Species_pair, y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Boxplots summarizing elapsed time between detections of any species and bobcat") +
    xlab("Species pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First Species - Bobcat")
  coy_response <- tbd_spp_pairs_all[tbd_spp_pairs_all$Focal_predator == "coyote",] %>%
    ggplot(aes(x = Species_pair, y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Boxplots summarizing elapsed time between detections of any species and coyote") +
    xlab("Species pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First Species - Coyote")
  lion_response <- tbd_spp_pairs_all[tbd_spp_pairs_all$Focal_predator == "mountain_lion",] %>%
    ggplot(aes(x = Species_pair, y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Boxplots summarizing elapsed time between detections of any species and lion") +
    xlab("Species pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First Species - Lion")
  wolf_response <- tbd_spp_pairs_all[tbd_spp_pairs_all$Focal_predator == "wolf",] %>%
    ggplot(aes(x = Species_pair, y=HoursSinceLastDet, fill=Species_pair)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Boxplots summarizing elapsed time between detections of any species and wolf") +
    xlab("Species pairings") + ylab("Hours between sequential detections") +
    labs(fill = "First Species - Wolf")
  
  bear_response
  bob_response
  coy_response
  lion_response
  wolf_response
  
  #'  Table of non-target mean tbd
  lago <- filter(avg_tbd, grepl("rabbit_hare", Species_pair))
  elk <- filter(avg_tbd, grepl("elk", Species_pair))
  moose <- filter(avg_tbd, grepl("moose", Species_pair))
  wtd <- filter(avg_tbd, grepl("whitetaileddeer", Species_pair))
  nt_pairs <- rbind(lago, elk, moose, wtd) %>%
    mutate(Species_pair = gsub("_", " ", Species_pair)) %>%
    arrange(Species_pair)
  