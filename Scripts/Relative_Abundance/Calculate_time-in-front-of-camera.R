  #'  ----------------------------------------
  #'  Calculate time-in-front-of-camera (TIFC)
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  September 2023
  #'  ----------------------------------------
  #'  Calculate time-in-front-of-camera for each species following methods 
  #'  described by Huggard et al. 2018, Warbington & Boyce 2020, and Becker et al. 2022.
  #'  Code modified from Becker et al. 2018 repository: https://github.com/mabecker89/tifc-method/tree/main
  #'  Last, organize field-of-view measurements at each cameras site based on 
  #'  IDFG's various methods of measuring the detection zone.
  #'  
  #'  Sequential probimgs generated in Detection_data_cleaning.R (Data_Formatting folder)
  #'  Sampling effort data generated in Detection_histories_for_occmod.R (MultiSpp_OccMod folder)
  #'  Leave probabilities taken directly from Becker et al. (2022)
  #'  ----------------------------------------
  
  #'  Load libraries
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  
  #'  -------------------------------
  ####  Load and format input data  ####
  #'  -------------------------------
  #'  Load detection data
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_2023-09-26.RData") 
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
  #'  Filter detection data to time period of interest and remove problem images
  thin_detections <- function(dets, seqprobs, start_date, end_date) {
    #'  Remove images with known problems (e.g., viewshed obscured)
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    
    #'  Filter images
    clean_dets <- skinny_dets %>%
      filter(Species != "none") %>%
      filter(Vehicle != "TRUE") %>%
      filter(OpState != "maintenance") %>%
      #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
      mutate(Count = ifelse(Count == 0, 1, Count),
             Count = ifelse(is.na(Count), 1, Count),
             Date = as.Date(Date, format = "%Y-%m-%d")) %>% #format = "%d-%b-%Y"
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove single mislabeled image from 2022 data (will not affect TIFC calculations)
      filter(NewLocationID != "GMU10A_U_50" | Date != "2022-08-27" | Species != "human") %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      arrange(NewLocationID, posix_date_time) %>%
      dplyr::select(c("NewLocationID", "posix_date_time", "Species", "Count")) %>%
      mutate(gap_class = NA)
    names(clean_dets) <- c("location", "date_detected", "common_name", "number_individuals", "gap_class")
    
    return(clean_dets)
  }
  df_all_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15") 
  df_all_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15")
  df_all_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-07-01", end_date = "2022-09-15") 
  
  #'  Update species names for "gap groups" based on Becker et al. (2022) classifications
  spp_groups_20s <- as.data.frame(unique(df_all_20s$common_name)); names(spp_groups_20s) <- "common_name"
  spp_groups_21s <- as.data.frame(unique(df_all_21s$common_name)); names(spp_groups_21s) <- "common_name"
  spp_groups_22s <- as.data.frame(unique(df_all_22s$common_name)); names(spp_groups_22s) <- "common_name"
  spp_groups_20s <- arrange(spp_groups_20s, common_name)
  spp_groups_21s <- arrange(spp_groups_21s, common_name)
  spp_groups_22s <- arrange(spp_groups_22s, common_name)
  spp_groups_20s$gap_group <- c("Small carnivores", "Bears", "Other", "Small carnivores", 
                                "Small carnivores", "Livestock", "Canids cougar", "Most ungulates", "Canids cougar", 
                                "Most ungulates", "Small carnivores", "Canids cougar", "Most ungulates", "Livestock", 
                                "Human", "Small carnivores", "Moose", "Canids cougar", "Most ungulates", "Small carnivores", 
                                "Small mammals", "Small mammals", "Small carnivores", "Small mammals", "Small mammals", 
                                "Small mammals", "Other", "Other", "Small carnivores", "Most ungulates", 
                                "Canids cougar", "Small carnivores")
  spp_groups_21s$gap_group <- c("Small carnivores", "Small mammals", "Bears", "Bears", "Other", "Small carnivores", 
                                "Small carnivores", "Livestock", "Canids cougar", "Most ungulates", "Canids cougar", 
                                "Most ungulates", "Small carnivores", "Canids cougar", "Livestock", "Most ungulates", 
                                "Livestock", "Human", "Other", "Small carnivores", "Small mammals", "Small carnivores", 
                                "Moose", "Canids cougar", "Most ungulates", "Small mammals", "Small mammals", 
                                "Small carnivores", "Small mammals", "Small mammals", "Small mammals", "Other", "Other", 
                                "Small carnivores", "Most ungulates", "Canids cougar")
  spp_groups_22s$gap_group <- c("Small carnivores", "Small mammals", "Bears", "Bears", "Other", "Small carnivores", 
                                "Small carnivores", "Livestock", "Other", "Canids cougar", "Most ungulates", "Canids cougar", 
                                "Most ungulates", "Small carnivores", "Canids cougar", "Livestock", "Most ungulates", 
                                "Small mammals", "Livestock", "Human", "Small carnivores", "Small mammals", "Moose", 
                                "Canids cougar", "Most ungulates", "Small mammals", "Small mammals", "Small carnivores",
                                "Small mammals", "Small mammals", "Small mammals", "Other", "Other", 
                                "Canids cougar", "Small carnivores", "Most ungulates", "Canids cougar")
  df_gap_groups_20s <- spp_groups_20s; df_gap_groups_21s <- spp_groups_21s; df_gap_groups_22s <- spp_groups_22s
  
  #'  Probabilistic gaps
  df_leave_prob_pred <- read_csv("./Data/Becker et al. data/gap-leave-prob_predictions.csv")
  
  #'  Sampling effort data (number of days cameras were operational)
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData") 
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe22s.RData") 
  
  df_tbd_20s <- effort_20s %>% dplyr::select(c("NewLocationID", "ndays")) 
  names(df_tbd_20s) <- c("location", "total_days")
  df_tbd_21s <- effort_21s %>% dplyr::select(c("NewLocationID", "ndays")) 
  names(df_tbd_21s) <- c("location", "total_days")
  df_tbd_22s <- effort_22s %>% dplyr::select(c("NewLocationID", "ndays")) 
  names(df_tbd_22s) <- c("location", "total_days")
  #'  Missing a camera in 2021 sampling effort data for some reason
  add_missing_cam <- c("GMU1_U_153", "18")
  df_tbd_21s <- rbind(df_tbd_21s, add_missing_cam) %>%
    mutate(total_days = as.numeric(total_days)) %>%
    arrange(location)
  
  #'  Seasonal start/end dates (julian day)
  july1.20 <- as.Date("2020-07-01", format = "%Y-%m-%d")
  sept15.20 <- as.Date("2020-09-15", format = "%Y-%m-%d")
  summer20.start.j <- format(july1.20, "%j") 
  summer20.end.j <- format(sept15.20, "%j") 
  
  july1.21 <- as.Date("2021-07-01", format = "%Y-%m-%d")
  sept15.21 <- as.Date("2021-09-15", format = "%Y-%m-%d")
  summer21.start.j <- format(july1.21, "%j") 
  summer21.end.j <- format(sept15.21, "%j") 
  
  july1.22 <- as.Date("2022-07-01", format = "%Y-%m-%d")
  sept15.22 <- as.Date("2022-09-15", format = "%Y-%m-%d")
  summer22.start.j <- format(july1.22, "%j") 
  summer22.end.j <- format(sept15.22, "%j") 
  
  
  #'  ------------------
  ####  Calculate TIFC  ####
  #'  ------------------
  ##### Step 1. Identify Detection Series  #####
  #'  ------------------------------------------
  img_series <- function(df_all, df_gap_groups) {
    df_series <- df_all %>%
      #' #' Join gap class (if have gap data)
      #' left_join(df_gap, by = c("location", "date_detected", "common_name")) %>%
      #'  Order observations
      arrange(location, date_detected, common_name) %>%
      #'  Identify series and gaps requiring probabilistic time assignment
      mutate(series_num = 0,
             #'  Lagged time stamp
             date_detected_previous = lag(date_detected),
             #'  Lead time stamp
             date_detected_next = lead(date_detected),
             #'  Calculate difference in time between ordered images
             diff_time_previous = as.numeric(date_detected - date_detected_previous),
             diff_time_next = as.numeric(abs(date_detected - date_detected_next)),
             #'  Lagged species
             common_name_previous = lag(common_name),
             #'  Was it a different species?
             diff_sp = ifelse(common_name != common_name_previous, TRUE, FALSE),
             #'  Lagged deployment
             location_previous = lag(location),
             #'  Was it a different deployment?
             diff_location = ifelse(location != location_previous, TRUE, FALSE),
             #'  Flag gaps that will need checking
             gap_check = ifelse(diff_location == FALSE & diff_sp == FALSE & (diff_time_previous <= 120 & diff_time_previous >= 20), 1, 0),
             #'  Lagged gap class (Not actually necessary since currently no gap class data) 
             gap_class_previous = replace_na(lag(gap_class)), #replace_na(lag(gap_class), ""),
             gap_class_previous = ifelse(is.na(gap_class_previous), "", gap_class_previous),
             #'  Identify new series, based on being different deployment, species, greater than 120 seconds, and approp gaps
             diff_series = ifelse(diff_location == TRUE | diff_sp == TRUE | diff_time_previous > 120 | (gap_class_previous == "L" | gap_class_previous == "N"), 1, 0),
             #'  Number series
             series_num = c(0, cumsum(diff_series[-1])),
             #'  Flag gaps that require probabilistic time assignment
             gap_prob = replace_na(ifelse(gap_check == 1 & (gap_class_previous == "" | gap_class_previous == "U"), 1, 0), 0)) %>%
      group_by(series_num) %>%
      mutate(diff_time_previous = ifelse(row_number() == 1, 0, diff_time_previous),
             diff_time_next = ifelse(row_number() == n(), 0, diff_time_next)) %>%
      ungroup() %>%
      #'  Join gap group classification
      left_join(df_gap_groups, by = "common_name") %>%
      #'  Drop detections from certain groups
      filter(gap_group != "Human" & gap_group != "Other" & gap_group != "Livestock") %>%
      #'  Join gap leaving predictions
      left_join(df_leave_prob_pred, by = c("gap_group", "diff_time_previous" = "diff_time")) %>%
      #'  Adjust time difference between ordered images that require probabilistic time assignment
      mutate(pred = replace_na(pred, 1),
             diff_time_previous_adj = ifelse(gap_prob == 1, diff_time_previous * (1 - pred), diff_time_previous),
             diff_time_next_adj = ifelse(lead(gap_prob == 1), diff_time_next * (1 - lead(pred)), diff_time_next),
             #'  On the rare occasion a single image of a species is last image at camera, 
             #'  force time_diff_next_adj = 0 so not left as NA (required for later calculations)
             #'  e.g., GMU6_U_97 Summer 2020 has 1 img of coyote & last photo taken by camera
             diff_time_next_adj = ifelse(is.na(diff_time_next_adj), 0, diff_time_next_adj))
    return(df_series)
  }
  df_series_20s <- img_series(df_all_20s, df_gap_groups = df_gap_groups_20s)
  df_series_21s <- img_series(df_all_21s, df_gap_groups = df_gap_groups_21s)
  df_series_22s <- img_series(df_all_22s, df_gap_groups = df_gap_groups_22s)
  
  
  ##### Step 2. Calculate time between photos (tbp), by species #####
  #'  ---------------------------------------------------------------
  tbp <- function(df_series) {
    df_tbp <- df_series %>%
      mutate(series_num_previous = lag(series_num)) %>%
      #'  Remove first image from each series
      filter(series_num == series_num_previous) %>%
      group_by(common_name) %>%
      #'  Calculate average tbp and number of images from each species
      summarise(tbp = mean(diff_time_previous),
                sample_size = n())
    return(df_tbp)
  }
  df_tbp_20s <- tbp(df_series_20s)
  df_tbp_21s <- tbp(df_series_21s)
  df_tbp_22s <- tbp(df_series_22s)
  
  
  #####  Step 3. Calculate total time in front of the camera, by series (tts = Total Time by Series)  #####
  #'  -----------------------------------------------------------------------------------------------------
  tts <- function(df_series, df_tbp) {
    df_tts <- df_series %>%
      left_join(df_tbp, by = "common_name") %>%
      group_by(series_num) %>%
      #'  Check whether the image was first or last in a series
      mutate(bookend = ifelse(row_number() == 1 | row_number() == n(), 1, 0),
             #'  Calculate time for each individual image
             image_time = ifelse(bookend == 1,
                                 ((diff_time_previous_adj + diff_time_next_adj) / 2) + (tbp / 2),
                                 (diff_time_previous_adj + diff_time_next_adj) / 2),
             #'   Multiply image time by the number of animals present
             image_time_ni = image_time * number_individuals) %>%
      #'  Group by common name as well to add it as a variable to output
      group_by(common_name, .add = TRUE) %>%
      #'  Calculate total time and number of images for each series
      summarise(n_images = n(),
                series_total_time = sum(image_time_ni)) %>%
      ungroup()
    return(df_tts)
  }
  df_tts_20s <- tts(df_series_20s, df_tbp_20s)
  df_tts_21s <- tts(df_series_21s, df_tbp_21s)
  df_tts_22s <- tts(df_series_22s, df_tbp_22s)
  
  
  #####  Step 4. Calculate total time in front of camera, by deployment and species (tt = total time)  #####
  #'  ------------------------------------------------------------------------------------------------------
  tt <- function(df_series, df_tts, df_tbd, summer.start.j, summer.end.j) {
    df_tt <- df_series %>%
      group_by(series_num) %>%
      arrange(date_detected, .by_group = TRUE) %>%
      filter(row_number() == 1) %>%
      left_join(df_tts, by = c("series_num", "common_name")) %>%
      dplyr::select(location, date_detected, common_name, series_num, series_total_time) %>%
      ungroup() %>%
      mutate(julian = as.numeric(format(date_detected, "%j")),
             season = ifelse(julian >= summer.start.j & julian <= summer.end.j, "summer", "winter")) %>%
      mutate_at(c("location", "common_name", "season"), factor) %>%
      group_by(location, common_name, season, .drop = FALSE) %>%
      summarise(total_duration = sum(series_total_time)) %>%
      ungroup() %>%
      mutate_if(is.factor, as.character) %>%
      left_join(df_tbd, by = "location")
    return(df_tt)
  }
  df_tt_20s <- tt(df_series_20s, df_tts_20s, df_tbd_20s, summer20.start.j, summer20.end.j)
  df_tt_21s <- tt(df_series_21s, df_tts_21s, df_tbd_21s, summer21.start.j, summer21.end.j)
  df_tt_22s <- tt(df_series_22s, df_tts_22s, df_tbd_22s, summer22.start.j, summer22.end.j)
  
  #'  For deployments with no images of native animals (nn = no natives):
  #'  Add 0 time for each species at these sites
  full_tt <- function(df_all_yrs, df_tt_yrs, df_tt, df_tbd) {
    #'  Vector of all deployments:
    dep <- df_all_yrs %>%
      dplyr::select(location) %>%
      distinct() %>%
      pull()
    
    #'  All possible species detected across years
    sp <- as.character(sort(unique(df_tt_yrs$common_name)))
    
    df_tt_nn <- df_tbd %>%
      #'  Retrieve only those that had no images of native species
      anti_join(df_tt, by = "location") %>%
      expand(location, season = "summer", common_name = sp) %>%
      #'  Re-join time-by-day information
      left_join(df_tbd, by = "location") %>%
      #'  Add total_duration column, which is zero in these cases
      mutate(total_duration = 0) 
    
    #'  Merge cameras w/o species detections with larger tt dataset
    df_tt_full <- df_tt %>%
      bind_rows(df_tt_nn) %>%
      arrange(location, common_name, season) %>%
      mutate(total_season_days = total_days) %>%
      dplyr::select(location:season, total_season_days, total_duration) %>% 
      #'  Remove "missing" grizzly bears at GMU6 & 10A cameras b/c not present in lower study areas
      filter(grepl("GMU1_", location) | common_name != "bear_grizzly")
    
    return(df_tt_full)
  }
  #'  Bind all possible sites & species across years
  df_all_yrs <- rbind(df_all_20s, df_all_21s, df_all_22s)
  df_tt_yrs <- rbind(df_tt_20s, df_tt_21s, df_tt_22s)
  
  df_tt_full_20s <- full_tt(df_all_yrs, df_tt_yrs = df_tt_yrs, df_tt = df_tt_20s, df_tbd = df_tbd_20s) %>%
    #'  Replace NAs for random species with no detections (needed for a couple wolverine & rodent detections)
    mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration)) %>%
    #'  Retain data from cameras where camera was operable for 30 days or more 
    filter(total_season_days >= 30)
  df_tt_full_21s <- full_tt(df_all_yrs, df_tt_yrs = df_tt_yrs, df_tt = df_tt_21s, df_tbd = df_tbd_21s)  %>%
    #'  Retain data from cameras where camera was operable for 30 days or more
    filter(total_season_days >= 30)
  df_tt_full_22s <- full_tt(df_all_yrs, df_tt_yrs = df_tt_yrs, df_tt = df_tt_22s, df_tbd = df_tbd_22s) %>%
    #'  Replace NAs for random species with no detections (needed for a bat detection)
    mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration)) %>%
    #'  Retain data from cameras where camera was operable for 30 days or more
    filter(total_season_days >= 30)
  
  tt_list <- list(df_tt_full_20s, df_tt_full_21s, df_tt_full_22s)
  names(tt_list) <- c("df_tt_full_20s", "df_tt_full_21s", "df_tt_full_22s")
  
  #' #'  SAVE!
  #' write_csv(df_tt_full_20s, "./Data/Relative abundance data/RAI Phase 2/eoe_all_20s_fov-time.csv")
  #' write_csv(df_tt_full_21s, "./Data/Relative abundance data/RAI Phase 2/eoe_all_21s_fov-time.csv")
  #' write_csv(df_tt_full_22s, "./Data/Relative abundance data/RAI Phase 2/eoe_all_22s_fov-time.csv")
  #' 
  #' save(tt_list, file = "./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time.RData")
  #' #'  This is input for Effective_detection_distance.R script
  
  
  #'  ------------------------------------
  ####  TIFC for bear investigation data  ####
  #'  ------------------------------------
  #'  Source similar script that runs bear investigation behavior data (bad_bear) 
  #'  through similar functions to calculate time spent investigating cameras and 
  #'  flags so this amount of time can be removed from bear TIFC data
  source("./Scripts/Relative_Abundance/Calculate_bad_bear_TIFC.R")
  
  #'  Reduce data set to just sites where bears interacted with camera
  bad_bear_times <- function(df_tt) {
    skinny_bad_bear <- filter(df_tt, bad_bear_behavior == TRUE) %>%
      dplyr::select(-c(season, bad_bear_behavior))
    return(skinny_bad_bear)
  }
  bear_behavior_tt <- lapply(bear_tt_list, bad_bear_times)
  
  #'  Subtract amount of time bears spent investigating cameras from larger TIFC data
  tst <- df_tt_full_22s %>%
    left_join(bear_behavior_tt[[3]], by = c("location", "common_name", "total_season_days")) %>%
    mutate(total_duration = ifelse(is.na(total_duration.y), total_duration.x, 
                                   total_duration.x - total_duration.y))
  
  
  
  #'  -------------------------
  ####  Area of Field-of-View  ####
  #'  -------------------------
  #'  Organize area of field-of-view for each site. Area (m2) measured by IDFG 
  #'  biologists using walk test, flagging, or other methods to estimate detection 
  #'  zone of each camera.   
  #'  -------------------------
  fov <- function(dets1, dets2, dets3) {
    #'  Filter detection data to just location and detection zone areas
    det_zone_20s <- dets1 %>% dplyr::select(c("NewLocationID", "Gmu", "Area_M2")) %>% unique() %>% mutate(Area_M2 = round(Area_M2, 2))
    det_zone_21s <- dets2 %>% dplyr::select(c("NewLocationID", "Gmu", "Area_M2")) %>% unique() %>% mutate(Area_M2 = round(Area_M2, 2))
    det_zone_22s <- dets3 %>% dplyr::select(c("NewLocationID", "Gmu", "Area_M2")) %>% unique() %>% mutate(Area_M2 = round(Area_M2, 2))
    
    #'  Bind across years and fill in gaps based on another year's measurement at same site when possible
    det_zone <- full_join(det_zone_20s, det_zone_21s, by = c("NewLocationID", "Gmu")) %>%
      full_join(det_zone_22s, by = c("NewLocationID", "Gmu")) %>%
      rename("Area_M2_20s" = "Area_M2.x") %>% 
      rename("Area_M2_21s" = "Area_M2.y") %>%
      rename("Area_M2_22s" = "Area_M2") %>%
      arrange(NewLocationID) %>% 
      #'  Fill in missing area measurements with another year's measurements from same site (excluding GMU 1 b/c not sampled in 2020)
      mutate(Area_M2_20s = ifelse(is.na(Area_M2_20s) & Gmu != "1" , Area_M2_21s, Area_M2_20s),
             Area_M2_21s = ifelse(is.na(Area_M2_21s), Area_M2_20s, Area_M2_21s),
             Area_M2_22s = ifelse(is.na(Area_M2_22s), Area_M2_21s, Area_M2_22s),
             #'  Fill in missing area measurements with nearby site if no measurement taken at site either year
             Area_M2_20s = ifelse(is.na(Area_M2_20s) & Gmu != "1", lead(Area_M2_20s), Area_M2_20s),
             Area_M2_21s = ifelse(is.na(Area_M2_21s), lead(Area_M2_21s), Area_M2_21s),
             Area_M2_22s = ifelse(is.na(Area_M2_22s), lead(Area_M2_22s), Area_M2_22s)) %>% #,
             #' #'  Repeat this process for the few straggler NAs
             #' Area_M2_20s = ifelse(is.na(Area_M2_20s) & Gmu != "1", lag(Area_M2_20s), Area_M2_20s),
             #' Area_M2_21s = ifelse(is.na(Area_M2_21s), lag(Area_M2_21s), Area_M2_21s),
             #' Area_M2_22s = ifelse(is.na(Area_M2_22s), lag(Area_M2_22s), Area_M2_22s),
             #' Area_M2_20s = ifelse(is.na(Area_M2_20s) & Gmu != "1", lag(Area_M2_20s), Area_M2_20s),
             #' Area_M2_21s = ifelse(is.na(Area_M2_21s), lag(Area_M2_21s), Area_M2_21s),
             #' Area_M2_22s = ifelse(is.na(Area_M2_22s), lag(Area_M2_22s), Area_M2_22s)) %>%
      #'  Remove duplicate camera where area is clearly wrong
      filter(NewLocationID != "GMU1_U_8" | Area_M2_22s != 347.20) %>%
      #'  Correct typos in area measurements
      mutate(Area_M2_22s = ifelse(NewLocationID == "GMU10A_P_14", 66.14, Area_M2_22s)) %>%
      filter(!is.na(NewLocationID))
    
    return(det_zone)
  }
  fov_area <- fov(eoe20s_allM, eoe21s_allM, eoe22s_allM)
  #'  Split field-of-view data by year and list together
  fov_area_20s <- fov_area %>% dplyr::select(-c(Area_M2_21s, Area_M2_22s)) %>% filter(!is.na(Area_M2_20s)) %>%
    rename("Area_M2" = "Area_M2_20s")
  fov_area_21s <- fov_area %>% dplyr::select(-c(Area_M2_20s, Area_M2_22s)) %>% rename("Area_M2" = "Area_M2_21s")
  fov_area_22s <- fov_area %>% dplyr::select(-c(Area_M2_20s, Area_M2_21s)) %>% rename("Area_M2" = "Area_M2_22s")
  
  area_m2 <- list(fov_area_20s, fov_area_21s, fov_area_22s)
  names(area_m2) <- c("area_m2_20s", "area_m2_21s", "area_m2_22s")
  
  save(area_m2, file = "./Data/Relative abundance data/RAI Phase 2/area_m2.RData")
  
  #'  Use these data as input for calculating density in Calculate_density.R
  