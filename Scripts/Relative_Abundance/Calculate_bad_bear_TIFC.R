  #'  ----------------------------------------
  #'  Calculate time-in-front-of-camera (TIFC) for bear behavior images
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  October 2023
  #'  ----------------------------------------
  #'  Script is sourced by Calculate_time-in-front-of-camera.R and calculates TIFC 
  #'  based on bear behavior - natural behaviors vs investigation behaviors.
  #'  Methods: Huggard et al. 2018, Warbington & Boyce 2020, and Becker et al. 2022.
  #'  Code modified from Becker et al. 2022 repository: https://github.com/mabecker89/tifc-method/tree/main
  #'  Leave probabilities taken directly from Becker et al. (2022)
  #'  ----------------------------------------
  
  #'  -------------------------------
  ####  Load and format input data  ####
  #'  -------------------------------
  #'  Load images labeled with bear investigation behavior 
  bears <- read_csv("./Data/IDFG camera data/bears.csv") 
  bad_bear <- bears %>%
    mutate(Date = as.Date(DateTime, format = "%Y-%m-%d"),
           Date = as.character(Date),
           Year = year(Date),
           Species = "bear_black") %>%
    relocate(Date, .before = DateTime) %>%
    relocate(Species, .after = DateTime) %>%
    rename("posix_date_time" = "DateTime") %>%
    dplyr::select(-c(RootFolder, RelativePath))
  bad_bear_smr20 <- filter(bad_bear, Year == "2020") %>% dplyr::select(-Year)
  bad_bear_smr21 <- filter(bad_bear, Year == "2021") %>% dplyr::select(-Year)
  bad_bear_smr22 <- filter(bad_bear, Year == "2022") %>% dplyr::select(-Year)
  
  #'  Connect bad bear behavior images to larger data set so looks like other data
  bad_bear_20s <- left_join(bad_bear_smr20, eoe20s_allM, by = c("File", "Date", "Species", "posix_date_time")) %>%
    filter(!is.na(NewLocationID))
  bad_bear_21s <- left_join(bad_bear_smr21, eoe21s_allM, by = c("File", "Date", "Species", "posix_date_time")) %>%
    filter(!is.na(NewLocationID))
  bad_bear_22s <- left_join(bad_bear_smr22, eoe22s_allM, by = c("File", "Date", "Species", "posix_date_time")) %>%
    filter(!is.na(NewLocationID)) 
  
  #'  Grab mislabeled bear images (mostly smooshy face pictures of other species)
  not_a_bear <- function(bear_pix) {
    no_bear <- filter(bear_pix, NoBear == "TRUE")
    return(no_bear)
  }
  notabear_20s <- not_a_bear(bad_bear_smr20)
  notabear_21s <- not_a_bear(bad_bear_smr21)
  notabear_22s <- not_a_bear(bad_bear_smr22)
  
  #'  Remove problem images from bad bear data and simplify df
  thin_bad_bears <- function(dets, seqprobs, notabear) {
    noprob_dets <- dets[!(dets$File %in% seqprobs$File),]
    skinny_bears <- noprob_dets[!(noprob_dets$File %in% notabear$File),]
    
    clean_bears <- skinny_bears %>%
      #'  Make sure Count column has a value
      mutate(Count = ifelse(Count == 0, 1, Count),
             Count = ifelse(is.na(Count), 1, Count),
             #'  Create column for all types of investigation behaviors
             bad_bear = ifelse(Cam_behavior == "TRUE", "TRUE", Flag_behavior),
             #'  Add gap class column
             gap_class = NA) %>% 
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      arrange(NewLocationID, posix_date_time) %>%
      dplyr::select(c("NewLocationID", "posix_date_time", "Species", "Count", "gap_class", "Cam_behavior", "Flag_behavior", "bad_bear")) %>%
      mutate(gap_class = NA)
    names(clean_bears) <- c("location", "date_detected", "common_name", "number_individuals", "gap_class", "cam_behavior", "flag_behavior", "bad_bear_behavior")
    
    return(clean_bears)
  }
  df_bad_bears_20s <- thin_bad_bears(bad_bear_20s, seqprobs = eoe_seqprob_20s, notabear = notabear_20s) 
  df_bad_bears_21s <- thin_bad_bears(bad_bear_21s, seqprobs = eoe_seqprob_21s, notabear = notabear_21s)
  df_bad_bears_22s <- thin_bad_bears(bad_bear_22s, seqprobs = eoe_seqprob_22s, notabear = notabear_22s)
  
  bad_bear_list <- list(df_bad_bears_20s, df_bad_bears_21s, df_bad_bears_22s)
  
  #'  Sampling effort data (number of days cameras were operational)
  df_tbd <- list(df_tbd_20s, df_tbd_21s, df_tbd_22s)
  
  #'  Seasonal start/end dates (julian day)
  smr_starts <- list(summer20.start.j, summer21.start.j, summer22.start.j)
  smr_ends <- list(summer20.end.j, summer21.end.j, summer22.end.j)
  
  #'  ------------------
  ####  Calculate TIFC  ####
  #'  ------------------
  ##### Step 1. Identify Detection Series  #####
  #'  ------------------------------------------
  bear_img_series <- function(df_all) {
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
             #'  Lagged deployment
             location_previous = lag(location),
             #'  Was it a different deployment?
             diff_location = ifelse(location != location_previous, TRUE, FALSE),
             #'  Flag investigation behavior images 
             behavior_check = ifelse(diff_location == FALSE & bad_bear_behavior == TRUE, 1, 0), 
             #'  Did behavior change?
             diff_behavior = ifelse(behavior_check != lag(behavior_check), TRUE, FALSE),
             #'  Flag gaps that will need checking
             gap_check = ifelse(diff_location == FALSE & (diff_time_previous <= 120 & diff_time_previous >= 20), 1, 0),
             #'  Lagged gap class (Not actually necessary since currently no gap class data) 
             gap_class_previous = replace_na(lag(gap_class)), 
             gap_class_previous = ifelse(is.na(gap_class_previous), "", gap_class_previous),
             #'  Identify new series based on being different deployment, behavior classification, or gaps longer than 120s
             diff_series = ifelse(diff_location == TRUE | diff_behavior == TRUE | diff_time_previous > 120, 1, 0), 
             #'  Number series
             series_num = c(0, cumsum(diff_series[-1])),
             #' #'  Flag gaps that require probabilistic time assignment
             #' behavior_prob = replace_na(ifelse(behavior_check == 1, 1, 0), 0),
             #'  Flag gaps that require probabilistic time assignment
             gap_prob = replace_na(ifelse(gap_check == 1 & (gap_class_previous == "" | gap_class_previous == "U"), 1, 0), 0)) %>% 
      group_by(series_num) %>%
      mutate(diff_time_previous = ifelse(row_number() == 1, 0, diff_time_previous),
             diff_time_next = ifelse(row_number() == n(), 0, diff_time_next),
             gap_group =  "Bears") %>%
      ungroup() %>%
      #'  Join gap leaving predictions
      left_join(df_leave_prob_pred, by = c("gap_group", "diff_time_previous" = "diff_time")) %>%
      #'  Even though rarely any real gaps in investigation behavior, still adjusting
      #'  for gap probability so amount of time subtracted from larger bear TIFC
      #'  is calculated the same way as the original TIFC data
      mutate(pred = replace_na(pred, 1),
             diff_time_previous_adj = ifelse(gap_prob == 1, diff_time_previous * (1 - pred), diff_time_previous),
             diff_time_next_adj = ifelse(lead(gap_prob == 1), diff_time_next * (1 - lead(pred)), diff_time_next),
             #'  On the rare occasion a single image of a species is last image at camera, 
             #'  force time_diff_next_adj = 0 so not left as NA (required for later calculations)
             #'  e.g., GMU6_U_97 Summer 2020 has 1 img of coyote & last photo taken by camera
             diff_time_next_adj = ifelse(is.na(diff_time_next_adj), 0, diff_time_next_adj))
    return(df_series)
  }
  bad_bear_series <- lapply(bad_bear_list, bear_img_series)
  
  
  ##### Step 2. Calculate time between photos (tbp), by species #####
  #'  ---------------------------------------------------------------
  bear_tbp <- function(df_series) {
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
  bear_tbp <- lapply(bad_bear_series, bear_tbp)
  
  
  #####  Step 3. Calculate total time in front of the camera, by series (tts = Total Time by Series)  #####
  #'  -----------------------------------------------------------------------------------------------------
  bear_tts <- function(df_series, df_tbp) {
    df_tts <- df_series %>%
      left_join(df_tbp, by = "common_name") %>%
      group_by(series_num) %>%
      #'  Check whether the image was first or last in a series
      mutate(bookend = ifelse(row_number() == 1 | row_number() == n(), 1, 0),
             #'  Calculate time for each individual image
             image_time = ifelse(bookend == 1,
                                 ((diff_time_previous_adj + diff_time_next_adj) / 2),
                                 (diff_time_previous_adj + diff_time_next_adj) / 2),
             #'   Multiply image time by the number of animals present
             image_time_ni = image_time * number_individuals) %>%
      #'  Group by common name as well to add it as a variable to output
      group_by(behavior_check, .add = TRUE) %>%
      #'  Calculate total time and number of images for each series
      summarise(n_images = n(),
                series_total_time = sum(image_time_ni)) %>%
      ungroup() %>%
      mutate(common_name = "bear_black")
    return(df_tts)
  }
  bear_tts <- mapply(bear_tts, df_series = bad_bear_series, df_tbp = bear_tbp, SIMPLIFY = FALSE)
  
  
  #####  Step 4. Calculate total time in front of camera, by deployment and species (tt = total time)  #####
  #'  ------------------------------------------------------------------------------------------------------
  bear_tt <- function(df_series, df_tts, df_tbd, summer.start.j, summer.end.j) {
    df_tt <- df_series %>%
      group_by(series_num) %>%
      arrange(date_detected, .by_group = TRUE) %>%
      filter(row_number() == 1) %>%
      left_join(df_tts, by = c("series_num", "common_name")) %>%
      dplyr::select(location, date_detected, common_name, series_num, series_total_time, bad_bear_behavior) %>%
      ungroup() %>%
      mutate(julian = as.numeric(format(date_detected, "%j")),
             season = ifelse(julian >= summer.start.j & julian <= summer.end.j, "summer", "winter")) %>%
      mutate_at(c("location", "common_name", "season"), factor) %>%
      group_by(location, common_name, season, bad_bear_behavior, .drop = FALSE) %>%
      summarise(total_duration = sum(series_total_time)) %>%
      ungroup() %>%
      mutate_if(is.factor, as.character) %>%
      left_join(df_tbd, by = "location")
    return(df_tt)
  }
  bear_tt <- mapply(bear_tt, df_series = bad_bear_series, df_tts = bear_tts, df_tbd = df_tbd, summer.start.j = smr_starts, summer.end.j = smr_ends, SIMPLIFY = FALSE)
  
  #'  For deployments with no images of native animals (nn = no natives):
  #'  Add 0 time for each species at these sites
  bear_full_tt <- function(df_all_yrs, df_tt_yrs, df_tt, df_tbd) {
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
      mutate(total_duration = 0,
             bad_bear_behavior = "FALSE") %>%
      relocate(season, .after = common_name) %>%
      relocate(total_days, .after = total_duration) %>%
      relocate(bad_bear_behavior, .after = season)
    
    #'  Merge cameras w/o species detections with larger tt dataset
    df_tt_full <- df_tt %>%
      bind_rows(df_tt_nn) %>%
      arrange(location, common_name, season) %>%
      mutate(total_season_days = total_days) %>%
      dplyr::select(-total_days) %>% 
      #'  Remove "missing" grizzly bears at GMU6 & 10A cameras b/c not present in lower study areas
      filter(grepl("GMU1_", location) | common_name != "bear_grizzly") %>%
      #'  Filter to season and species of interest
      filter(season == "summer") %>%
      filter(common_name == "bear_black") %>%
      mutate(bad_bear_behavior = ifelse(is.na(bad_bear_behavior), FALSE, bad_bear_behavior))
    
    return(df_tt_full)
  }
  #'  Bind all possible sites & species across years
  df_bad_bear_all_yrs <- rbind(df_bad_bears_20s, df_bad_bears_21s, df_bad_bears_22s)
  df_bear_tt_yrs <- rbind(bear_tt[[1]], bear_tt[[2]], bear_tt[[3]])
  
  df_bear_tt_20s <- bear_full_tt(df_bad_bear_all_yrs, df_tt_yrs = df_bear_tt_yrs, df_tt = bear_tt[[1]], df_tbd = df_tbd[[1]]) %>%
    #'  Replace NAs for random species with no detections (needed for a couple wolverine & rodent detections)
    mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration)) %>%
    #'  Retain summer data from cameras where camera was operable for 30 days or more 
    filter(total_season_days >= 30) 
  df_bear_tt_21s <- bear_full_tt(df_bad_bear_all_yrs, df_tt_yrs = df_bear_tt_yrs, df_tt = bear_tt[[2]], df_tbd = df_tbd[[2]])  %>%
    #'  Retain summer data from cameras where camera was operable for 30 days or more
    filter(total_season_days >= 30) 
  df_bear_tt_22s <- bear_full_tt(df_bad_bear_all_yrs, df_tt_yrs = df_bear_tt_yrs, df_tt = bear_tt[[3]], df_tbd = df_tbd[[3]]) %>%
    #'  Replace NAs for random species with no detections (needed for a bat detection)
    mutate(total_duration = ifelse(is.na(total_duration), 0, total_duration)) %>%
    #'  Retain summer data from cameras where camera was operable for 30 days or more
    filter(total_season_days >= 30) 
  
  bear_tt_list <- list(df_bear_tt_20s, df_bear_tt_21s, df_bear_tt_22s)
  names(bear_tt_list) <- c("df_bear_tt_20s", "df_bear_tt_21s", "df_bear_tt_22s")
  
  #' #'  SAVE!
  #' write_csv(df_bear_tt_20s, "./Data/Relative abundance data/RAI Phase 2/eoe_bear_20s_fov-time.csv")
  #' write_csv(df_bear_tt_21s, "./Data/Relative abundance data/RAI Phase 2/eoe_bear_21s_fov-time.csv")
  #' write_csv(df_bear_tt_21s, "./Data/Relative abundance data/RAI Phase 2/eoe_bear_22s_fov-time.csv")
  #' 
  #' save(bear_tt_list, file = "./Data/Relative abundance data/RAI Phase 2/eoe_bear_fov-time.RData")
  
  
  
