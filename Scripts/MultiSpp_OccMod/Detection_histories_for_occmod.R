  #'  -------------------------------
  #'  Camera trap detection histories
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  -------------------------------
  #'  Script to combine species detection data with camera operational data using
  #'  camtrapR package to generate species-specific encounter histories & sampling
  #'  effort data to be used in occupancy models.
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  ----------------
  ####  Read in data  ####
  #'  ----------------
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID.RData") 
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe20w_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
    
  #'  ------------------------
  ###  Filter detection data  ####
  #'  ------------------------
  #'  1) Remove sequential problem images
  thin_detections <- function(dets, seqprobs) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    return(skinny_dets)
  }
  eoe20s_allM_skinny <- thin_detections(eoe20s_allM, eoe_seqprob_20s)
  eoe20w_allM_skinny <- thin_detections(eoe20w_allM, eoe_seqprob_20w)
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s)
  eoe22s_allM_skinny <- thin_detections(eoe22s_allM, eoe_seqprob_22s)
  
  #'  2) Filter detection data to time period of interest
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Count") %>%
      mutate(
        Date = as.Date(Date, format = "%Y-%m-%d"), #format = "%d-%b-%Y"
        Time = chron(times = Time)
      ) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      arrange(NewLocationID)
    return(dets)
  }
  eoe20s_dets <- detections(eoe20s_allM_skinny, start_date = "2020-07-01", end_date = "2020-09-15")
  eoe20w_dets <- detections(eoe20w_allM_skinny, start_date = "2020-12-01", end_date = "2021-02-01")
  eoe21s_dets <- detections(eoe21s_allM_skinny, start_date = "2021-07-01", end_date = "2021-09-15")
  eoe22s_dets <- detections(eoe22s_allM_skinny, start_date = "2022-07-01", end_date = "2022-09-15")

  
  #'  -----------------------------------------
  ####  Generate independent detection events  ####
  #'  -----------------------------------------
  #'  Create a column grouping images into "independent" detection events based 
  #'  on defined amount of time that must elapse between images of the same species 
  #'  at the same camera site. Filter to the first image of each event and then
  #'  split into species-specific detection events.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
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
  eoe20s_det_events <- unique_detections(eoe20s_dets, elapsed_time = 300)
  eoe20w_det_events <- unique_detections(eoe20w_dets, elapsed_time = 300)
  eoe21s_det_events <- unique_detections(eoe21s_dets, elapsed_time = 300)
  eoe22s_det_events <- unique_detections(eoe22s_dets, elapsed_time = 300)
    
  #'  Save for making summary tables
  save(eoe20s_det_events, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/eoe20s_det_events.RData")
  save(eoe20w_det_events, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/eoe20w_det_events.RData")
  save(eoe21s_det_events, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/eoe21s_det_events.RData")
  save(eoe22s_det_events, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/eoe22s_det_events.RData")
    
  
  #'  --------------------------
  ####  Camera Operation Table  ####
  #'  -------------------------- 
  #'  Create a matrix with each camera & dates it deployed
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camera_operation_tbl <- function(cams) {
    cams <- arrange(cams, NewLocationID)
    #'  Add one day to retrieval date when setup & retrieval dates are the same
    #'  Necessary for cameraOperation function below... so annoying
    same_startend <- cams[cams$Setup_date == cams$Retrieval_date,]
    same_startend <- mutate(same_startend, Retrieval_date = as.Date(Retrieval_date) + 1)
    same_startend <- mutate(same_startend, Problem1_to = as.Date(Problem1_to) +1)
    cams[match(same_startend$NewLocationID, cams$NewLocationID),] <- same_startend
    
    #'  Make sure camera data are organized by NewLocationID
    cams <- arrange(cams, NewLocationID)
    
    #'  Create camera operation table
    camop_problem <- cameraOperation(CTtable = cams,
                                     stationCol = "NewLocationID",
                                     setupCol = "Setup_date",
                                     retrievalCol = "Retrieval_date",
                                     hasProblems = TRUE,
                                     dateFormat = "%Y-%m-%d", 
                                     writecsv = FALSE) 
    
    return(camop_problem)
  }
  eoe20s_probs <- camera_operation_tbl(eoe_probcams_20s) 
  eoe20w_probs <- camera_operation_tbl(eoe_probcams_20w) 
  eoe21s_probs <- camera_operation_tbl(eoe_probcams_21s)
  eoe22s_probs <- camera_operation_tbl(eoe_probcams_22s)
  
  #'  Make sure two data streams match up
  #'  OK to have cameras in camera setup data that aren't in detection data but 
  #'  can't have cameras in detection data that aren't in camera operation table
  match_datsets <- function(det_events, stations) {
    det_cams <- as.data.frame(unique(det_events$NewLocationID))
    det_cams$Source1 <- "Detection Cameras"
    colnames(det_cams) <- c("NewLocationID", "Source1")
    cam_stations <- as.data.frame(rownames(stations))
    cam_stations$Source2 <- "Operational Cameras"
    colnames(cam_stations) <- c("NewLocationID", "Source2")
    match_cams <- full_join(det_cams, cam_stations, by = "NewLocationID")
    return(match_cams)
  } 
  #'  NAs in Source1 column are OK, NAs in Source2 column are a problem
  eoe20s_match_cams <- match_datsets(det_events = eoe20s_det_events, stations = eoe20s_probs)
  eoe20w_match_cams <- match_datsets(det_events = eoe20w_det_events, stations = eoe20w_probs)
  eoe21s_match_cams <- match_datsets(det_events = eoe21s_det_events, stations = eoe21s_probs)
  eoe22s_match_cams <- match_datsets(det_events = eoe22s_det_events, stations = eoe22s_probs)

  
  #'  -----------------------
  ####  Detection Histories  ####
  #'  -----------------------
  #'  Function to create species-specific detection histories for each species 
  #'  and create a sampling effort matrix
  #'  
  #'  Key arguments:
  #'  -occasionLength: currently using 7 day sampling occasions 
  #'  -day1: sampling occasion begins upon station setup date ("station"), 
  #'   first day of survey ("survey"), or a specific date (e.g., "2020-07-01")
  #'   FYI this defines start date but NOT end date so DH goes until camera pulls
  #'  -includeEffort: compute # active trap days/station/occasion- effects DH
  #'   if FALSE then when camera is not set up or malfunctioning (NA or 0 in
  #'   camOp) during part of sampling occasion DH will be 0 for that occasion
  #'  -scaleEffort: center & scale effort (can also do this later in model)
  #'  -dateAsOccasionNames: if day1 = "survey" then uses 1st and last day of 
  #'   occasion as name for each sampling occasion
  #'  -output: return binary detections or counts of detections; don't want to
  #'   use "count" right now b/c would count each image not independent events
  #'  -occasionStartTime: time of day (full hour) at which to begin occasions
  #'   default is midnight (0)
  #'  -unmarkedMultFrameInput: create input for multi-season occmod in unmarked
  #'   ONLY use if running multi-season models & need to add more info to camop
  #'   
  #'   FYI: July 1 - Sept 15 = 11 1-wk sampling periods
  #'        Dec 1 - Feb 1 = 9 1-wk sampling periods

  DH <- function(dets, cam_probs, spp, start_date, y, rm_rows, oc) {
    det_hist <- detectionHistory(recordTable = dets,
                                 camOp = cam_probs,
                                 stationCol = "NewLocationID",
                                 speciesCol = "Species",
                                 recordDateTimeCol = "posix_date_time",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 7,
                                 day1 = start_date, 
                                 datesAsOccasionNames = TRUE,
                                 # occasionStartTime = 12, # starts at noon
                                 timeZone = "America/Edmonton",
                                 output = y,
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 # writecsv = TRUE,
                                 outDir = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories")
    
    #'  Reduce detection histories to sampling occasions of interest (drop extra
    #'  occasions after focal period of interest)
    #'  This will vary depending on length of primary period and sampling occasions!
    short_dh <- det_hist[[1]][,1:oc]
    short_effort <- det_hist[[2]][,1:oc]
    
    short_dh <- short_dh[-rm_rows,]
    short_effort <- short_effort[-rm_rows,]
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create season-specific detection histories for species listed below
  spp_smr <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf")
  spp_wtr <- list("bobcat", "coyote", "mountain_lion", "wolf")
  
  rm_rows_eoe20s <- c(61, 79, 82, 98, 125, 157, 171, 177, 178, 181, 186, 192, 200, 214, 228, 235, 236, 259, 311, 334, 346, 361, 371, 379, 380, 385, 419, 433, 437, 439, 458, 493)
  DH_eoe20s_predators <- lapply(spp_smr, DH, dets = eoe20s_det_events, cam_probs = eoe20s_probs, start_date = "2020-07-01", y = "binary", rm_rows = rm_rows_eoe20s, oc = 11) 
    
  rm_rows_eoe20w <- c(7, 16, 32, 43, 123, 138, 195, 215, 227, 242, 252, 268) 
  DH_eoe20w_predators <- lapply(spp_wtr, DH, dets = eoe20w_det_events, cam_probs = eoe20w_probs, start_date = "2020-12-01", y = "binary", rm_rows = rm_rows_eoe20w, oc = 9) 
    
  rm_rows_eoe21s <- c(6, 106, 112, 116, 127, 145, 147, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682)
  DH_eoe21s_predators <- lapply(spp_smr, DH, dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01", y = "binary", rm_rows = rm_rows_eoe21s, oc = 11)
    
  rm_rows_eoe22s <- c(13, 29, 85, 108, 113, 154, 173, 175, 199, 208, 215, 216, 238, 242, 256, 257, 269, 286, 292, 349, 352, 357, 398, 427, 429, 450, 452, 455, 489, 501, 512, 516, 584, 632, 637, 656, 659, 661, 687, 689)
  DH_eoe22s_predators <- lapply(spp_smr, DH, dets = eoe22s_det_events, cam_probs = eoe22s_probs, start_date = "2022-07-01", y = "binary", rm_rows = rm_rows_eoe22s, oc = 11)
  # rownmbr <- seq(1:nrow(short_dh))
  # tst <- cbind(short_dh, rownmbr)
   
  #'  -----------------------
  ####  Max Count of Wolves  ####
  #'  -----------------------
  #'  Count number of wolf detections per sampling occasion
  count_eoe20s_wolf <- DH(spp = "wolf", dets = eoe20s_det_events, cam_probs = eoe20s_probs, start_date = "2020-07-01", y = "count", rm_rows = rm_rows_eoe20s,oc = 11)  
  count_eoe20w_wolf <- DH(spp = "wolf", dets = eoe20w_det_events, cam_probs = eoe20w_probs, start_date = "2020-12-01", y = "count", rm_rows = rm_rows_eoe20w,oc = 9)  
  count_eoe21s_wolf <- DH(spp = "wolf", dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01", y = "count", rm_rows = rm_rows_eoe21s, oc = 11)
  count_eoe22s_wolf <- DH(spp = "wolf", dets = eoe22s_det_events, cam_probs = eoe22s_probs, start_date = "2022-07-01", y = "count", rm_rows = rm_rows_eoe22s, oc = 11)
  
  #'  Count maximum number of wolves detected in a single image per detection event
  #'  at each camera, then find the average to represent the average minimum group
  #'  size detected at each camera
  avg_min_group_size <- function(dets, stations, elapsed_time) {
    dets <- arrange(dets, NewLocationID)
    min_group_size <- dets %>%
      filter(Species == "wolf") %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Identify unique detection events based on elapsed amount of time
      group_by(NewLocationID) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds! 
      ungroup() %>%
      #'  Find the maximum number of individuals detected in a single image per
      #'  detection event --> this is the minimum group size per detection
      group_by(NewLocationID, det_events) %>%
      summarise(max_count = max(Count)) %>%
      ungroup() %>%
      #'  Take the average of the minimum group size for each camera station
      group_by(NewLocationID) %>%
      summarise(avg_min_group_size = mean(max_count)) %>%
      ungroup() %>%
      mutate(avg_min_group_size = round(avg_min_group_size, 3)) %>%
      #'  Join with camera stations - fills in NAs where no wolves were detected
      full_join(stations, by = "NewLocationID") %>%
      dplyr::select(c("NewLocationID", "Lat", "Long", "avg_min_group_size")) %>%
      arrange(NewLocationID) %>%
      #'  Replace NAs with 0s - no wolves detected at these sizes so average minimum
      #'  group size is 0 (to the best of our knowledge)
      mutate(avg_min_group_size = ifelse(is.na(avg_min_group_size), 0, avg_min_group_size))
    
    return(min_group_size)
  }
  min_group_size_eoe20s <- avg_min_group_size(eoe20s_dets, stations = eoe_probcams_20s, elapsed_time = 300)
  min_group_size_eoe20w <- avg_min_group_size(eoe20w_dets, stations = eoe_probcams_20w, elapsed_time = 300)
  min_group_size_eoe21s <- avg_min_group_size(eoe21s_dets, stations = eoe_probcams_21s, elapsed_time = 300)
  min_group_size_eoe22s <- avg_min_group_size(eoe22s_dets, stations = eoe_probcams_22s, elapsed_time = 300)
  
  
  #'  -------------------
  ####  Sampling effort  ####
  #'  -------------------
  #'  Count number of operational days and hours per camera
  sampling_effort <- function(dh) {
    #' Number of days per sampling occasion each camera was operational
    #' Second list generated by camtrapR's detection histories function
    camdays <- dh[[1]][[2]]
    effort <- as.data.frame(camdays) %>%
      rownames_to_column(var = "NewLocationID") %>%
      #'  Count total operational days and what that is in hours
      mutate(ndays = rowSums(.[2:12], na.rm = T),
             nhrs = ndays*24) 
    return(effort)
  }
  effort_20s <- sampling_effort(DH_eoe20s_predators)
  effort_21s <- sampling_effort(DH_eoe21s_predators)
  effort_22s <- sampling_effort(DH_eoe22s_predators)
  
  camdays_20w <- DH_eoe20w_predators[[1]][[2]]
  effort_20w <- as.data.frame(camdays_20w) %>%
    rownames_to_column(var = "NewLocationID") %>%
    #'  Count total operational days and what that is in hours
    mutate(ndays = rowSums(.[2:10], na.rm = T),
           nhrs = ndays*24) 
  
  
  #' #'  --------
  #' ####  SAVE  ####
  #' #'  --------
  #' #'  Seasonal detection histories
  #' save(DH_eoe20s_predators, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  #' save(DH_eoe20w_predators, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20w_predators.RData")
  #' save(DH_eoe21s_predators, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  #' save(DH_eoe22s_predators, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe22s_predators.RData")
  #' #'  Seasonal wolf detection data
  #' save(count_eoe20s_wolf, file = "./Data/Wolf count data/count_eoe20s_wolf.RData")
  #' save(count_eoe20w_wolf, file = "./Data/Wolf count data/count_eoe20w_wolf.RData")
  #' save(count_eoe21s_wolf, file = "./Data/Wolf count data/count_eoe21s_wolf.RData")
  #' save(count_eoe22s_wolf, file = "./Data/Wolf count data/count_eoe22s_wolf.RData")
  #' #'  Minimum group size counts
  #' save(min_group_size_eoe20s, file = "./Data/Wolf count data/min_group_size_eoe20s.RData")
  #' save(min_group_size_eoe20w, file = "./Data/Wolf count data/min_group_size_eoe20w.RData")
  #' save(min_group_size_eoe21s, file = "./Data/Wolf count data/min_group_size_eoe21s.RData")
  #' save(min_group_size_eoe22s, file = "./Data/Wolf count data/min_group_size_eoe22s.RData")
  #' #' Sampling effort
  #' save(effort_20s, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  #' save(effort_20w, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20w.RData")
  #' save(effort_21s, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData")
  #' save(effort_22s, file = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe22s.RData")
  
  #'  Save for publication
  DetectionHist_Smr20 <- list(DH_eoe20s_predators[[1]][[1]], DH_eoe20s_predators[[2]][[1]], DH_eoe20s_predators[[3]][[1]], DH_eoe20s_predators[[4]][[1]], DH_eoe20s_predators[[5]][[1]])
  DetectionHist_Smr21 <- list(DH_eoe21s_predators[[1]][[1]], DH_eoe21s_predators[[2]][[1]], DH_eoe21s_predators[[3]][[1]], DH_eoe21s_predators[[4]][[1]], DH_eoe21s_predators[[5]][[1]])
  DetectionHist_Smr22 <- list(DH_eoe22s_predators[[1]][[1]], DH_eoe22s_predators[[2]][[1]], DH_eoe22s_predators[[3]][[1]], DH_eoe22s_predators[[4]][[1]], DH_eoe22s_predators[[5]][[1]])
  SamplingEffort_Smr20 <- effort_20s
  SamplingEffort_Smr21 <- effort_21s
  SamplingEffort_Smr22 <- effort_22s
  save(DetectionHist_Smr20, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr20.RData")
  save(DetectionHist_Smr21, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr21.RData")
  save(DetectionHist_Smr22, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr22.RData")
  save(SamplingEffort_Smr20, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr20.RData")
  save(SamplingEffort_Smr21, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr21.RData")
  save(SamplingEffort_Smr22, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr22.RData")
  