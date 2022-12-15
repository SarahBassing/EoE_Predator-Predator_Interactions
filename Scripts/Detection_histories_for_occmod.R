  #'  -------------------------------
  #'  Camera trap detection histories
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  -------------------------------
  #'  Script to combine species detection data with camera operational data using
  #'  camtrapR package to generate species-specific encounter histoies and sampling
  #'  effort data to be used in occupancy models.
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
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  #'  Fix a couple of capitalization issues in the GMU part of NewLocationID
  eoe_probcams_21s <- eoe_probcams_21s %>%
    mutate(NewLocationID = toupper(NewLocationID))
  
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  #'  Fix a couple of capitalization issues in the GMU part of NewLocationID
  eoe21s_allM <- eoe21s_allM %>%
    mutate(NewLocationID = toupper(NewLocationID))
  
  #'  Filter detection data to focal species and time period of interest
  detections <- function(dets, start_date, end_date) {
    dets <- eoe21s_allM %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Count") %>%
      # filter(Species == "bear_black" | Species == "bear_grizzly" | Species == "bobcat" | 
      #          Species == "coyote" | Species == "mountain_lion" | Species == "wolf") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time)
      ) %>%
      #'  Filter to images between July 1 and Sept 15
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID))
    return(dets)
  }
  eoe21s_dets <- detections(eoe21s_allM, start_date = "2021-07-01", end_date = "2021-09-15")
  
  
  #'  ----------------------------------------
  ####  Extract independent detection events  ####
  #'  ----------------------------------------
  #'  Create a column grouping images into "independent" detection events based 
  #'  on defined amount of time that must elapse between images of the same species 
  #'  at the same camera site. Filter to the first image of each event and then
  #'  split into species-specific detection events.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    det_events <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      group_by(NewLocationID, Species) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds! 
      ungroup() %>%
      #'  Retain only the first image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup()
    
    return(det_events)
    
    #' #'  Filter by focal species
    #' bear <- det_events[det_events$Species == "bear_black",]
    #' bob <- det_events[det_events$Species == "bobcat",]
    #' coy <- det_events[det_events$Species == "coyote",]
    #' lion <- det_events[det_events$Species == "mountain_lion",]
    #' wolf <- det_events[det_events$Species == "wolf",]
    #' 
    #' #'  List unique predator detection events
    #' pred_list <- list(bear, bob, coy, lion, wolf)
    #' 
    #' return(pred_list)
  }
  eoe21s_det_events <- unique_detections(eoe21s_dets, elapsed_time = 300)
  
  
  #'  --------------------------
  ####  Camera Operation Table  ####
  #'  -------------------------- 
  #'  Create a matrix with each camera & dates it deployed
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camera_operation_tbl <- function(cams) {
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
    
    # probs <- as.data.frame(camop_problem)
    return(camop_problem)
  }
  eoe21s_probs <- camera_operation_tbl(eoe_probcams_21s)
  
  
  #'  Make sure two data streams match up
  #'  OK to have cameras in camera setup data that aren't ind detection data but 
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
  eoe21s_match_cams <- match_datsets(det_events = eoe21s_det_events, stations = eoe21s_probs)

  
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
  #'   FYI: July 1 - Sept 15, 2021 = 11 1-wk sampling periods

  DH <- function(dets, cam_probs, spp, start_date, y, oc) {
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
                                 outDir = "./Data/Detection_Histories")
    
    #'  Reduce detection histories to sampling occasions of interest (drop extra
    #'  occasions after focal period of interest)
    #'  This will vary depending on length of primary period and sampling occasions!
    short_dh <- det_hist[[1]][,1:oc]
    short_effort <- det_hist[[2]][,1:oc]
    
    short_dh <- short_dh[-c(6, 106, 112, 116, 127, 145, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682),]
    short_effort <- short_effort[-c(6, 106, 112, 116, 127, 145, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682),]
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create species and season-specific detection histories
  spp <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf")
  DH_eoe21s_predators <- lapply(spp, DH, dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01", y = "binary", oc = 11)
    
  # save(DH_eoe21s_predators, file = "./Data/Detection_Histories/DH_eoe21s_predators.RData")
  
  
  #'  Count number of wolf detections per sampling occasion
  count_eoe21s_wolf <- DH(spp = "wolf", dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01", y = "count", oc = 11)
  
  # save(count_eoe21s_wolf, file = "./Data/Wolf count data/count_eoe21s_wolf.RData")
  
 
  ####  NEED A BETTER WAY OF DUMPING ROWS WHERE CAMERA WAS COMPLETELY INOPERABLE  ####
  #'  These are problem cameras for eoe smr21
  #'  c(6, 106, 112, 116, 127, 145, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682)
    
  
  
  #'  -----------------------
  ####  Max Count of Wolves  ####
  #'  -----------------------
  #'  Count maximum number of wolves detected in a single image per detection event
  #'  at each camera, then find the average to represent the average minimum group
  #'  size detected at each camera
  avg_min_group_size <- function(dets, stations, elapsed_time) {
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
  min_group_size_eoe21s <- avg_min_group_size(eoe21s_dets, stations = eoe_probcams_21s, elapsed_time = 300)
  
  #'  Save
  # save(min_group_size_eoe21s, file = "./Data/Wolf count data/min_group_size_eoe21s.RData")
  
  
  
  
  