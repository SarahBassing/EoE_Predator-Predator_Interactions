  #'  -------------------------------
  #'  Camera trap detection histories
  #'  IDCFWRU - Predator Interactions
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
  
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  
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
    
    
    
    #'  Create camera operation table
    camop_problem <- cameraOperation(CTtable = cams,
                                     stationCol = "NewLocationID",
                                     setupCol = "Setup_date",
                                     retrievalCol = "Retrieval_date",
                                     hasProblems = TRUE,
                                     dateFormat = "%Y-%m-%d", 
                                     writecsv = FALSE) 
    
    # 2021-07-01probs <- as.data.frame(camop_problem)
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

  DH <- function(dets, cam_probs, spp, start_date) {
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
                                 output = "binary",
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 # writecsv = TRUE,
                                 outDir = "./Data/Detection_Histories")
    
    return(det_hist)
  }
  #'  Create species and season-specific detection histories
  spp <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf")
  DH_eoe21s_predators <- lapply(spp, DH, dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01")
  
  save(DH_eoe21s_predators, file = "./Data/Detection_Histories/DH_eoe21s_predators.RData")
  
 
    
  
  
  
  
  
  
  
  
  
  
  
  