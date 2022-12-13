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
  predators <- eoe21s_allM %>%
    dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                  "OpState", "Species", "Count") %>%
    filter(Species == "bear_black" | Species == "bear_grizzly" | Species == "bobcat" | 
             Species == "coyote" | Species == "mountain_lion" | Species == "wolf") %>%
    mutate(
      Date = as.Date(Date, format = "%d-%b-%Y"),
      Time = chron(times = Time)
    ) %>%
    #'  Filter to images between July 1 and Sept 15
    filter(Date >= "2021-07-01" & Date <= "2021-09-15")
  
  
  #'  ----------------------------------------
  ####  Extract independent detection events  ####
  #'  ----------------------------------------
  #'  Create a column grouping images into "independent" event based on defined
  #'  amount of time that must elapse between images of the same species at the
  #'  same camera site.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
  dets <- predators %>%
    arrange(NewLocationID, posix_date_time) %>%
    group_by(NewLocationID, Species) %>%
    mutate(det_events = cumsum(c(1, diff(posix_date_time) > 300))) %>% # units in seconds! 
    ungroup() %>%
    #'  Retain only the first image from each unique detection event
    group_by(NewLocationID, Species, det_events) %>%
    slice(1L) %>%
    ungroup()
  
  #'  Split up by species
  bear <- dets[dets$Species == "bear_black",]
  bob <- dets[dets$Species == "bobcat",]
  coy <- dets[dets$Species == "coyote",]
  lion <- dets[dets$Species == "mountain_lion",]
  wolf <- dets[dets$Species == "wolf",]
  
 
  
  
 
    
  
  
  
  
  
  
  
  
  
  
  
  