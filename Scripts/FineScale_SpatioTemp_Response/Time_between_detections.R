  #'  ---------------------------------------
  #'  Calculate time between detections
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------------
  #'  Calculate time-between-detections of one predator species following the 
  #'  detection of another predator species at the same camera site.
  #'  ---------------------------------------
  

  #'  Load libraries
  library(data.table)
  library(lubridate)
  library(chron)
  # library(sp)
  # library(raster)
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
      filter(Species == "elk" | Species == "moose" | Species == "muledeer" | Species == "whitetaileddeer" |
               Species == "rabbit_hare" | Species == "human" | Species == "dog_domestic" | 
               Species == "horse" |  Species == "cattle_cow" | Species == "cat_domestic" |
               Vehicle == "TRUE") %>%
      dplyr::select("NewLocationID", "CamID", "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time),
        Species = ifelse(Vehicle == "TRUE", "human_motorized", Species),
        mutate(Category = ifelse(Species == "Bobcat" | Species == "Black bear" | 
                                   Species == "Coyote" | Species == "Mountain lion" | 
                                   Species == "Wolf", "Predator", "Other")) %>%
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
  
  #'  ------------------------------------------
  ####  Summarize independent detection events  ####
  #'  ------------------------------------------
  #'  Generate "independent" detection events based on defined amount of time 
  #'  elapsing between sequential images of same species at a camera, then count
  #'  the number of unique detection events of a given species at each camera.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    #'  Generate unique detection events
    det_events <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Flag images of same species at same camera as being a different detection event
      #'  when time since last image of that species is greater than defined time interval
      group_by(Species, NewLocationID) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds!
      ungroup() #%>%
      #' #'  Retain only the first image from each unique detection event
      #' group_by(NewLocationID, Species, det_events) %>%
      #' slice(1L) %>%
      #' ungroup()

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
    firstimg <- capdata[capdata$Category == "Predator",] %>%
      group_by(det_events) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Det_type = "first")
    #'  Last image of each detection event
    lastimg <- capdata[capdata$Category == "Predator",] %>%
      group_by(det_events) %>%
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last")
    #'  Last image of all OTHER detection events
    lastother <- capdata[capdata$Category == "Other",] %>%
      group_by(det_events) %>%
      slice_tail() %>%
      ungroup() %>%
      mutate(Det_type = "last")
    #'  Merge last image of each predator/other and first image of each predator
    firstlast_img <- rbind(firstimg, lastimg, lastother)
    return(firstlast_img)
  }
  firstlast_img <- lapply(eoe_5min_list, first_last_image)
  
  #'  Detections with only one image get duplicated b/c its the first & last image
  #'  Identify duplicate images (ignoring last column where they differ) and filter
  #'  to just one image per detection event
  dups <- firstlast_img[[1]] %>%
    group_by_at(vars(-Det_type)) %>% 
    filter(n() > 1) %>%
    filter(Det_type == "last")
  #'  Remove the duplicate images from the larger data set & arrange by date/time/location
  firstlast_img_20s <- anti_join(firstlast_img[[1]], dups) %>%
    arrange(CameraLocation, DateTime)
  
  
  
  