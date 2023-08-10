  #'  --------------------------------
  #'  Time In Front of Camera (TIFC)
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  August 2023
  #'  --------------------------------
  #'  Calculate relative abundance at each site using the Time In Front of Camera
  #'  (TIFC) approach (Huggard et al. 2018, Warbington & Boyce 2020, Becker et al. 2022).
  #'  Code adapted from Becker et al. 2018 repository: https://github.com/mabecker89/tifc-method/tree/main
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R (Data_Formatting folder)
  #'  Sampling effort data generated in Detection_histories_for_occmod.R (MultiSpp_OccMod folder)
  #'  -------------------------------
  
  #'  Load libraries
  library(data.table)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  --------------------------------
  ####  Read & format detection data  ####
  #'  --------------------------------
  #'  Detection data (motion trigger observations only)
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_2023-08-09.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_2023-08-09.RData")
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
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
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s) 
  
  #'  2) Filter  data to time period and species of interest
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "File", "Location_Relative_Project", 
                    "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      filter(Species != "none") %>%
      filter(Vehicle != "TRUE") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time)) %>%
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "elk" | Species == "human" | Species == "moose" |
               Species == "mountain_lion" | Species == "muledeer" | Species == "rabbit_hare" |
               Species == "whitetaileddeer" | Species == "wolf" | Species == "cattle_cow") %>%
      #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
      mutate(Count = ifelse(Count == 0, 1, Count)) %>%
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
  eoe21s_dets <- detections(eoe21s_allM_skinny, start_date = "2021-07-01", end_date = "2021-09-15")
  eoe_dets <- list(eoe20s_dets, eoe21s_dets)
  
  #'  ---------------------------
  ####  TIME IN FRONT of CAMERA  ####
  #'  ---------------------------
  #'  Calculate Time In Front of Camera (TIFC) 
  #'  Function identifies series of sequential images and gaps between sequential
  #'  images, adjusts series by the probability an animal left the FoV, and 
  #'  calculates the time in field-of-view (Ta) for each species and camera.
  #'  
  #'  Code adapted from Becker et al. 2018 repository
  #'  https://github.com/mabecker89/tifc-method/blob/main/R/base/01_probabilistic-gaps.R
  #'  ---------------------------
  tifc <- function(dets) {
    
    #'  Identify series of consecutive images of the same species and flag gaps 
    #'  in sequential images where an animal may have left & returned OR remains 
    #'  but did not trigger camera for short period of time (gaps <= 120 sec).
    #'  "Series" defined as consecutive images with intervals <120 seconds between
    #'  consecutive images. Intermediate image without the species ends the series.
    series_and_gaps <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Identify image series and gaps between sequential images
      mutate(#'  Lagged time stamp
             previous_img_time = lag(posix_date_time),
             #'  Lead time stamp
             next_img_time = lead(posix_date_time),
             #'  Calculate difference in time between current and previous image
             time_diff = as.numeric(posix_date_time - previous_img_time),
             #'  Calculate difference in time between current and next image
             time_diff_next = as.numeric(abs(posix_date_time - next_img_time)),
              #'  Flag whether species differs in next image
             previous_species = lag(Species),
             spp_diff = ifelse(Species != previous_species, TRUE, FALSE),
             #'  Flag whether data are from a different camera site
             previous_site_loc = lag(NewLocationID),
             site_diff = ifelse(NewLocationID != previous_site_loc, TRUE, FALSE),
             #'  Flag gaps between sequential images that need to be evaluated
             gap_check = ifelse(site_diff == FALSE & spp_diff == FALSE & (time_diff <= 120 & time_diff >=120), 1, 0),
             #'  Identify new series of images based on a change in location, species, or time gap > 120 s
             new_series = ifelse(site_diff ==  TRUE | spp_diff == TRUE | time_diff > 120, 1, 0),
             series_ID = c(0, cumsum(new_series[-1]))) %>%
      group_by(series_ID) %>%
      #'  Force time_diff to 0 if start of new series & time_diff_next to 0 if end of a series
      mutate(time_diff = ifelse(row_number() == 1, 0, time_diff),
             time_diff_next = ifelse(row_number() == n(), 0, time_diff_next)) %>%
      ungroup()
    #'  Eventually adjust time differences by probability of leaving if gaps 20 - 120 sec
    
    #'  Calculate time between sequential images per species
    time_btwn_imgs <- series_and_gaps %>%
      mutate(series_ID_previous = lag(series_ID)) %>%
      #'  Drop first image from each series b/c time_diff = 0
      filter(series_ID == series_ID_previous) %>%
      group_by(Species) %>%
      #'  Calculate average time between images and number of images per species
      summarise(time_btwn_imgs = mean(time_diff),
                sample_size = n())
    
    #'  Calculate total time in front of the camera by series (tts)
    total_time_per_series <- series_and_gaps %>%
      left_join(time_btwn_imgs, by = "Species") %>%
      group_by(series_ID) %>%
      mutate(#'  Check if image is first or last in series
        bookend = ifelse(row_number() == 1 | row_number() == n(), 1, 0),
        
        #'  Calculate time for each image            ####  UPDATE time_diff & time_diff_next WITH PROB. LEAVE ADJUSTED TIMES ONCE CALCULATED  ####
        image_time = ifelse(bookend == 1,
                            #'  Sum time elapsed between sequential images
                            #'  If first/last image in series, add times between images t-1 to t+1, then divide by two, 
                            #'  and add half the average time btwn images for a given species 
                            #'  to account for time before/after first/last image to account for animal entering/leaving trigger zone
                            ((time_diff + time_diff_next) / 2) + (time_btwn_imgs / 2),
                            #'  Else, add times between images t-1 to t+1, then divide by two
                            (time_diff + time_diff_next) / 2),
        
        #'  Multiply image time by the number of animals present
        image_time_ni = image_time * Count) %>%      ####  NEED TO FIGURE OUT HOW TO GIVE IT AVERAGE # ANIMALS PER SERIES  ####
      
      #'  Group by species and add to df
      group_by(Species, .add = TRUE) %>%
      #'  Calculate total time and number of images for each series
      summarise(n_images = n(),
                series_total_time = sum(image_time_ni)) %>%
      ungroup()
    
    #'  Calculate total time in front of camera per site and species
    total_time <- series_and_gaps %>%
      group_by(series_ID) %>%
      arrange(posix_date_time, .by_group = TRUE) %>%
      #'  Snag first image per series
      filter(row_number() == 1) %>%
      #'  Bind with total time in front of camera
      left_join(total_time_per_series, by = c("series_ID", "Species")) %>%
      #'  Reduce to only necessary columns
      select(NewLocationID, posix_date_time, Species, series_total_time) %>%
      ungroup() %>%
      #'  Bind average time between images for each species
      left_join(time_btwn_imgs, by = "Species") %>%
      relocate(NewLocationID, .before = series_ID) %>%
      relocate(posix_date_time, .after = NewLocationID) %>%
      rename("avg_time_btwn_imgs_spp" = "time_btwn_imgs") %>%
      rename("total_imgs_spp" = "sample_size")
    
    #'  List all data frames together
    tifc_list <- list(series_and_gaps, time_btwn_imgs, total_time_per_series, total_time)
    
    return(tifc_list)
  }
  eoe_total_time_in_FoV <- lapply(eoe_dets, tifc)
  checkit <- eoe_total_time_in_FoV[[2]][[4]]
  
  #'  Save
  save(eoe_total_time_in_FoV, file = "./Data/Time_In_Front_of_Camera/eoe_total_time_in_FoV.RData")
  
  
  
  
  
  
  