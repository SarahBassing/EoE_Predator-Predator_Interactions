  #'  -------------------------------
  #'  Relative Abundance Indices
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  July 2023
  #'  -------------------------------
  #'  Format detection data and generate relative abundance indices based on 
  #'  different counts of detection data from camera traps. 
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R (Data_Formatting folder)
  #'  Sampling effort data generated in Detection_histories_for_occmod.R (MultiSpp_OccMod folder)
  #'  TIFC density estimates generated in Time_In_Front_of_Camera.R (Relative_Abundance folder)
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
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  Sampling effort data (number of days cameras operation)
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData") 
  
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
  
  #' #'  Save cleaned data
  #' save(eoe20s_allM, file = paste0("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_", Sys.Date(), ".RData"))
  #' save(eoe21s_allM, file = paste0("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_", Sys.Date(), ".RData"))
  
  #'  ------------------------
  ####  Filter detection data  ####
  #'  ------------------------
  #'  1) Remove sequential problem images from larger image set
  thin_detections <- function(dets, seqprobs) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    return(skinny_dets)
  }
  eoe20s_allM_skinny <- thin_detections(eoe20s_allM, eoe_seqprob_20s) #eoe20s_allM_2023-08-09 is most recent version
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s) #eoe21s_allM_2023-08-09 is most recent version
  
  #'  2) Filter detection data to time period and species of interest
  #'  Species of interest are ungulates, humans/vehicles, and domestic animals
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dets, start_date, end_date, days_operable) {
    dets <- dets %>%
      dplyr::select("NewLocationID", "CamID", "File", "Location_Relative_Project", 
                    "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      filter(Species != "none") %>%
      mutate(
        Date = as.Date(Date, format = "%d-%b-%Y"),
        Time = chron(times = Time),
        Species = ifelse(Vehicle == "TRUE", "human_motorized", Species),
        Category = ifelse(Species == "bobcat" | Species == "bear_black" | 
                            Species == "coyote" | Species == "mountain_lion" | 
                            Species == "wolf", "Predator", "Other")) %>% 
      #'  Remove observations of "unknown deer species" - not including these in analyses
      #'  AND often one image of "unknown" deer mixed in with sequential images
      #'  of a known deer species, skewing detection rates & RAI
      filter(Species != "deer_speciesunknown") %>%
      dplyr::select(-Vehicle) %>%
      #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
      mutate(Count = ifelse(Count == 0, 1, Count),
             Count = ifelse(is.na(Count), 1, Count)) %>%
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      #'  Remove maintenance images so not included in human/motorized image sets
      filter(OpState != "maintenance") %>%
      arrange(NewLocationID)
    
    #'  Add sampling effort to detections data frame
    effort <- dplyr::select(days_operable, c(NewLocationID, ndays))
    dets <- left_join(dets, effort, by = "NewLocationID") %>%
      #'  Retain data from cameras where camera was operable for 30 days or more
      filter(ndays >= 30)
    
    return(dets)
  }
  eoe20s_dets <- detections(eoe20s_allM_skinny, start_date = "2020-07-01", end_date = "2020-09-15", days_operable = effort_20s)
  eoe21s_dets <- detections(eoe21s_allM_skinny, start_date = "2021-07-01", end_date = "2021-09-15", days_operable = effort_21s)
  
  #'  Double check all cameras have sampling effort data
  unique(eoe20s_dets$NewLocationID[is.na(eoe20s_dets$ndays)])
  unique(eoe21s_dets$NewLocationID[is.na(eoe21s_dets$ndays)])
  
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
  eoe21s_5min_dets <- unique_detections(eoe21s_dets, elapsed_time = 300)
  #'  List 5min-elapsed detection events
  eoe_5min_list <- list(eoe20s_5min_dets, eoe21s_5min_dets)
  
  #'  -------------------------------------------------  
  ####  Calculating simple Relative Abundance Indices  ####
  #'  -------------------------------------------------
  #'  Multiple approaches to count animal detections which are then used to 
  #'  generate different metrics representing relative abundance indices
  #'  RAI = (total counts/trap nights) (Palmer et al. 2018, O'Brien et al. 2003)
  
  RAI_metrics <- function(dets, effort) {
    #'  Pull out number of operable camera trap days from sampling effort data
    trap_nights <- dplyr::select(effort, c(NewLocationID, ndays)) %>%
      #'  A few cameras were operable 0.5 days which makes for funky RAI values
      #'  based on how they're calculated below so rounding up to 1 day
      mutate(ndays = ifelse(ndays < 1.0, 1.0, ndays))
    
    #'  Count 1: total number of images
    #'  -------------------------------  
    #'  Sum all images of a given species per camera site (Palmer et al. 2018)
    all_imgs <- dets %>%
      group_by(NewLocationID, Species) %>%
      summarize(n_imgs = n()) %>%
      ungroup()
    
    #'  Count 2: total number of unique detection events
    #'  ------------------------------------------------
    #'  Sum unique detection events per species per camera site (Palmer et al. 2018)
    all_dets <- dets %>%
        group_by(NewLocationID, Species, det_events) %>%
        slice(1L) %>%
        ungroup() %>%
        group_by(NewLocationID, Species) %>%
        summarize(n_dets = n()) %>%
        ungroup()
    
    #'  Count 3: total number of hourly detection events
    #'  ------------------------------------------------
    #'  Sum number of hours where at least 1 detection event occurred per species 
    #'  per camera site (Aubsand et al. 2023)
    sum_det_hr <- dets %>%
      #'  Retain single image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup() %>%
      #'  Floor time to the hour of each image
      mutate(floor_DT = floor_date(posix_date_time, unit = "hour")) %>%
      #'  Thin to just the first image per hour for a given species per camera
      group_by(NewLocationID, Species, floor_DT) %>%
      slice(1L) %>%
      ungroup() %>%
      group_by(NewLocationID, Species) %>%
      #'  Count total number of hours each species was detected per camera
      summarize(n_det_hrs = n()) %>%
      ungroup()
    
    #'  Merge all counts into single data frame and calculate RAI metrics
    count_dat <- full_join(all_imgs, all_dets, by = c("NewLocationID", "Species")) %>%
      full_join(sum_det_hr, by = c("NewLocationID", "Species")) %>%
      full_join(trap_nights, by = "NewLocationID") %>%
      mutate(RAI_nimgs = round((n_imgs/ndays) , 3),
             RAI_ndets = round((n_dets/ndays) , 3),
             RAI_nhrs = round((n_det_hrs/ndays) , 3))
    
    return(count_dat)
  }
  eoe20s_RAI <- RAI_metrics(eoe20s_5min_dets, effort_20s)
  eoe21s_RAI <- RAI_metrics(eoe21s_5min_dets, effort_21s)
  eoe_RAI <- list(eoe20s_RAI, eoe21s_RAI)
  
  #'  Save
  save(eoe_RAI, file = "./Data/Relative abundance data/RAI Phase 2/eoe_RAI.RData")
  
  #'  --------------------
  ####  Correlation test  ####
  #'  --------------------
  #'  Read in TIFC density estimates for comparison
  load("./Outputs/Relative_Abundance/TIFC/eoe_density_list.RData")
  
  #'  Test for correlation between different RAI metrics
  compare_counts <- function(dets, tifc) {
    tifc <- dplyr::select(tifc, c("NewLocationID", "Species", "density_km2")) 
    
    all_RAIs <- dets %>%
      #'  Bind tifc density measure to larger RAI data set
      left_join(tifc, by = c("NewLocationID", "Species")) %>%
      #'  Reduce to species of interest and remove sites with all NAs
      filter(!is.na(RAI_nimgs)) %>%
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "elk" | Species == "human" | Species == "moose" |
               Species == "mountain_lion" | Species == "muledeer" | Species == "rabbit_hare" |
               Species == "whitetaileddeer" | Species == "wolf" | Species == "cattle_cow") 
    
    #'  Make sure there are no sites with missing TIFC data 
    print(unique(all_RAIs$NewLocationID[is.na(all_RAIs$density_km2)]))
    
    #'  Calculate correlation coefficient for all combinations
    pearsons_cor <- all_RAIs %>%
      dplyr::select(-NewLocationID) %>%
      group_by(Species) %>%
      #'  Calculate correlation coefficient for each pairwise combo of counts
      summarize(img_dets = round(cor(RAI_nimgs, RAI_ndets), 3), #, use = "complete.obs"
                img_hrs = round(cor(RAI_nimgs, RAI_nhrs), 3),
                dets_hrs = round(cor(RAI_nhrs, RAI_ndets), 3),
                img_tifc = round(cor(RAI_nimgs, density_km2), 3),
                dets_tifc = round(cor(RAI_ndets, density_km2), 3),
                hrs_tifc = round(cor(RAI_nhrs, density_km2), 3)) %>%
      ungroup()
    print(pearsons_cor)
    return(pearsons_cor)
  }
  eoe20s_corr <- compare_counts(eoe20s_RAI, eoe_density_list[[1]])
  eoe21s_corr <- compare_counts(eoe21s_RAI, eoe_density_list[[2]])
  
  
  
  
  