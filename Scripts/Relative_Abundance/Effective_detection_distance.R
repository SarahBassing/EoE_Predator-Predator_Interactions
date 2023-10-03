  #'  -------------------------------
  #'  Effective Detection Distance
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  October 2023
  #'  -------------------------------
  #'  Summarize effective detection distance (EDD) for each predator and camera setup 
  #'  (predator vs ungulate) based on measuring distance to animal (m) in motion 
  #'  triggered camera trap images. Append EDD to time-in-front-of-camera (TIFC) 
  #'  measurements for each camera, species, and year, filling in missing observations 
  #'  with species and camera setup-specific average EDDs.
  #'  -------------------------------
  
  #'  Load libraries
  library(data.table)
  library(tidyverse)
  
  #'  ----------------------
  ####  Read & format data  ####
  #'  ----------------------
  #'  Distance measurements by camera setup type
  pred_mt <- read_csv("./Data/Relative abundance data/RAI Phase 2/Pred_MT_distance.csv") %>%
    filter(ReferencePics == FALSE) %>%
    filter(Species != "marten") %>%
    mutate(CamID = ifelse(File == "EOE2022S_IDFG2508_20220725_171547_MD_1.JPG", "IDFG2508", CamID),
           Species = ifelse(Species == "black_bear", "bear_black", Species)) 
  
  ung_mt <- read_csv("./Data/Relative abundance data/RAI Phase 2/Ung_MT_distance.csv") %>%
    filter(ReferencePics == FALSE) %>%
    filter(Species != "whitetailed_deer") %>%
    mutate(Species = ifelse(Species == "black_bear", "bear_black", Species))
  
  #'  Camera location data
  cams <- read_csv("./Data/IDFG camera data/cams_eoe_long_Smr2020-2022.csv") %>%
    dplyr::select(c("NewLocationID", "CamID", "Setup", "Season", "Lat", "Long", "CameraHeight_M", "CameraFacing")) %>%
    filter(Season != "Wtr20") %>%
    filter(NewLocationID != "UNKNOWN") %>%
    filter(CamID != "unknown") %>%
    group_by(NewLocationID) %>%
    arrange(Season) %>%
    slice(n()) %>%
    ungroup() 
  
  #'  Calculated time-in-front-of-camera for each species and camera 
  load("./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time.RData")
  
  format_tifc <- function(tifc) {
    skinny_tifc <- tifc %>%
      filter(common_name == "bear_black" | common_name == "bobcat" | common_name == "coyote" | 
             common_name == "mountain_lion" | common_name == "wolf") %>%
      mutate(Setup = ifelse(grepl("P", location), "predator", "ungulate"),
      ) %>%
      #'  Change column names to match distance and camera location data
      rename("NewLocationID" = "location") %>%
      rename("Species" = "common_name")
    return(skinny_tifc)
  }
  tifc_skinny <- lapply(tt_list, format_tifc)
  
  #'  Hang onto non-predator species data
  prey_tifc <- function(tifc) {
    tifc_no_pred <- tifc %>%
      filter(common_name == "elk" | common_name == "moose" | common_name == "muledeer" | 
               common_name == "rabbit_hare" | common_name == "whitetaileddeer")
    return(tifc_no_pred)
  }
  tifc_prey <- lapply(tt_list, prey_tifc)
  
  #'  ---------------------------------
  ####  Format and summarize distance  ####
  #'  ---------------------------------
  #'  Format so each detected individual per image has a unique row where its
  #'  unique detection distance is recorded
  detection_dis_per_individual <- function(cam_dat, dist_dat) {
    #'  Bind camera and distance data together
    individual_distances <- full_join(cam_dat, dist_dat, by = "CamID") %>%
      dplyr::select(c("NewLocationID", "CamID", "Setup", "Season", "Lat", "Long", "File", "CameraHeight_M", "CameraFacing",
                      "DateTime", "Species", "Count", "Ind1_dist", "Ind2_dist", "Ind3_dist", "Ind4_dist")) %>%
      #'  Replace Count NAs with 1 (required for uncount function below to work)
      mutate(Count = ifelse(is.na(Count), 1, Count)) %>%
      #'  Repeat rows based on count data (i.e., if 2 animals detected, repeat observation once)
      uncount(Count, .remove = FALSE) %>%
      mutate(Count = ifelse(is.na(Species), NA, Count)) %>%
      #'  Group and number each observation within group data so repeat rows have a unique identifier
      group_by(NewLocationID, File, DateTime) %>%
      mutate(ind_animal = row_number()) %>%
      ungroup() %>%
      #'  Create column where unique detection distance is recorded for each individual 
      #'  based on how many times each detection was repeated and the number of 
      #'  individual distances measured for that image
      mutate(det_dist = ifelse(ind_animal == 1, Ind1_dist, Ind2_dist),
             det_dist = ifelse(ind_animal == 3, Ind3_dist, det_dist),
             det_dist = ifelse(ind_animal == 4, Ind4_dist, det_dist)) %>%
      dplyr::select(-c(File, Ind1_dist, Ind2_dist, Ind3_dist, Ind4_dist))
    return(individual_distances)
  }
  
  pred_dist <- detection_dis_per_individual(cam_dat = cams[cams$Setup == "predator",], dist_dat = pred_mt)
  ung_dist <- detection_dis_per_individual(cam_dat = cams[cams$Setup == "ungulate",], dist_dat = ung_mt)
  
  #'  Summarize detection distances
  #'  Note- warnings due to infinite values being created when finding min/max of NAs
  avg_dist <- function(dist) {
    #'  Summarize detection distances per species per camera site
    avg_dist_per_spp_cam <- dist %>%
      group_by(NewLocationID, Species) %>%
      summarise(avg_dist = mean(det_dist, na.rm = TRUE),
                min_dist = min(det_dist, na.rm = TRUE),
                max_dist = max(det_dist, na.rm = TRUE),
                n_obs = n()) %>%
      ungroup() %>%
      mutate(avg_dist = round(avg_dist, 2), 
             avg_dist = ifelse(is.na(avg_dist), NA, avg_dist),
             min_dist = ifelse(is.infinite(min_dist), NA, min_dist),
             max_dist = ifelse(is.infinite(max_dist), NA, max_dist)) 
      
    #'  Summarize detection distances per species
    avg_dist_per_spp <- dist %>%
      filter(!is.na(Species)) %>%
      group_by(Species) %>%
      summarise(avg_dist = mean(det_dist, na.rm = TRUE),
                se_dist = sd(det_dist, na.rm = TRUE)/sqrt(nrow(.)),
                min_dist = min(det_dist, na.rm = TRUE),
                max_dist = max(det_dist, na.rm = TRUE),
                n_obs = n()) %>%
      ungroup() %>%
      mutate(avg_dist = round(avg_dist, 2),
             se_dist = round(se_dist, 3))
    print(as.data.frame(avg_dist_per_spp))
    
    avg_dist <- list(avg_dist_per_spp_cam, avg_dist_per_spp)
    return(avg_dist)
  }
  pred_dist_avg <- avg_dist(pred_dist)
  ung_dist_avg <- avg_dist(ung_dist)
  
  #'  Generate single data frame containing TIFC and detection distance data
  tifc_edd <- function(tifc, dist, setup) {
    #'  Bind detection distance data to TIFC data
    #'  Note- can have distances for species at sites where total_duration = 0 in TIFC df
    #'  b/c distance measurements include observations from outside July 1 - Sept 15 time frame
    merge_dat <- left_join(tifc[tifc$Setup == setup,], dist[[1]], by = c("NewLocationID", "Species"))
    
    #'  Identify which TIFC observations are missing detection data and fill in 
    #'  with species averages. Replacing observations where mean detection distance = 0
    #'  (owing to bears investigating camera) b/c messes up density calculations
    tifc_with_dist <- merge_dat %>%
      filter(!is.na(avg_dist)) %>%
      filter(avg_dist != 0)
    
    tifc_missing_dist <- merge_dat %>%
      filter(is.na(avg_dist) | avg_dist == 0) %>%
      dplyr::select(-c(avg_dist, min_dist, max_dist, n_obs)) %>%
      left_join(dist[[2]], by = "Species") %>%
      dplyr::select(-se_dist)
    
    #'  Bind so each species- & site-specific TIFC observation has detection distance data
    tifc_edd_full <- rbind(tifc_with_dist, tifc_missing_dist) %>%
      #'  Re-organize to data are ready for Calculate_density.R script
      arrange(NewLocationID, Species) %>%
      dplyr::select(-Setup) %>%
      #'  Change column names back to original TIFC names
      rename("location" = "NewLocationID") %>%
      rename("common_name" = "Species")
    
    return(tifc_edd_full)
  }
  pred_tifc_edd <- lapply(tifc_skinny, tifc_edd, dist = pred_dist_avg, setup = "predator")
  ung_tifc_edd <- lapply(tifc_skinny, tifc_edd, dist = ung_dist_avg, setup = "ungulate")
  
  #'  Merge predator and ungulate camera setups back together for final dataset
  all_together_now <- function(pred, ung) {
    full_dat <- bind_rows(pred, ung) %>%
      arrange(location, common_name)
    return(full_dat)
  }
  tifc_edd_final <- mapply(all_together_now, pred = pred_tifc_edd, ung = ung_tifc_edd, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  #'  Rename and save for use in Calculate_density.R script
  tt_list <- tifc_edd_final
  save(tt_list, file = "./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time_edd.RData")
  
  
