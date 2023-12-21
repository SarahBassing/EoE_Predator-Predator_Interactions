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
  #'  TIFC density estimates generated in Calculate_density.R (Relative_Abundance folder)
  #'  -------------------------------
  
  #'  Load libraries
  library(data.table)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  --------------------------------
  ####  Read & format detection data  ####
  #'  --------------------------------
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  #'  Sampling effort data (number of days cameras operation)
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData") 
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe22s.RData") 
  
  #' #'  -----------------------------------------------------------------
  #' ####  Initial cleaning to correct known species misidentifications   ####
  #' #'  -----------------------------------------------------------------
  #' #'  Detection data (motion trigger observations only)
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID.RData")
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID.RData")
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID.RData")
  #' 
  #' #'  Load corrected species ID and update larger datasets
  #' newSppID <- read.csv("./Data/IDFG camera data/questionable_images_doublecheck_SBBupdated.csv") %>%
  #'   mutate(posix_date_time = as.POSIXct(posix_date_time, format="%m/%d/%Y %H:%M", tz="UTC"), 
  #'          Year = as.numeric(format(posix_date_time, "%Y"))) %>%
  #'   distinct() %>%
  #'   dplyr::select(c("CamID", "NewLocationID", "File", "Location_Relative_Project", "posix_date_time", "Species", "NewSpecies", "Year"))
  #' newnewSppID <- read.csv("./Data/IDFG camera data/More questionable IDs_ST_SBB.csv") %>% 
  #'   mutate(posix_date_time = gsub("<a0>", " ", posix_date_time),
  #'          posix_date_time = as.POSIXct(posix_date_time, format="%Y-%m-%d %H:%M", tz="UTC"),
  #'          Year = as.numeric(format(posix_date_time, "%Y"))) %>%
  #'   dplyr::select(c("CamID", "NewLocationID", "File", "Location_Relative_Project", "posix_date_time", "Species", "NewSpecies", "Year"))
  #' newSppID <- rbind(newSppID, newnewSppID) %>% distinct() 
  #' #'  Filter by year and only observations with misidentified species
  #' newSppID_20s <- filter(newSppID, Year == "2020") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  #' newSppID_21s <- filter(newSppID, Year == "2021") %>% dplyr::select(-Year) %>% filter(Species != NewSpecies)
  #' #'  Snag unique ID for each image
  #' change_sppID_20s <- as.vector(newSppID_20s$Location_Relative_Project)
  #' change_sppID_21s <- as.vector(newSppID_21s$Location_Relative_Project)
  #' 
  #' #'  Function to correct species misclassifications
  #' remove_obs <- function(dets, prob_images, correctSppID) {
  #'   #'  Grab all rows that have a misID's species
  #'   obs_to_remove <- dets[dets$Location_Relative_Project %in% prob_images,]
  #'   #'  Switch original species ID with correct one
  #'   newspp <- left_join(obs_to_remove, correctSppID, by = c("CamID", "NewLocationID", "File", "Species")) %>%
  #'     dplyr::select(-c("Location_Relative_Project.y", "posix_date_time.y")) %>%
  #'     relocate(NewSpecies, .after = Species) %>%
  #'     dplyr::select(-Species) %>%
  #'     rename(Species = "NewSpecies", Location_Relative_Project = "Location_Relative_Project.x", posix_date_time = "posix_date_time.x")
  #'   #'  Remove these rows from the larger df
  #'   dets_reduced <- anti_join(dets, obs_to_remove)
  #'   #'  Replace with observation containing correct species ID
  #'   dets_updated <- bind_rows(dets_reduced, newspp) %>%
  #'     #'  Arrange everything back in order
  #'     arrange(NewLocationID, posix_date_time, File)
  #'   return(dets_updated)
  #' }
  #' eoe20s_allM_new <- remove_obs(eoe20s_allM, prob_images = change_sppID_20s, correctSppID = newSppID_20s)
  #' eoe21s_allM_new <- remove_obs(eoe21s_allM, prob_images = change_sppID_21s, correctSppID = newSppID_21s)
  #' 
  #' #'  Re-assign
  #' eoe20s_allM <- eoe20s_allM_new
  #' eoe21s_allM <- eoe21s_allM_new
  #' 
  #' #'  Save cleaned data
  #' save(eoe20s_allM, file = paste0("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_", Sys.Date(), ".RData"))
  #' save(eoe21s_allM, file = paste0("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_", Sys.Date(), ".RData"))
  #' save(eoe22s_allM, file = paste0("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_", Sys.Date(), ".RData"))
  
  #'  ------------------------
  ####  Filter detection data  ####
  #'  ------------------------
  #'  Load cleaned data (with corrected species IDs)
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_2023-09-26.RData") #2023-09-26 is most recent version
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_2023-09-26.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_2023-09-26.RData")
  
  #'  1) Remove sequential problem images from larger image set
  thin_detections <- function(dets, seqprobs) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    return(skinny_dets)
  }
  eoe20s_allM_skinny <- thin_detections(eoe20s_allM, eoe_seqprob_20s) 
  eoe21s_allM_skinny <- thin_detections(eoe21s_allM, eoe_seqprob_21s) 
  eoe22s_allM_skinny <- thin_detections(eoe22s_allM, eoe_seqprob_21s) 
  
  #'  2) Filter detection data to time period and species of interest
  #'  Species of interest are ungulates, humans/vehicles, and domestic animals
  #'  Time periods of interest determined by IDFG and plotting histograms of when 
  #'  most cameras were active each season (Detection_data_cleaning.R script)
  detections <- function(dat, start_date, end_date, days_operable) {
    dets <- dat %>%
      dplyr::select("NewLocationID", "CamID", "File", "Location_Relative_Project", 
                    "Date", "Time", "posix_date_time", "TriggerMode",
                    "OpState", "Species", "Vehicle", "Count") %>%
      filter(Species != "none") %>%
      mutate(
        Date = as.Date(Date, format = "%Y-%m-%d"), #format = "%d-%b-%Y"
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
  eoe22s_dets <- detections(eoe22s_allM_skinny, start_date = "2022-07-01", end_date = "2022-09-15", days_operable = effort_22s)
  
  #'  Double check all cameras have sampling effort data
  unique(eoe20s_dets$NewLocationID[is.na(eoe20s_dets$ndays)])
  unique(eoe21s_dets$NewLocationID[is.na(eoe21s_dets$ndays)])
  unique(eoe22s_dets$NewLocationID[is.na(eoe22s_dets$ndays)])
  
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
  eoe22s_5min_dets <- unique_detections(eoe22s_dets, elapsed_time = 300)
  #'  List 5min-elapsed detection events
  eoe_5min_list <- list(eoe20s_5min_dets, eoe21s_5min_dets ,eoe22s_5min_dets)
  
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
             RAI_nhrs = round((n_det_hrs/ndays) , 3)) %>%
      #'  Retain data from cameras where camera was operable for 30 days or more
      filter(ndays >= 30) %>%
      #'  Remove sites where camera was operable 30+ days but no wildlife detections occurred
      filter(!is.na(n_imgs))
    
    return(count_dat)
  }
  eoe20s_RAI <- RAI_metrics(eoe20s_5min_dets, effort_20s)
  eoe21s_RAI <- RAI_metrics(eoe21s_5min_dets, effort_21s)
  eoe22s_RAI <- RAI_metrics(eoe22s_5min_dets, effort_22s)
  eoe_RAI <- list(eoe20s_RAI, eoe21s_RAI, eoe22s_RAI)
  names(eoe_RAI) <- c("eoe20s_RAI", "eoe21s_RAI", "eoe22s_RAI")
  
  #'  Save
  save(eoe_RAI, file = "./Data/Relative abundance data/RAI Phase 2/eoe_RAI.RData")
  
  #'  Summary stats for RAI
  rai_stats <- function(rai) {
    summary_dat <- rai %>%
      #'  Add GMU to dataset
      mutate(gmu = sub("_.*", "", NewLocationID),
             setup = ifelse(grepl("P", NewLocationID), "predator", "ungulate")) %>%
      relocate(gmu, .after = "NewLocationID") %>%
      relocate(setup, .after = "gmu") %>%
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" | Species == "mountain_lion" | Species == "wolf") %>%
      group_by(gmu, Species) %>%
      summarise(mean_RAI = round(mean(RAI_nhrs, na.rm = TRUE), 3),
                se_RAI = round((sd(RAI_nhrs, na.rm = TRUE)/sqrt(nrow(.))), 3)) %>%
      ungroup() %>%
      arrange(gmu, mean_RAI)
    print(summary_dat)
    return(summary_dat)
  }
  rai_summary_stats <- lapply(eoe_RAI, rai_stats)
  rai_summary_stats_20s <- rai_summary_stats[[1]] %>% mutate(Year = "2020")
  rai_summary_stats_21s <- rai_summary_stats[[2]] %>% mutate(Year = "2021")
  rai_summary_stats_22s <- rai_summary_stats[[3]] %>% mutate(Year = "2022")
  rai_summary_stats_all <- rbind(rai_summary_stats_20s, rai_summary_stats_21s, rai_summary_stats_22s)
  
  write_csv(rai_summary_stats_all, "./Data/Relative abundance data/RAI Phase 2/RAI_summary_stats.csv")
  
  
  #####  Visualize relative abundance indices  #####
  #'  -----------------------------------------
  #'  Map relative density data per species, study area and year
  library(sf)
  library(ggplot2)
  library(patchwork)
  
  #'  Load spatial data
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  cams_20s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_20s_wgs84.shp")
  cams_21s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_21s_wgs84.shp")
  cams_22s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_22s_wgs84.shp")
  
  #'  List camera spatial data
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84)
  
  #'  Load and reorganize RAI data by species and year
  load("./Data/Relative abundance data/RAI Phase 2/eoe_RAI.RData")
  RAI_bear <- list(eoe_RAI[[1]][eoe_RAI[[1]]$Species == "bear_black",], eoe_RAI[[2]][eoe_RAI[[2]]$Species == "bear_black",], eoe_RAI[[3]][eoe_RAI[[3]]$Species == "bear_black",])
  RAI_bob <- list(eoe_RAI[[1]][eoe_RAI[[1]]$Species == "bobcat",], eoe_RAI[[2]][eoe_RAI[[2]]$Species == "bobcat",], eoe_RAI[[3]][eoe_RAI[[3]]$Species == "bobcat",])
  RAI_coy <- list(eoe_RAI[[1]][eoe_RAI[[1]]$Species == "coyote",], eoe_RAI[[2]][eoe_RAI[[2]]$Species == "coyote",], eoe_RAI[[3]][eoe_RAI[[3]]$Species == "coyote",])
  RAI_lion <- list(eoe_RAI[[1]][eoe_RAI[[1]]$Species == "mountain_lion",], eoe_RAI[[2]][eoe_RAI[[2]]$Species == "mountain_lion",], eoe_RAI[[3]][eoe_RAI[[3]]$Species == "mountain_lion",])
  RAI_wolf <- list(eoe_RAI[[1]][eoe_RAI[[1]]$Species == "wolf",], eoe_RAI[[2]][eoe_RAI[[2]]$Species == "wolf",], eoe_RAI[[3]][eoe_RAI[[3]]$Species == "wolf",])

  #'  Append RAI to spatial data
  spatial_rai <- function(rai, spp, cams) {
    #'  Filter data to single species
    single_spp_rai <- rai %>%
      filter(Species == spp) %>%
      #'  Rename camera location column to match spatial data
      rename("NwLctID" = "NewLocationID") 
    
    #'  Join spatial data with rai data
    rai_shp <- full_join(cams, single_spp_rai, by = "NwLctID") %>%
      mutate(Species = ifelse(is.na(Species), spp, Species),
             RAI_nimgs = ifelse(is.na(RAI_nimgs), 0, RAI_nimgs),
             RAI_ndets = ifelse(is.na(RAI_ndets), 0, RAI_ndets),
             RAI_nhrs = ifelse(is.na(RAI_nhrs), 0, RAI_nhrs)) %>%
      #'  Change camera location column back to something less awkward
      rename("NewLocationID" = "NwLctID") 
    
    return(rai_shp)
  }
  spatial_rai_bear <- mapply(rai = RAI_bear, spatial_rai, spp = "bear_black", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rai_bob <- mapply(rai = RAI_bob, spatial_rai, spp = "bobcat", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rai_coy <- mapply(rai = RAI_coy, spatial_rai, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rai_lion <- mapply(rai = RAI_lion, spatial_rai, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rai_wolf <- mapply(rai = RAI_wolf, spatial_rai, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  List spatial RAI data and save
  spatial_RAI_list <- list(spatial_rai_bear, spatial_rai_bob, spatial_rai_coy, spatial_rai_lion, spatial_rai_wolf)
  save(spatial_RAI_list, file = "./Shapefiles/IDFG spatial data/Camera_locations/spatial_RAI_list.RData")
  
  year_list <- list("2020", "2021", "2022")
  
  #'  Function to map RAI
  map_rai <- function(sf_rai, yr, spp) {
    #'  Create column for GMU
    sf_rai <- sf_rai %>%
      mutate(locs = NewLocationID) %>%
      separate(NewLocationID, c("GMU", "Setup", "site"), sep = "_") %>%
      dplyr::select(-site) %>%
      rename("NewLocationID" = "locs")
    
    #'  Filter spatial RN data by study area
    sf_rai_gmu1 <- sf_rai[sf_rai$GMU == "GMU1",]
    sf_rai_gmu6 <- sf_rai[sf_rai$GMU == "GMU6",]
    sf_rai_gmu10a <- sf_rai[sf_rai$GMU == "GMU10A",]
    
    #'  GMU 1 plot
    gmu1_rai <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "1",], fill = NA) +
      geom_sf(data = sf_rai_gmu1, aes(size = RAI_nhrs), shape  = 21, 
              col = "darkred", fill = "darkred", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "Average hourly \ndetections per day", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 6 plot
    gmu6_rai <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "6",], fill = NA) +
      geom_sf(data = sf_rai_gmu6, aes(size = RAI_nhrs), shape  = 21, 
              col = "darkgreen", fill = "darkgreen", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "Average hourly \ndetections per day", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 10A plot
    gmu10a_rai <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
      geom_sf(data = sf_rai_gmu10a, aes(size = RAI_nhrs), shape = 21, 
              col = "darkblue", fill = "darkblue", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "Average hourly \ndetections per day", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  Plot each map
    print(gmu1_rai); print(gmu6_rai); print(gmu10a_rai)
    
    #'  List GMU RAI maps together
    gmu_maps <- list(gmu1_rai, gmu6_rai, gmu10a_rai)
    
    return(gmu_maps)
  }
  rai_maps_bear <- mapply(map_rai, sf_rai = spatial_rai_bear, yr = year_list, spp = "Black bear", SIMPLIFY = FALSE)
  rai_maps_bob <- mapply(map_rai, sf_rai = spatial_rai_bob, yr = year_list, spp = "Bobcat", SIMPLIFY = FALSE)
  rai_maps_coy <- mapply(map_rai, sf_rai = spatial_rai_coy, yr = year_list, spp = "Coyote", SIMPLIFY = FALSE)
  rai_maps_lion <- mapply(map_rai, sf_rai = spatial_rai_lion, yr = year_list, spp = "Mountain lion", SIMPLIFY = FALSE)
  rai_maps_wolf <- mapply(map_rai, sf_rai = spatial_rai_wolf, yr = year_list, spp = "wolf", SIMPLIFY = FALSE)
  
  #'  Plot same species and study area back to back across years
  gmu_by_yr_plots <- function(fig, spp) {
    #'  Note: list order is [[i]][[j]] 
    #'  where i = 2020, 2021, or 2022 and j = GMU1, GMU6, or GMU10a 
    #'  GMU 1 plots
    gmu1_patch <- fig[[1]][[1]] + fig[[2]][[1]] + fig[[3]][[1]] +
      plot_annotation(paste("GMU 1", spp, "relative abundance index (detection rate)"))
    
    #'  GMU 6 plots
    gmu6_patch <- fig[[1]][[2]] + fig[[2]][[2]] + fig[[3]][[2]] +
      plot_annotation(paste("GMU 6", spp, "relative abundance index (detection rate)"))
    
    #'  GMU 10A plots
    gmu10a_patch <- fig[[1]][[3]] + fig[[2]][[3]] + fig[[3]][[3]] +
      plot_annotation(paste("GMU 10A", spp, "relative abundance index (detection rate)"))
    
    #'  Print figure panels
    print(gmu1_patch); print(gmu6_patch); print(gmu10a_patch)
    
    #'  List
    gmu_patchwork_list <- list(gmu1_patch, gmu6_patch, gmu10a_patch)
    
    return(gmu_patchwork_list)
  }
  rai_gmu_bear <- gmu_by_yr_plots(rai_maps_bear, spp = "black bear")
  rai_gmu_bob <- gmu_by_yr_plots(rai_maps_bob, spp = "bobcat")
  rai_gmu_coy <- gmu_by_yr_plots(rai_maps_coy, spp = "coyote")
  rai_gmu_lion <- gmu_by_yr_plots(rai_maps_lion, spp = "mountain lion")
  rai_gmu_wolf <- gmu_by_yr_plots(rai_maps_wolf, spp = "wolf")
  
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu1_bear.tiff", rai_gmu_bear[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu6_bear.tiff", rai_gmu_bear[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu10A_bear.tiff", rai_gmu_bear[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu1_bob.tiff", rai_gmu_bob[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu6_bob.tiff", rai_gmu_bob[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu10A_bob.tiff", rai_gmu_bob[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu1_coy.tiff", rai_gmu_coy[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu6_coy.tiff", rai_gmu_coy[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu10A_coy.tiff", rai_gmu_coy[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu1_lion.tiff", rai_gmu_lion[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu6_lion.tiff", rai_gmu_lion[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu10A_lion.tiff", rai_gmu_lion[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu1_wolf.tiff", rai_gmu_wolf[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu6_wolf.tiff", rai_gmu_wolf[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/Figures/RAI_gmu10A_wolf.tiff", rai_gmu_wolf[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")

    
  #'  --------------------
  ####  Correlation test  ####
  #'  --------------------
  #'  Load all relative abundance metrics
  load("./Data/Relative abundance data/RAI Phase 2/eoe_RAI.RData")
  tifc_density <- read_csv("./Data/Relative abundance data/RAI Phase 2/eoe_all_years_density_edd_predonly.csv")
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  
  #'  Test for correlation between different RAI metrics
  compare_counts <- function(dets, tifc, rn) {
    tifc <- dplyr::select(tifc, c("location", "common_name", "cpue", "cpue_km2", "cpue_100km2")) 
    names(tifc) <- c("NewLocationID", "Species", "cpue", "cpue_km2", "cpue_100km2")
    
    rn <- dplyr::select(rn, -c("GMU", "Setup", "season", "RN.sd"))
    
    all_RAIs <- dets %>%
      #'  Bind tifc density measure to larger RAI data set
      left_join(tifc, by = c("NewLocationID", "Species")) %>%
      left_join(rn, by = c("NewLocationID", "Species")) %>%
      #'  Filter out observations of non-focal species included in RAI but not TIFC
      filter(!is.na(cpue_100km2)) %>%
      #'  Reduce to species of interest and remove sites with all NAs
      filter(!is.na(RAI_nimgs)) %>%
      # filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
      #          Species == "elk" | Species == "human" | Species == "moose" |
      #          Species == "mountain_lion" | Species == "muledeer" | Species == "rabbit_hare" |
      #          Species == "whitetaileddeer" | Species == "wolf" | Species == "cattle_cow")
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "mountain_lion" | Species == "wolf")
    
    #'  Make sure there are no sites with missing TIFC data 
    print(unique(all_RAIs$NewLocationID[is.na(all_RAIs$cpue_100km2)]))
    
    #'  Calculate correlation coefficient for all combinations
    pearsons_cor <- all_RAIs %>%
      dplyr::select(-NewLocationID) %>%
      group_by(Species) %>%
      #'  Calculate correlation coefficient for each pairwise combo of counts
      summarize(img_hrs = round(cor(RAI_nimgs, RAI_nhrs), 3),
                img_tifc = round(cor(RAI_nimgs, cpue_km2), 3), 
                hrs_tifc = round(cor(RAI_nhrs, cpue_km2), 3),
                img_rn = round(cor(RAI_nimgs, RN.n), 3),
                hrs_rn = round(cor(RAI_nhrs, RN.n), 3),
                tifc_rn = round(cor(cpue_km2, RN.n), 3)) %>%
      ungroup()
    print(pearsons_cor)
    
    #'  List full RAI & TIFC density estimates, along with correlation test
    RAI_TIFC <- list(all_RAIs, pearsons_cor)
    
    return(RAI_TIFC)
  }
  eoe20s_abund <- compare_counts(eoe_RAI[[1]], tifc = tifc_density[tifc_density$season == "Smr20",], rn = RN_abundance[[1]]) 
  eoe21s_abund <- compare_counts(eoe_RAI[[2]], tifc = tifc_density[tifc_density$season == "Smr21",], rn = RN_abundance[[2]]) 
  eoe22s_abund <- compare_counts(eoe_RAI[[3]], tifc = tifc_density[tifc_density$season == "Smr22",], rn = RN_abundance[[3]]) 
  
  #'  List all relative abundance indices per year
  all_RAI_metrics <- list(eoe20s_abund[[1]], eoe21s_abund[[1]], eoe22s_abund[[1]])
  
  #'  Save
  save(all_RAI_metrics, file = "./Data/Relative abundance data/RAI Phase 2/all_RAI_metrics_edd_predonly.RData")
  
  #'  -----------------------------------
  #####  Visualize mean RAI, TIFC, & RN  #####
  #'  -----------------------------------
  #'  Read in summary data and reformat to be consistent
  #'  Mean detection rate
  mean_DR <- read_csv("./Data/Relative abundance data/RAI Phase 2/RAI_summary_stats.csv") %>%
    mutate(RAI_method = "DR")
  names(mean_DR) <- c("GMU", "Species", "Mean", "Variability", "Year", "RAI_method")
  
  #'  TIFC mean density
  mean_TIFC <- read_csv("./Data/Relative abundance data/RAI Phase 2/tifc_density_stats_avg_edd_predonly.csv") %>%
    mutate(season = ifelse(season == "Smr20", "2020", season),
           season = ifelse(season == "Smr21", "2021", season),
           season = ifelse(season == "Smr22", "2022", season),
           season = as.numeric(season),
           RAI_method = "TIFC")
  names(mean_TIFC) <- c("Year", "GMU", "Species", "mean_density_km2", "se_density_km2", "Mean", "Variability", "RAI_method")
  
  #'  RN model GMU-wide average DENSITY estimates
  mean_RN_gmu1 <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu1.csv") %>% mutate(GMU = "GMU1")
  mean_RN_gmu6 <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu6.csv") %>% mutate(GMU = "GMU6")
  mean_RN_gmu10A <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu10A.csv") %>% mutate(GMU = "GMU10A")
  mean_RNden <- rbind(mean_RN_gmu1, mean_RN_gmu6, mean_RN_gmu10A) %>% arrange(Year, GMU, Species) %>%
    mutate(Species = ifelse(Species == "black bear", "bear_black", Species),
           Species = ifelse(Species == "mountain lion", "mountain_lion", Species), 
           RAI_method = "RN model")
  names(mean_RNden) <- c("Species", "Year", "Mean", "SD", "GMU", "RAI_method")
  
  #'  RN model local ABUNDANCE estimates (to be averaged)
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  mean_RNn_20s <- RN_abundance[[1]] %>%
    group_by(GMU, Species) %>%
    summarise(mean_RN_abundance = round(mean(RN.n), 3),
              se_RN_abundance = round((sd(RN.n, na.rm = TRUE)/sqrt(nrow(.))), 3)) %>%
    ungroup() %>%
    mutate(Year = as.numeric("2020"), 
           RAI_method = "RN model")
  mean_RNn_21s <- RN_abundance[[2]] %>%
    group_by(GMU, Species) %>%
    summarise(mean_RN_abundance = round(mean(RN.n), 3),
              se_RN_abundance = round((sd(RN.n, na.rm = TRUE)/sqrt(nrow(.))), 3)) %>%
    ungroup() %>%
    mutate(Year = as.numeric("2021"),
           RAI_method = "RN model")
  mean_RNn_22s <- RN_abundance[[3]] %>%
    group_by(GMU, Species) %>%
    summarise(mean_RN_abundance = round(mean(RN.n), 3),
              se_RN_abundance = round((sd(RN.n, na.rm = TRUE)/sqrt(nrow(.))), 3)) %>%
    ungroup() %>%
    mutate(Year = as.numeric("2022"), 
           RAI_method = "RN model")
  mean_RNn <- rbind(mean_RNn_20s, mean_RNn_21s, mean_RNn_22s)
  names(mean_RNn) <- c("GMU", "Species", "Mean", "Variability", "Year", "RAI_method")
  
  #'  Bind summary data together
  all_metrics <- full_join(mean_DR, mean_TIFC, by = c("Year", "GMU", "Species", "RAI_method", "Mean", "Variability")) %>%
    full_join(mean_RNn, by = c("Year", "GMU", "Species", "RAI_method", "Mean", "Variability")) %>%
    mutate(Year = as.factor(as.character(Year))) %>%
    relocate(Year, .before = GMU) %>%
    relocate(RAI_method, .after = "Species") %>%
    dplyr::select(-c("mean_density_km2", "se_density_km2"))
  
  #'  Facet_grid labels
  RAI_names <- c("Detection rate", "RN abundance", "TIFC density")
  names(RAI_names) <- c("DR", "RN model", "TIFC")
  
  #'  Visualize differences in DR, RN abundance, and TIFC density
  RAI_comparison_plot <- ggplot(all_metrics, aes(x = Species, y = Mean, color = GMU)) +
    geom_errorbar(aes(ymin = Mean-Variability, ymax = Mean+Variability, width = 0.3), position = position_dodge(width=0.9)) +
    geom_point(stat = "identity", position = position_dodge(width=0.9), size = 2) +
    facet_grid(RAI_method ~ Year, scales = "free", labeller = labeller(RAI_method = RAI_names)) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Predator relative abundance indices") +
    ylab("Mean RAI")
  RAI_comparison_plot_gmu1 <- ggplot(all_metrics[all_metrics$GMU == "GMU1",], aes(x = Species, y = Mean)) + 
    geom_errorbar(aes(ymin = Mean-Variability, ymax = Mean+Variability, width = 0.3), position = position_dodge(width=0.9), color = "#D55E00") +
    geom_point(stat = "identity", position = position_dodge(width=0.9), size = 2, color = "#D55E00") +
    facet_grid(RAI_method ~ Year, scales = "free", labeller = labeller(RAI_method = RAI_names)) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Predator relative abundance indices, GMU1")
  RAI_comparison_plot_gmu6 <- ggplot(all_metrics[all_metrics$GMU == "GMU6",], aes(x = Species, y = Mean)) + 
    geom_errorbar(aes(ymin = Mean-Variability, ymax = Mean+Variability, width = 0.3), position = position_dodge(width=0.9), color = "#0072B2") +
    geom_point(stat = "identity", position = position_dodge(width=0.9), size = 2, color = "#0072B2") +
    facet_grid(RAI_method ~ Year, scales = "free", labeller = labeller(RAI_method = RAI_names)) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Predator relative abundance indices, GMU6")
  RAI_comparison_plot_gmu10a <- ggplot(all_metrics[all_metrics$GMU == "GMU10A",], aes(x = Species, y = Mean)) + 
    geom_errorbar(aes(ymin = Mean-Variability, ymax = Mean+Variability, width = 0.3), position = position_dodge(width=0.9), color = "#CC79A7") +
    geom_point(stat = "identity", position = position_dodge(width=0.9), size = 2, color = "#CC79A7") +
    facet_grid(RAI_method ~ Year, scales = "free", labeller = labeller(RAI_method = RAI_names)) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Predator relative abundance indices, GMU10A")
  
  ggsave("./Outputs/Figures/RAI_comparison_plot.tiff", RAI_comparison_plot, units = "in",
         width = 12, height = 14, dpi = 600, device = "tiff", compression = 'lzw')
  ggsave("./Outputs/Figures/RAI_comparison_plot_gmu1.tiff", RAI_comparison_plot_gmu1, units = "in",
         width = 10, height = 14, dpi = 600, device = "tiff", compression = 'lzw')
  ggsave("./Outputs/Figures/RAI_comparison_plot_gmu6.tiff", RAI_comparison_plot_gmu6, units = "in",
         width = 10, height = 14, dpi = 600, device = "tiff", compression = 'lzw')
  ggsave("./Outputs/Figures/RAI_comparison_plot_gmu10a.tiff", RAI_comparison_plot_gmu10a, units = "in",
         width = 10, height = 14, dpi = 600, device = "tiff", compression = 'lzw')
  
  