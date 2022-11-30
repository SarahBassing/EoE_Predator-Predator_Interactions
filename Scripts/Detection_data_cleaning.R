  #'  ---------------------------------
  #'  Cleaning raw detection data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean camera detection data, make sure annual data sets are 
  #'  consistent, identify cameras with potential problems (non-operational or 
  #'  missing data), and filter into more manageable data sets.
  #'  
  #'  Camera deployment data prepared in Camera_cleaning.R
  #'  Massive detection data sets separated/thinned in Camera_split_raw_data.R
  #'  ---------------------------------
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(sf)
  library(raster)
  library(ggplot2)
  library(tidyverse)
  library(lubridate)
  library(chron)
  
  #'  Camera locations
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(-X)
  cams_wolf_long <- read.csv("./Data/IDFG camera data/cams_wolf_long.csv") %>%
    dplyr::select(-X)
  
  #' #'  Full detection data sets - MASSIVE so only load 1-2 at a time
  #' #'  Motion triggered images
  #' load("./Data/IDFG camera data/Split datasets/eoe20s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe20w_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe21s_allM.RData")
  #' 
  #' load("./Data/IDFG camera data/Split datasets/wolf19s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf20s_allM.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf21s_allM.RData")
  #' 
  #' #'  Time triggered images
  #' load("./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  #' 
  #' load("./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  #'  Keeper data sets (all animal/human/vehicle images & noon timelapse images)
  #'  Much easier to manage...
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  #'  Add camera setup info (P/U/O/A) and consistent naming structure to match 
  #'  camera location data
  #'  EoE cameras
  eoe_deploy_info <- function(dets, season, pred) { 
    #'  Filter to specific season and predator setup cameras
    sub_cams <- cams_eoe_long %>%
      filter(Season == season) %>%
      filter(Setup == pred)
    #'  Identify which observations came from predator setup cameras (TRUE)
    pcams <- dets$CamID %in% sub_cams$CamID
    #'  Add column to detection data
    dets$pcams <- pcams
    #'  Create new camera location name that matches camera deployment data
    dets <- dets %>%
      mutate(Setup = ifelse(pcams == "TRUE", "P", "U"),
             NewLocationID = paste0("GMU", Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-pcams) %>%
      relocate(NewLocationID, .before = File) %>%
      relocate(Setup, .after = NewLocationID)
    return(dets)
  }
  #'  Run keeper eoe data sets through function
  eoe_seasons <- list("Smr20", "Wtr20", "Smr21")
  eoe_motion_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_motion_skinny, pred = "predator", SIMPLIFY = FALSE) 
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  eoe_motion_list[[3]] <- dplyr::select(eoe_motion_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  #'  Save for later use
  # save(eoe_motion_list, file = "./Data/IDFG camera data/Split datasets/eoe_motion_skinny_NewLocationID.RData")
  
  eoe_noon_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  eoe_noon_list[[3]] <- dplyr::select(eoe_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  #'  Save for later use
  eoe_time_list <- eoe_noon_list
  # save(eoe_time_list, file = "./Data/IDFG camera data/Split datasets/eoe_time_skinny_NewLocationID.RData")
  
  #'  Double check it worked
  eoe21s_noon <- eoe_noon_list[[3]]
  
  #'  Run full eoe motion trigger data sets through function
  #' eoe_motion_smr20 <- eoe_deploy_info(dets = eoe20s_allM, season = "Smr20", pred = "predator")
  #' eoe_motion_wtr20 <- eoe_deploy_info(dets = eoe20w_allM, season = "Wtr20", pred = "predator")
  #' eoe_motion_smr21 <- eoe_deploy_info(dets = eoe21s_allM, season = "Smr21", pred = "predator") %>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' eoe20s_allM <- eoe_deploy_info(dets = eoe20s_allM, season = "Smr20", pred = "predator")
  #' eoe20w_allM <- eoe_deploy_info(dets = eoe20w_allM, season = "Wtr20", pred = "predator")
  #' eoe21s_allM <- eoe_deploy_info(dets = eoe21s_allM, season = "Smr21", pred = "predator")
  #' eoe21s_allM <- mutate(eoe21s_allM, NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' save(eoe20s_allM, file = "./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  #' save(eoe20w_allM, file = "./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  #' save(eoe21s_allM, file = "./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  
  #' #'  Run full eoe timelapse data sets through function
  #' eoe_time_smr20 <- eoe_deploy_info(dets = eoe20s_allT, season = "Smr20", pred = "predator")
  #' eoe_time_wtr20 <- eoe_deploy_info(dets = eoe20w_allT, season = "Wtr20", pred = "predator")
  #' eoe_time_smr21 <- eoe_deploy_info(dets = eoe21s_allT, season = "Smr21", pred = "predator") %>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' eoe20s_allT <- eoe_deploy_info(season = "Smr20", dets = eoe20s_allT, pred = "predator")
  #' eoe20w_allT <- eoe_deploy_info(season = "Wtr20", dets = eoe20w_allT, pred = "predator")
  #' eoe21s_allT <- eoe_deploy_info(season = "Smr21", dets = eoe21s_allT, pred = "predator")
  #' eoe21s_allT <- mutate(eoe21s_allT, NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' save(eoe20s_allT, file = "./Data/IDFG camera data/Split datasets/eoe20s_allT_NewLocationID.RData")
  #' save(eoe20w_allT, file = "./Data/IDFG camera data/Split datasets/eoe20w_allT_NewLocationID.RData")
  #' save(eoe21s_allT, file = "./Data/IDFG camera data/Split datasets/eoe21s_allT_NewLocationID.RData")
  
  
  #'  Wolf cameras
  wolf_deploy_info <- function(dets, season, abund, abund_occu) { 
    #'  Filter to specific season and Abundance cameras
    sub_abund_cams <- cams_wolf_long %>%
      filter(Season == season) %>%
      filter(Target == abund)
    #'  Identify which observations came from Abundance cameras (TRUE)
    acams <- dets$CamID %in% sub_abund_cams$CamID
    #'  Add column to detection data
    dets$acams <- acams
    #'  Filter to specific season and Abund_Occu cameras
    sub_abundoccu_cams <- cams_wolf_long %>%
      filter(Season == season) %>%
      filter(Target == abund_occu)
    #'  Identify which observations came from Abund_Occu cameras (TRUE)
    aocams <- dets$CamID %in% sub_abundoccu_cams$CamID
    #'  Add column to detection data
    dets$aocams <- aocams
    #'  Create new camera location name that matches camera deployment data
    dets <- dets %>%
      mutate(Setup = ifelse(acams == "TRUE", "A", "O"),
             Setup = ifelse(aocams == "TRUE", "B", Setup),
             NewLocationID = paste0(Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-c(acams, aocams)) %>%
      relocate(NewLocationID, .before = File) %>%
      relocate(Setup, .after = NewLocationID)
    return(dets)
  }
  #'  Run keeper wolf data sets through function
  wolf_seasons <- list("Smr19", "Smr20", "Smr21")
  wolf_motion_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_motion_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)
  #'  Fix NewLocationID info- this was recorded differently for Smr20 & Smr21 
  #'  images in original data set so unnecessary info gets added in this function - 
  #'  need to revert back to original information
  wolf_motion_list[[2]] <- mutate(wolf_motion_list[[2]], NewLocationID = paste0("GMU", NewLocationID))
  wolf_motion_list[[3]] <- dplyr::select(wolf_motion_list[[3]], -c(NewLocationID)) %>%
    mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
           NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  #'  Save for later use
  # save(wolf_motion_list, file = "./Data/IDFG camera data/Split datasets/wolf_motion_skinny_NewLocationID.RData")
  
  wolf_noon_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_time_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)  
  wolf_noon_list[[2]] <- mutate(wolf_noon_list[[2]], NewLocationID = paste0("GMU", NewLocationID))
  wolf_noon_list[[3]] <- dplyr::select(wolf_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
           NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  wolf_motion_list[[3]] <- dplyr::select(wolf_motion_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  wolf_noon_list[[3]] <- dplyr::select(wolf_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  #'  Save for later use
  wolf_time_list <- wolf_noon_list
  # save(wolf_time_list, file = "./Data/IDFG camera data/Split datasets/wolf_time_skinny_NewLocationID.RData")
  
  #'  Double check it worked
  wolf21s_noon <- wolf_noon_list[[3]]
  
  #' #'  Run full wolf timelapse data sets through function
  #' wolf_motion_smr20 <- wolf_deploy_info(dets = wolf19s_allM, season = "Smr19", abund = "Abundance", abund_occu = "Abund_Occu")
  #' wolf_motion_wtr20 <- wolf_deploy_info(dets = wolf20s_allM, season = "Smr20", abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   mutate(NewLocationID = paste0("GMU", NewLocationID))
  #' wolf_motion_smr21 <- wolf_deploy_info(dets = wolf21s_allM, season = "Smr21", abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
  #'          NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' wolf19s_allM <- wolf_deploy_info(season = "Smr19", dets = wolf19s_allM, abund = "Abundance", abund_occu = "Abund_Occu")
  #' wolf20s_allM <- wolf_deploy_info(season = "Smr20", dets = wolf20s_allM, abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   mutate(NewLocationID = paste0("GMU", NewLocationID))
  #' wolf21s_allM <- wolf_deploy_info(season = "Smr21", dets = wolf21s_allM, abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
  #'          NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' save(wolf19s_allM, file = "./Data/IDFG camera data/Split datasets/wolf19s_allM_NewLocationID.RData")
  #' save(wolf20s_allM, file = "./Data/IDFG camera data/Split datasets/wolf20s_allM_NewLocationID.RData")
  #' save(wolf21s_allM, file = "./Data/IDFG camera data/Split datasets/wolf21s_allM_NewLocationID.RData")
  
  #' #'  Run full wolf timelapse data sets through function
  #' wolf_time_smr20 <- wolf_deploy_info(dets = wolf19s_allT, season = "Smr19", abund = "Abundance", abund_occu = "Abund_Occu")
  #' wolf_time_wtr20 <- wolf_deploy_info(dets = wolf20s_allT, season = "Smr20", abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   mutate(NewLocationID = paste0("GMU", NewLocationID))
  #' wolf_time_smr21 <- wolf_deploy_info(dets = wolf21s_allT, season = "Smr21", abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
  #'          NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' wolf19s_allT <- wolf_deploy_info(season = "Smr19", dets = wolf19s_allT, abund = "Abundance", abund_occu = "Abund_Occu")
  #' wolf20s_allT <- wolf_deploy_info(season = "Smr20", dets = wolf20s_allT, abund = "Abundance", abund_occu = "Abund_Occu") %>%
  #'   mutate(NewLocationID = paste0("GMU", NewLocationID))
  #' wolf21s_allT <- wolf_deploy_info(season = "Smr21", dets = wolf21s_allT, abund = "Abundance", abund_occu = "Abund_Occu")%>%
  #'   dplyr::select(-c(NewLocationID)) %>%
  #'   mutate(LocationID = paste0("GMU", Gmu, "_", LocationID),
  #'          NewLocationID = LocationID) %>%
  #'   relocate(NewLocationID, .after = LocationID)
  #' save(wolf19s_allT, file = "./Data/IDFG camera data/Split datasets/wolf19s_allT_NewLocationID.RData")
  #' save(wolf20s_allT, file = "./Data/IDFG camera data/Split datasets/wolf20s_allT_NewLocationID.RData")
  #' save(wolf21s_allT, file = "./Data/IDFG camera data/Split datasets/wolf21s_allT_NewLocationID.RData")
  
  
  ####  Set time zone  ####
  #'  ------------------
  #'  Detection data was recorded in MST (America/Edmonton (UTC-07:00); tz="America/Edmonton")
  set_tzone <- function(dat) {
    tz(dat$posix_date_time) <- "America/Edmonton"
    print(attr(dat$posix_date_time, "tzone"))
    print(range(dat$posix_date_time, na.rm = T))
    # dat$Date <- as.Date(dat$Date, format = "%d-%b-%Y")
    return(dat)
  }
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny_NewLocationID.RData")
  eoe_motion_list <- lapply(eoe_motion_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny_NewLocationID.RData")
  eoe_noon_list <- lapply(eoe_time_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny_NewLocationID.RData")
  wolf_motion_list <- lapply(wolf_motion_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny_NewLocationID.RData")
  wolf_noon_list <- lapply(wolf_time_list, set_tzone)
  
  #' #'  Set time zone on full data sets
  #' load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  #' eoe20s_allM <- set_tzone(eoe20s_allM)
  #' load("./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  #' eoe20w_allM <- set_tzone(eoe20w_allM)
  #' load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  #' eoe21s_allM <- set_tzone(eoe21s_allM)
  
  load("./Data/IDFG camera data/Split datasets/eoe20s_allT_NewLocationID.RData")
  eoe20s_allT <- set_tzone(eoe20s_allT)
  load("./Data/IDFG camera data/Split datasets/eoe20w_allT_NewLocationID.RData")
  eoe20w_allT <- set_tzone(eoe20w_allT)
  load("./Data/IDFG camera data/Split datasets/eoe21s_allT_NewLocationID.RData")
  eoe21s_allT <- set_tzone(eoe21s_allT)
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allM_NewLocationID.RData")
  wolf19s_allM <- set_tzone(wolf19s_allM)
  load("./Data/IDFG camera data/Split datasets/wolf20s_allM_NewLocationID.RData")
  wolf20s_allM <- set_tzone(wolf20s_allM)
  load("./Data/IDFG camera data/Split datasets/wolf21s_allM_NewLocationID.RData")
  wolf21s_allM <- set_tzone(wolf21s_allM)
  
  #'  NOTE only cameras set for abundance monitoring took time-trigger images -
  #'  this is only half the cameras!
  load("./Data/IDFG camera data/Split datasets/wolf19s_allT_NewLocationID.RData")
  wolf19s_allT <- set_tzone(wolf19s_allT)
  load("./Data/IDFG camera data/Split datasets/wolf20s_allT_NewLocationID.RData")
  wolf20s_allT <- set_tzone(wolf20s_allT)
  load("./Data/IDFG camera data/Split datasets/wolf21s_allT_NewLocationID.RData")
  wolf21s_allT <- set_tzone(wolf21s_allT)
  
  #' eoe_motion_smr20 <- set_tzone(eoe_motion_smr20)
  #' eoe_motion_wtr20 <- set_tzone(eoe_motion_wtr20)
  #' eoe_motion_smr21 <- set_tzone(eoe_motion_smr21)
  #' 
  #' eoe_time_smr20 <- set_tzone(eoe_time_smr20)
  #' eoe_time_wtr20 <- set_tzone(eoe_time_wtr20)
  #' eoe_time_smr21 <- set_tzone(eoe_time_smr21)
  #' 
  #' wolf_motion_smr20 <- set_tzone(wolf_motion_smr20)
  #' wolf_motion_wtr20 <- set_tzone(wolf_motion_wtr20)
  #' wolf_motion_smr21 <- set_tzone(wolf_motion_smr21)
  #' 
  #' wolf_time_smr20 <- set_tzone(wolf_time_smr20)
  #' wolf_time_wtr20 <- set_tzone(wolf_time_wtr20)
  #' wolf_time_smr21 <- set_tzone(wolf_time_smr21)

  
  ####  Visualize observations over time  ####
  #'  ------------------------------------
  #'  https://r4ds.had.co.nz/dates-and-times.html
  #'  Note that when you use date-times in a numeric context (like in a histogram), 
  #'  1 means 1 second, so a binwidth of 86400 means one day. For dates, 1 means 1 day.
  #'  Noon timelapse images only- represents camera function, not animal activity/detections
  EoE_noon_cams_Smr20 <- eoe_noon_list[[1]] %>%  #  Summer 2020
    filter(posix_date_time < ymd(20210101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) + # 86400 seconds = 1 day
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-09-01")), linetype = 4, color = "blue") +
    ggtitle("Number of active EoE cameras, Summer 2020")
  EoE_noon_cams_Smr20
  ggsave("./Outputs/Figures/EoE_noon_cams_Smr20.png", EoE_noon_cams_Smr20, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  EoE_noon_cams_Wtr2021 <- eoe_noon_list[[2]] %>%  #  Winter 2020-2021
    filter(posix_date_time > ymd(20201101) & posix_date_time < ymd(20210531)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-12-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-02-01")), linetype = 4, color = "blue") +
    ggtitle("Number of active EoE cameras, Winter 2020-2021")
  EoE_noon_cams_Wtr2021
  ggsave("./Outputs/Figures/EoE_noon_cams_Wtr2021.png", EoE_noon_cams_Wtr2021, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  EoE_noon_cams_Smr21 <- eoe_noon_list[[3]] %>%  #  Summer 2021
    filter(posix_date_time > ymd(20210501) & posix_date_time < ymd(20211101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-09-01")), linetype = 4, color = "blue") +
    ggtitle("Number of active EoE cameras, Summer 2021")
  EoE_noon_cams_Smr21
  ggsave("./Outputs/Figures/EoE_noon_cams_Smr21.png", EoE_noon_cams_Smr21, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  Wolf_noon_cams_Smr19 <- wolf_noon_list[[1]] %>%  # Summer 2019
    filter(posix_date_time < ymd(20200101)) %>% 
    filter(Time == "12:0:0") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) + # 86400 seconds = 1 day
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2019-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2019-09-15")), linetype = 4, color = "blue") +
    ggtitle("Number of active wolf cameras, Summer 2019")
  Wolf_noon_cams_Smr19
  ggsave("./Outputs/Figures/Wolf_noon_cams_Smr19.png", Wolf_noon_cams_Smr19, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  Wolf_noon_cams_Smr20 <- wolf_noon_list[[2]] %>%  #  Summer 2020
    filter(posix_date_time < ymd(20210101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-09-15")), linetype = 4, color = "blue") +
    ggtitle("Number of active wolf cameras, Summer 2020")
  Wolf_noon_cams_Smr20
  ggsave("./Outputs/Figures/Wolf_noon_cams_Smr20.png", Wolf_noon_cams_Smr20, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  Wolf_noon_cams_Smr21 <- wolf_noon_list[[3]] %>%  #  Summer 2021
    filter(posix_date_time < ymd(20220101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-09-15")), linetype = 4, color = "blue") +
    ggtitle("Number of active wolf cameras, Summer 2021")
  Wolf_noon_cams_Smr21
  ggsave("./Outputs/Figures/Wolf_noon_cams_Smr21.png", Wolf_noon_cams_Smr21, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  
  #'  By GMU - are there temporal trends in deployment by GMU?
  EoE_noon_cams_Smr20GMU <- eoe_noon_list[[1]] %>%  #  Summer 2020
    filter(posix_date_time < ymd(20210101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time, group = Gmu)) + 
    geom_freqpoly(binwidth = 86400) + # 86400 seconds = 1 day
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-09-01")), linetype = 4, color = "blue") +
    facet_wrap(~Gmu, scales = "free_y") +
    ggtitle("Number of active EoE cameras by GMU, Summer 2020")
  EoE_noon_cams_Smr20GMU
  ggsave("./Outputs/Figures/EoE_noon_cams_Smr20GMU.png", EoE_noon_cams_Smr20GMU, units = "in", 
         width = 12, height = 6, dpi = 600, device = "png")
  
  EoE_noon_cams_Wtr2021GMU <- eoe_noon_list[[2]] %>%  #  Winter 2020-2021
    filter(posix_date_time > ymd(20201101) & posix_date_time < ymd(20210531)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time, group = Gmu)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2020-12-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-02-01")), linetype = 4, color = "blue") +
    facet_wrap(~Gmu, scales = "free_y") +
    ggtitle("Number of active EoE cameras by GMU, Winter 2020-2021")
  EoE_noon_cams_Wtr2021GMU
  ggsave("./Outputs/Figures/EoE_noon_cams_Wtr2021GMU.png", EoE_noon_cams_Wtr2021GMU, units = "in", 
         width = 12, height = 6, dpi = 600, device = "png")
  
  EoE_noon_cams_Smr21GMU <- eoe_noon_list[[3]] %>%  #  Summer 2021
    filter(posix_date_time > ymd(20210501) & posix_date_time < ymd(20211101)) %>% 
    filter(Time == "12:00:00") %>%
    filter(!is.na(Gmu)) %>%
    ggplot(aes(posix_date_time, group = Gmu)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-07-01")), linetype = 4, color = "blue") +
    geom_vline(xintercept = as.POSIXct(as.Date("2021-09-01")), linetype = 4, color = "blue") +
    facet_wrap(~Gmu, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle("Number of active EoE cameras, Summer 2021")
  EoE_noon_cams_Smr21GMU
  ggsave("./Outputs/Figures/EoE_noon_cams_Smr21GMU.png", EoE_noon_cams_Smr21GMU, units = "in", 
         width = 12, height = 6, dpi = 600, device = "png")
  
  #'  Not plotting by GMU for wolf data b/c basically state-wide... TOO MANY to plot
  
  
  ####  How many cameras were operating  ####
  #'  ------------------------------------
  #'  How many cameras were deployed each season?
  length(unique(cams_eoe_long$NewLocationID[cams_eoe_long$Season == "Smr20"]))
  length(unique(cams_eoe_long$NewLocationID[cams_eoe_long$Season == "Wtr20"]))
  length(unique(cams_eoe_long$NewLocationID[cams_eoe_long$Season == "Smr21"]))
  
  length(unique(cams_wolf_long$NewLocationID[cams_wolf_long$Season == "Smr19"]))
  length(unique(cams_wolf_long$NewLocationID[cams_wolf_long$Season == "Smr20"]))
  length(unique(cams_wolf_long$NewLocationID[cams_wolf_long$Season == "Smr21"]))
  
  #'  How many cameras were taking motion triggered images?
  length(unique(eoe_motion_list[[1]]$NewLocationID))
  length(unique(eoe_motion_list[[2]]$NewLocationID))
  length(unique(eoe_motion_list[[3]]$NewLocationID))
  
  length(unique(wolf_motion_list[[1]]$NewLocationID))
  length(unique(wolf_motion_list[[2]]$NewLocationID))
  length(unique(wolf_motion_list[[3]]$NewLocationID))
  
  #'  How many cameras were taking time triggered images?
  length(unique(eoe_noon_list[[1]]$NewLocationID))
  length(unique(eoe_noon_list[[2]]$NewLocationID))
  length(unique(eoe_noon_list[[3]]$NewLocationID))
  
  length(unique(wolf_noon_list[[1]]$NewLocationID))
  length(unique(wolf_noon_list[[2]]$NewLocationID))
  length(unique(wolf_noon_list[[3]]$NewLocationID))
  
  #'  Which cameras were taking motion but not time triggered images?
  mt_mismatch <- function(motion_cams, noon_cams) {
    mcams <- unique(motion_cams$NewLocationID)
    tcams <- unique(noon_cams$NewLocationID)
    print("Cameras with motion but no noon time triggered images")
    print(m_not_in_t <- setdiff(mcams, tcams))
    print("Cameras with noon time but no motion triggered images")
    print(t_not_in_m <- setdiff(tcams, mcams))
    mismatch_cams <- cbind(m_not_in_t, t_not_in_m)
    return(mismatch_cams)
  }
  eoe_mismatch <- mapply(mt_mismatch, motion_cams = eoe_motion_list, noon_cams = eoe_noon_list)
  wolf_mismatch <- mapply(mt_mismatch, motion_cams = wolf_motion_list, noon_cams = wolf_noon_list)
  #'  Looks like wolf occupancy cameras did not take time-triggered images which
  #'  explains big differences in # of cams taking motion vs time triggered images
  
  #' #'  Same thing but with FULL data sets, not keepers only
  #' eoe_mismatch <- mt_mismatch(motion_cams = eoe_motion_smr20, noon_cams = eoe_time_smr20)
  #' eoe_mismatch <- mt_mismatch(motion_cams = eoe_motion_wtr20, noon_cams = eoe_time_wtr20)
  #' eoe_mismatch <- mt_mismatch(motion_cams = eoe_motion_smr21, noon_cams = eoe_time_smr21)
  #' 
  #' wolf_mismatch <- mt_mismatch(motion_cams = wolf19s_allM, noon_cams = wolf19s_allT)
  #' wolf_mismatch <- mt_mismatch(motion_cams = wolf20s_allM, noon_cams = wolf20s_allT)
  #' wolf_mismatch <- mt_mismatch(motion_cams = wolf21s_allM, noon_cams = wolf21s_allT)
  
  #'  Review motion triggered images
  eoe_m_smr20 <- eoe_motion_list[[1]]
  # eoe_m_smr20 <- eoe_motion_smr20 # full data set
  prob_eoe_m_smr20 <- eoe_m_smr20[eoe_m_smr20$NewLocationID == "GMU6_U_11" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_81" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_P_11" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_36" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_129" | 
                                    eoe_m_smr20$NewLocationID == "GMU10A_U_160" | 
                                    eoe_m_smr20$NewLocationID == "GMU10A_U_8",]
  eoe_m_wtr20 <- eoe_motion_list[[2]]
  # eoe_m_wtr20 <- eoe_motion_wtr20 # full data set
  prob_eoe_m_wtr20 <- eoe_m_wtr20[eoe_m_wtr20$NewLocationID == "GMU6_U_164" |
                                    eoe_m_wtr20$NewLocationID == "GMU6_P_99" |
                                    eoe_m_wtr20$NewLocationID == "GMU6_U_99" |
                                    eoe_m_wtr20$NewLocationID == "GMU6_U_36" |
                                    eoe_m_wtr20$NewLocationID == "GMU10A_U_46" |
                                    eoe_m_wtr20$NewLocationID == "GMU10A_U_87" |
                                    eoe_m_wtr20$NewLocationID == "GMU6_U_69" |
                                    eoe_m_wtr20$NewLocationID == "GMU6_U_105" |
                                    eoe_m_wtr20$NewLocationID == "GMU10A_U_11",]
  eoe_m_smr21 <- eoe_motion_list[[3]]
  # eoe_m_smr21 <- eoe_motion_smr21 # full data set
  prob_eoe_m_smr21 <- eoe_m_smr21[eoe_m_smr21$NewLocationID == "GMU6_P_19" |
                                    eoe_m_smr21$NewLocationID == "GMU10A_U_75" |
                                    eoe_m_smr21$NewLocationID == "GMU10A_U_159" |
                                    eoe_m_smr21$NewLocationID == "GMU6_P_29" |
                                    eoe_m_smr21$NewLocationID == "GMU6_P_92",]
  
  wolf_m_smr19 <- wolf_motion_list[[1]]
  # wolf_m_smr19 <- wolf19s_allM # full data set
  prob_wolf_m_smr19 <- wolf_m_smr19[wolf_m_smr19$NewLocationID == "GMU48_A_2777",]
  
  wolf_m_smr20 <- wolf_motion_list[[2]]
  # wolf_m_smr20 <- wolf20s_allM # full data set
  prob_wolf_m_smr20 <- wolf_m_smr20[wolf_m_smr20$NewLocationID == "GMU29_A_1559" |
                                      wolf_m_smr20$NewLocationID == "GMU13_B_1026" |
                                      wolf_m_smr20$NewLocationID == "GMU12_B_947" |
                                      wolf_m_smr20$NewLocationID == "GMU12_B_49" |
                                      wolf_m_smr20$NewLocationID == "GMU10A_B_45" |
                                      wolf_m_smr20$NewLocationID == "GMU12_B_54" |
                                      wolf_m_smr20$NewLocationID == "GMU8A_B_36" |
                                      wolf_m_smr20$NewLocationID == "GMU10_B_48",]
  
  wolf_m_smr21 <- wolf_motion_list[[3]]
  # wolf_m_smr21 <- wolf21s_allM # full data set
  prob_wolf_m_smr21 <- wolf_m_smr21[wolf_m_smr21$NewLocationID == "GMU3_A_310" |
                                      wolf_m_smr21$NewLocationID == "GMU4A_A_201" |
                                      wolf_m_smr21$NewLocationID == "GMU8A_A_576" |
                                      wolf_m_smr21$NewLocationID == "GMU4A_A_202" |
                                      wolf_m_smr21$NewLocationID == "GMU39_A_3101" |
                                      wolf_m_smr21$NewLocationID == "GMU59_A_2057" |
                                      wolf_m_smr21$NewLocationID == "GMU51_A_2606" |
                                      wolf_m_smr21$NewLocationID == "GMU51_A_2735" |
                                      wolf_m_smr21$NewLocationID == "GMU61_A_1995",]
  
  #'  Review noontime triggered images
  eoe_t_smr20 <- eoe_noon_list[[1]] 
  # eoe_t_smr20 <- eoe_motion_wtr20 # full data set
  prob_eoe_t_smr20 <- eoe_t_smr20[eoe_t_smr20$NewLocationID == "GMU6_U_11" | 
                                    eoe_t_smr20$NewLocationID == "GMU6_U_81" | 
                                    eoe_t_smr20$NewLocationID == "GMU6_P_11" | 
                                    eoe_t_smr20$NewLocationID == "GMU6_U_36" | 
                                    eoe_t_smr20$NewLocationID == "GMU6_U_129" | 
                                    eoe_t_smr20$NewLocationID == "GMU10A_U_160" | 
                                    eoe_t_smr20$NewLocationID == "GMU10A_U_8",]
  eoe_t_wtr20 <- eoe_noon_list[[2]]
  # eoe_t_wtr20 <- eoe_time_wtr20 # full data set
  prob_eoe_t_wtr20 <- eoe_t_wtr20[eoe_t_wtr20$NewLocationID == "GMU6_U_164" |
                                    eoe_t_wtr20$NewLocationID == "GMU6_P_99" |
                                    eoe_t_wtr20$NewLocationID == "GMU6_U_99" |
                                    eoe_t_wtr20$NewLocationID == "GMU6_U_36" |
                                    eoe_t_wtr20$NewLocationID == "GMU10A_U_46" |
                                    eoe_t_wtr20$NewLocationID == "GMU10A_U_87" |
                                    eoe_t_wtr20$NewLocationID == "GMU6_U_69" |
                                    eoe_t_wtr20$NewLocationID == "GMU6_U_105" |
                                    eoe_t_wtr20$NewLocationID == "GMU10A_U_11",]
  eoe_t_smr21 <- eoe_noon_list[[3]]
  # eoe_t_smr21 <- eoe_time_smr21 # full data set
  prob_eoe_t_smr21 <- eoe_t_smr21[eoe_t_smr21$NewLocationID == "GMU6_P_19" |
                                    eoe_t_smr21$NewLocationID == "GMU10A_U_75" |
                                    eoe_t_smr21$NewLocationID == "GMU10A_U_159" |
                                    eoe_t_smr21$NewLocationID == "GMU6_P_29" |
                                    eoe_t_smr21$NewLocationID == "GMU6_P_92",]
  
  wolf_t_smr19 <- wolf_noon_list[[1]]
  # wolf_t_smr19 <- wolf19s_allT # full data set
  prob_wolf_t_smr19 <- wolf_t_smr19[wolf_t_smr19$NewLocationID == "GMU48_A_2777",]
  
  wolf_t_smr20 <- wolf_noon_list[[2]]
  # wolf_t_smr20 <- wolf20s_allT # full data set
  prob_wolf_t_smr20 <- wolf_t_smr20[wolf_t_smr20$NewLocationID == "GMU29_A_1559" |
                                      wolf_t_smr20$NewLocationID == "GMU13_B_1026" |
                                      wolf_t_smr20$NewLocationID == "GMU12_B_947" |
                                      wolf_t_smr20$NewLocationID == "GMU12_B_49" |
                                      wolf_t_smr20$NewLocationID == "GMU10A_B_45" |
                                      wolf_t_smr20$NewLocationID == "GMU12_B_54" |
                                      wolf_t_smr20$NewLocationID == "GMU8A_B_36" |
                                      wolf_t_smr20$NewLocationID == "GMU10_B_48",]
  
  wolf_t_smr21 <- wolf_noon_list[[3]]
  # wolf_t_smr21 <- wolf21s_allT # full data set
  prob_wolf_t_smr21 <- wolf_t_smr21[wolf_t_smr21$NewLocationID == "GMU3_A_310" |
                                      wolf_t_smr21$NewLocationID == "GMU4A_A_201" |
                                      wolf_t_smr21$NewLocationID == "GMU8A_A_576" |
                                      wolf_t_smr21$NewLocationID == "GMU4A_A_202" |
                                      wolf_t_smr21$NewLocationID == "GMU39_A_3101" |
                                      wolf_t_smr21$NewLocationID == "GMU59_A_2057" |
                                      wolf_t_smr21$NewLocationID == "GMU51_A_2606" |
                                      wolf_t_smr21$NewLocationID == "GMU51_A_2735" |
                                      wolf_t_smr21$NewLocationID == "GMU61_A_1995",]
  
  #'  Notes regarding missing motion vs time trigger images recorded in "Problem cameras - notes" Excel file in Data folder
  
  
  ####  Cameras with known operational issues  ####
  #'  -----------------------------------------
  #'  Pull out any images from cameras with a noted misdirected/obscured viewshed
  #'  Looking for images labeled: "severely misdirected", "partially obscured", 
  #'  "completely obscured", "malfunction", "nightbad__dayok
  #'  Images marked "minorly misdirected" are still OK to use per S.Thompson
  #'  "partially obscured": 25-75% obscured but probably can still see some animals,
  #'  generally associated with animal messing with camera but other images from 
  #'  trigger still usable to ID species - currently assuming OK to use
  problem_children <- function(dat) {
    #'  Retain NewLocationID info of cameras with known problems
    prob_cams <- as.data.frame(unique(dat$NewLocationID[dat$OpState == "severely misdirected" |
                                                          dat$OpState == "completely obscured" |
                                                          dat$OpState == "malfunction" |
                                                          dat$OpState == "nightbad__dayok"]))
    colnames(prob_cams) <- "NewLocationID"
    #'  Retain all images from just problem cameras
    prob_pix <- semi_join(dat, prob_cams, by = "NewLocationID")
    #'  Double check I have the same number of cameras
    print(nrow(prob_cams))
    print(length(unique(prob_pix$NewLocationID)))
    return(prob_pix)
  }
  #'  Run keeper data sets through to help flag potentially problematic cameras
  eoe_m_probs <- lapply(eoe_motion_list, problem_children)
  eoe_t_probs <- lapply(eoe_noon_list, problem_children)   
  wolf_m_probs <- lapply(wolf_motion_list, problem_children)
  wolf_t_probs <- lapply(wolf_noon_list, problem_children)   
  
  #'  Focus on the full time trigger data sets - these provide a lot more information
  #'  about whether visual obstructions and misalignment are a long-term problem
  #'  that need to be addressed or short-lived and not a real issue.
  eoe_t_20s_probs <- problem_children(eoe20s_allT)
  eoe_t_20w_probs <- problem_children(eoe20w_allT)
  eoe_t_21s_probs <- problem_children(eoe21s_allT)
  
  wolf19s_all <- list(wolf19s_allM, wolf19s_allT)
  wolf19s_probs <- lapply(wolf19s_all, problem_children)
  wolf20s_all <- list(wolf20s_allM, wolf20s_allT)
  wolf20s_probs <- lapply(wolf20s_all, problem_children)
  wolf21s_all <- list(wolf21s_allM, wolf21s_allT)
  wolf21s_probs <- lapply(wolf21s_all, problem_children)
  
  # wolf_t_19s_probs <- problem_children(wolf19s_allT)
  # wolf_t_20s_probs <- problem_children(wolf20s_allT)
  # wolf_t_21s_probs <- problem_children(wolf21s_allT)
  

  #'  Filter images to series where camera was obscured or misdirected for 1+ hour
  #'  Using this to represent days when camera wasn't full operational
  sequential_probs <- function(dat) {
    #'  If camera's view is completely obscured
    cond1 <- expr(dat$OpState == "completely obscured")
    bad_view_pix1 <- dat %>% 
      mutate(problem_pix = 
               rep(rle(!!cond1)$values & rle(!!cond1)$lengths >= 6, 
                   rle(!!cond1)$lengths)) %>%
      filter(problem_pix)
    #'  If camera's angle is severely misdirected from initial deployment
    cond2 <- expr(dat$OpState == "severely misdirected")
    bad_view_pix2 <- dat %>% 
      mutate(problem_pix = 
               rep(rle(!!cond2)$values & rle(!!cond2)$lengths >= 6, 
                   rle(!!cond2)$lengths)) %>%
      filter(problem_pix)
    #'  If camera's view is ok during the day but compromised at night
    cond3 <- expr(dat$OpState == "nightbad__dayok")
    bad_view_pix3 <- dat %>% 
      mutate(problem_pix = 
               rep(rle(!!cond3)$values & rle(!!cond3)$lengths >= 6, 
                   rle(!!cond3)$lengths)) %>%
      filter(problem_pix)
    #'  Merge all images with problems
    bad_view_pix <- rbind(bad_view_pix1, bad_view_pix2, bad_view_pix3) %>%
      arrange(NewLocationID, posix_date_time)
    return(bad_view_pix)
  }
  #'  Flag problematic images from each full time-trigger data set
  eoe_1hr_20s <- sequential_probs(eoe_t_20s_probs)
  eoe_1hr_20w <- sequential_probs(eoe_t_20w_probs)
  eoe_1hr_21s <- sequential_probs(eoe_t_21s_probs)
  
  wolf_1hr_19s <- lapply(wolf19s_probs, sequential_probs)
  wolf_1hr_20s <- lapply(wolf20s_probs, sequential_probs)
  wolf_1hr_21s <- lapply(wolf21s_probs, sequential_probs)
  
  # wolf_1hr_19s <- sequential_probs(wolf_t_19s_probs)
  # wolf_1hr_20s <- sequential_probs(wolf_t_20s_probs)
  # wolf_1hr_21s <- sequential_probs(wolf_t_21s_probs)
  
  
  #'  Pull out problem dates based on images when camera is obscured/misdirected
  #'  for 1+ hour on a given day
  prob_days <- function(dat) {
    prob_dates <- dat %>%
      dplyr::select(c(NewLocationID, Date, OpState)) %>%
      distinct() %>%
      mutate(Date = as.Date(Date, format = "%d-%b-%Y")) %>%
      #'  Label sequential dates from same camera as a single burst
      group_by(NewLocationID) %>%
        mutate(Burst = cumsum(c(1, diff(Date) > 1))) %>%
      ungroup() %>%
      #'  filter to just the start and end date of each burst
      group_by(NewLocationID, Burst) %>%
        mutate(StartEnd = ifelse(Date == min(Date), "First", NA),
               StartEnd = ifelse(Date == max(Date), "Last", StartEnd)) %>%
        filter(!is.na(StartEnd)) %>%
      # filter(Date == min(Date) | Date == max(Date)) %>%
      ungroup()
    #'  Filter down to rows where only one day at a camera site had issues - leads
    #'  to only "Last" date, no "First" date
    lone_date <- prob_dates %>%
      filter((NewLocationID != lag(NewLocationID)) & (StartEnd == "Last" & lag(StartEnd == "Last"))) %>%
      mutate(StartEnd = "First")
    #'  Filter down to cameras where there were multiple date ranges with problems, 
    #'  including single days with problems that need a "First" date duplicated
    extra_single_day <- prob_dates %>%
      filter((NewLocationID == lag(NewLocationID)) & (StartEnd == "Last" & lag(StartEnd == "Last"))) %>%
      mutate(StartEnd = "First")
    #'  Merge with larger data set so cameras with only one problem day have a
    #'  start and end date for that problem day
    prob_dates <- rbind(prob_dates, lone_date, extra_single_day) %>%
      arrange(NewLocationID, Date, StartEnd)
    return(prob_dates)
  }
  eoe_prob_dates_20s <- prob_days(eoe_1hr_20s) %>%
    #'  Drop this one observation - same camera labeled w/ 2 problems on same day
    filter(NewLocationID != "GMU10A_U_67" | OpState != "completely obscured")
  eoe_prob_dates_20w <- prob_days(eoe_1hr_20w) %>%
    #'  Drop these observation - same camera labeled w/ 2 problems on same day
    filter(NewLocationID != "GMU10A_U_23" | OpState != "completely obscured") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "completely obscured" | Date != "2020-11-01") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "nightbad__dayok" | Date != "2021-01-02") %>%
    filter(NewLocationID != "GMU6_P_34" | OpState != "nightbad__dayok" | Date != "2020-12-15") 
  eoe_prob_dates_21s <- prob_days(eoe_1hr_21s)
  
  #'  Merge motion and trigger datasets together for wolf cameras
  wolf_prob_dates_19s <- lapply(wolf_1hr_19s, prob_days) %>%
    do.call(rbind.data.frame, .) %>%
    #'  Remove unnamed camera with no coordinate information
    filter(NewLocationID != "NA_O_NA")
  wolf_prob_dates_20s <- lapply(wolf_1hr_20s, prob_days)
  wolf_prob_dates_21s <- lapply(wolf_1hr_21s, prob_days) %>%
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, Date, StartEnd)
  #'  Need to add a "First" row for this data set since lag in function above doesn't
  #'  work if first observation is labeled "Last"
  first_row <- wolf_prob_dates_21s[1,] %>%
    mutate(StartEnd = ifelse(StartEnd == "Last", "First", StartEnd))
  wolf_prob_dates_21s <- bind_rows(first_row, wolf_prob_dates_21s)

  
  #' wolf_prob_dates_19s <- prob_days(wolf_1hr_19s) %>%
  #'   #'  Remove unnamed camera with no coordinate information
  #'   filter(NewLocationID != "NA_O_NA")
  #' wolf_prob_dates_20s <- prob_days(wolf_1hr_20s)
  #' wolf_prob_dates_21s <- prob_days(wolf_1hr_21s)
  
  #'  Save dates of first/last image and problem dates for all cameras
  write.csv(eoe_prob_dates_20s, "./Data/IDFG camera data/Problem cams/eoe20s_OpState_probs.csv")
  write.csv(eoe_prob_dates_20w, "./Data/IDFG camera data/Problem cams/eoe20w_OpState_probs.csv")
  write.csv(eoe_prob_dates_21s, "./Data/IDFG camera data/Problem cams/eoe21s_OpState_probs.csv")
  
  write.csv(wolf_prob_dates_19s, "./Data/IDFG camera data/Problem cams/wolf19s_OpState_probs.csv")
  write.csv(wolf_prob_dates_20s, "./Data/IDFG camera data/Problem cams/wolf20s_OpState_probs.csv")
  write.csv(wolf_prob_dates_21s, "./Data/IDFG camera data/Problem cams/wolf21s_OpState_probs.csv")
  
  #'  Create wide data frame based on each burst per camera
  #'  Start and end dates will be listed horizontally by camera burst
  problem_df <- function(dat) {
    wide_format <- dat %>% 
      mutate(Problem = ifelse(StartEnd == "First", paste0("Problem", Burst, "_from"), paste0("Problem", Burst, "_to"))) %>% 
             # Problem = paste0(Burst, "_", StartEnd)) %>%
      dplyr::select(-c(OpState, Burst, StartEnd)) %>%
      spread(Problem, Date)    
    return(wide_format)
  }
  eoe_wide_probs_20s <- problem_df(eoe_prob_dates_20s) %>%
    relocate(c(`Problem10_from`, `Problem10_to`, `Problem11_from`, `Problem11_to`, `Problem12_from`, `Problem12_to`, `Problem13_from`, `Problem13_to`, `Problem14_from`, `Problem14_to`), .after = `Problem9_to`)
  eoe_wide_probs_20w <- problem_df(eoe_prob_dates_20w) %>%
    relocate(c(`Problem1_from`, `Problem1_to`, `Problem2_from`, `Problem2_to`, `Problem3_from`, `Problem3_to`, `Problem4_from`, `Problem4_to`, `Problem5_from`, `Problem5_to`, `Problem6_from`, `Problem6_to`, `Problem7_from`, `Problem7_to`, `Problem8_from`, `Problem8_to`, `Problem9_from`, `Problem9_to`), .after = `NewLocationID`)
  eoe_wide_probs_21s <- problem_df(eoe_prob_dates_21s) %>%
    relocate(c(`Problem10_from`, `Problem10_to`, `Problem11_from`, `Problem11_to`, `Problem12_from`, `Problem12_to`), .after = `Problem9_to`)
  
  wolf_wide_probs_19s <- problem_df(wolf_prob_dates_19s) %>%
    relocate(c(`Problem10_from`, `Problem10_to`, `Problem11_from`, `Problem11_to`, `Problem12_from`, `Problem12_to`), .after = `Problem9_to`)
  wolf_wide_probs_20s <- problem_df(wolf_prob_dates_20s) %>%
    relocate(c(`Problem10_from`, `Problem10_to`, `Problem11_from`, `Problem11_to`, `Problem12_from`, `Problem12_to`, `Problem13_from`, `Problem13_to`), .after = `Problem9_to`)
  wolf_wide_probs_21s <- problem_df(wolf_prob_dates_21s)
  
  
  #'  Extract dates of first and last images taken by camera per season
  start_stop_dates <- function(dat) {
    startend <- dat %>%
      #'  Group by each camera site
      group_by(NewLocationID) %>%
      #'  Arrange by date and time
      arrange(posix_date_time) %>%
      #'  Filter to just the first and last observation per group
      slice(c(1, n())) %>%
      ungroup() %>%
      dplyr::select(c(NewLocationID, Date, Lat, Long)) %>%
      #'  Remove camera if missing NewLocationID - currently only one EoE Smr21 
      #'  camera with this problem - need to find out where these data came from!
      filter(!is.na(NewLocationID)) %>%
      #'  Add a column representing whether the date is the start or end of photos
      mutate(Date = as.Date(Date, format = "%d-%b-%Y"),
             StartEnd = ifelse((NewLocationID == lag(NewLocationID)), "Retrieval_date", "Setup_date"), 
             StartEnd = ifelse(is.na(StartEnd), "Setup_date", StartEnd)) %>%
      spread(StartEnd, Date) %>%
      relocate(Retrieval_date, .after = Setup_date)
    return(startend)
  }
  #'  Merge deployment data with problem dates
  eoe_probcams_20s <- start_stop_dates(eoe20s_allT) %>%
    full_join(eoe_wide_probs_20s, by = "NewLocationID")
  eoe_probcams_20w <- start_stop_dates(eoe20w_allT) %>%
    full_join(eoe_wide_probs_20w, by = "NewLocationID")
  eoe_probcams_21s <- start_stop_dates(eoe21s_allT) %>%
    full_join(eoe_wide_probs_21s, by = "NewLocationID")
    
  wolf_probcams_19s <- start_stop_dates(wolf19s_allT) %>%
    full_join(wolf_wide_probs_19s, by = "NewLocationID") %>%
    #'  Remove unnamed camera with no coordinate information
    filter(NewLocationID != "NA_O_NA")
  wolf_probcams_20s <- start_stop_dates(wolf20s_allT) %>%
    full_join(wolf_wide_probs_20s, by = "NewLocationID")
  wolf_probcams_21s <- start_stop_dates(wolf21s_allT) %>%
    full_join(wolf_wide_probs_21s, by = "NewLocationID")
    
  #'  Save dates of first/last image and problem dates for all cameras
  save(eoe_probcams_20s, file = "./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  write.csv(eoe_probcams_20s, "./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.csv")
  save(eoe_probcams_20w, file = "./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  write.csv(eoe_probcams_20w, "./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.csv")
  save(eoe_probcams_21s, file = "./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  write.csv(eoe_probcams_21s, "./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.csv")
  
  save(wolf_probcams_19s, file = "./Data/IDFG camera data/Problem cams/wolf19s_problem_cams.RData")
  write.csv(wolf_probcams_19s, "./Data/IDFG camera data/Problem cams/wolf19s_problem_cams.csv")
  save(wolf_probcams_20s, file = "./Data/IDFG camera data/Problem cams/wolf20s_problem_cams.RData")
  write.csv(wolf_probcams_20s, "./Data/IDFG camera data/Problem cams/wolf20s_problem_cams.csv")
  save(wolf_probcams_21s, file = "./Data/IDFG camera data/Problem cams/wolf21s_problem_cams.RData")
  write.csv(wolf_probcams_21s, "./Data/IDFG camera data/Problem cams/wolf21s_problem_cams.csv")
  
  
    #'  From here, add full start and end dates to this df,
  #'  review partial misdirection data to see if most of is based on animal disturbance
  #'  double check sever misdirection is an all season long issue
  #'  look for date ranges where cameras just didn't take any pictures and add these problem date ranges to the above df
  
  
  #'  Pull out images from any day were the sequential_probs function indicated
  #'  potential problem - includes normal images to help assess the extent of problem
  bad_day_cams <- function(dat, prob_condition) {
    bad_days <- semi_join(dat, prob_condition, by = c("Date", "NewLocationID")) %>%
      # dplyr::select(c(NewLocationID, Date, Time, posix_date_time, OpState, Species, File)) %>%
      arrange(NewLocationID, posix_date_time)
    return(bad_days)
  }
  #'  Save all images from camera on problem days - to be used to record start 
  #'  and stop problem days and indicated which images should be filtered out
  eoe_badday_t20s <- bad_day_cams(eoe_t_20s_probs, prob_condition = eoe_1hr_20s)
  save(eoe_badday_t20s, file = "./Data/IDFG camera data/Problem cams/eoe20s_t_bad_cam_days.RData")
  write.csv(eoe_badday_t20s, "./Data/IDFG camera data/Problem cams/eoe20s_t_bad_cam_days.csv")
  eoe_badday_t20w <- bad_day_cams(eoe_t_20w_probs, prob_condition = eoe_1hr_20w)
  save(eoe_badday_t20w, file = "./Data/IDFG camera data/Problem cams/eoe20w_t_bad_cam_days.RData")
  write.csv(eoe_badday_t20w, "./Data/IDFG camera data/Problem cams/eoe20w_t_bad_cam_days.csv")
  eoe_badday_t21s <- bad_day_cams(eoe_t_21s_probs, prob_condition = eoe_1hr_21s)
  save(eoe_badday_t21s, file = "./Data/IDFG camera data/Problem cams/eoe21s_t_bad_cam_days.RData")
  write.csv(eoe_badday_t21s, "./Data/IDFG camera data/Problem cams/eoe21s_t_bad_cam_days.csv")
  
  wolf_badday_t19s <- bad_day_cams(wolf_t_19s_probs, prob_condition = wolf_1hr_19s)
  save(wolf_badday_t19s, file = "./Data/IDFG camera data/Problem cams/wolf19s_t_bad_cam_days.RData")
  write.csv(wolf_badday_t19s, "./Data/IDFG camera data/Problem cams/wolf19s_t_bad_cam_days.csv")
  wolf_badday_t20s <- bad_day_cams(wolf_t_20s_probs, prob_condition = wolf_1hr_20s)
  save(wolf_badday_t20s, file = "./Data/IDFG camera data/Problem cams/wolf20s_t_bad_cam_days.RData")
  write.csv(wolf_badday_t20s, "./Data/IDFG camera data/Problem cams/wolf20s_t_bad_cam_days.csv")
  wolf_badday_t21s <- bad_day_cams(wolf_t_21s_probs, prob_condition = wolf_1hr_21s)
  save(wolf_badday_t21s, file = "./Data/IDFG camera data/Problem cams/wolf21s_t_bad_cam_days.RData")
  write.csv(wolf_badday_t21s, "./Data/IDFG camera data/Problem cams/wolf21s_t_bad_cam_days.csv")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  