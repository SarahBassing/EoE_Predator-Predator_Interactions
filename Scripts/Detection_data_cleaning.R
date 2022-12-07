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
  #'  Set time zone on keeper data sets
  load("./Data/IDFG camera data/Split datasets/eoe_motion_skinny_NewLocationID.RData")
  eoe_motion_list <- lapply(eoe_motion_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/eoe_time_skinny_NewLocationID.RData")
  eoe_noon_list <- lapply(eoe_time_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/wolf_motion_skinny_NewLocationID.RData")
  wolf_motion_list <- lapply(wolf_motion_list, set_tzone)
  load("./Data/IDFG camera data/Split datasets/wolf_time_skinny_NewLocationID.RData")
  wolf_noon_list <- lapply(wolf_time_list, set_tzone)
  
  #'  Set time zone on full data sets
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM_NewLocationID.RData")
  eoe20s_allM <- set_tzone(eoe20s_allM)
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM_NewLocationID.RData")
  eoe20w_allM <- set_tzone(eoe20w_allM)
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM_NewLocationID.RData")
  eoe21s_allM <- set_tzone(eoe21s_allM) %>%
    mutate(Setup = "U or P whatever") %>%
    relocate(Setup, .after = NewLocationID)
  
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
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allT_NewLocationID.RData")
  wolf19s_allT <- set_tzone(wolf19s_allT)
  load("./Data/IDFG camera data/Split datasets/wolf20s_allT_NewLocationID.RData")
  wolf20s_allT <- set_tzone(wolf20s_allT)
  load("./Data/IDFG camera data/Split datasets/wolf21s_allT_NewLocationID.RData")
  wolf21s_allT <- set_tzone(wolf21s_allT)
  
  
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
  
  
  ####  Flag gaps in camera data  ####
  #'  ----------------------------
  #'  Calculate length of gaps in data due to NO pictures being taken
  no_pix <- function(dat) {
    gaps <- dat %>%
      dplyr::select(NewLocationID, Lat, Long, Date, Time, posix_date_time, OpState) %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Floor each time to nearest minute (drops seconds to 00)
      mutate(Floordt = floor_date(posix_date_time, "minutes")) %>%
      group_by(NewLocationID) %>%
      #'  Calculate time since last image at each camera site
      mutate(nmin = as.numeric(difftime(Floordt, lag(Floordt), units = "mins"))) %>%
      ungroup() %>%
      #'  Report in hours and days
      mutate(nhrs = round(nmin/60, 2),
             ndays = round(nhrs/24, 2))
    return(gaps)
  }
  eoe_nopix_20s <- no_pix(eoe20s_allT)
  eoe_nopix_20w <- no_pix(eoe20w_allT)
  eoe_nopix_21s <- no_pix(eoe21s_allT)
  
  
  ########## FROM HERE integrate into functions below? Get start/end date based on 1, 6, 12 hr rule and incorporate into problem cams
  
 
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
  eoe20s_all <- list(eoe20s_allM, eoe20s_allT)
  eoe20s_probs <- lapply(eoe20s_all, problem_children)
  eoe20w_all <- list(eoe20w_allM, eoe20w_allT)
  eoe20w_probs <- lapply(eoe20w_all, problem_children)
  eoe21s_all <- list(eoe21s_allM, eoe21s_allT)
  eoe21s_probs <- lapply(eoe21s_all, problem_children)
  
  wolf19s_all <- list(wolf19s_allM, wolf19s_allT)
  wolf19s_probs <- lapply(wolf19s_all, problem_children)
  wolf20s_all <- list(wolf20s_allM, wolf20s_allT)
  wolf20s_probs <- lapply(wolf20s_all, problem_children)
  wolf21s_all <- list(wolf21s_allM, wolf21s_allT)
  wolf21s_probs <- lapply(wolf21s_all, problem_children)
  

  #'  Filter images to series where camera was obscured or misdirected for a defined amount of time
  #'  Using this to represent days when camera wasn't full operational
  sequential_probs <- function(dat, ntime) {
    #'  Create condition for data to meet: camera's view is completely obscured, 
    #'  severely misdirected, or the images are nightbad_dayOK
    cond <- expr(dat$OpState == "completely obscured" | dat$OpState == "severely misdirected" | 
                   dat$OpState == "nightbad__dayok")
    #'  Look through sequential images and if they meet this condition for the
    #'  defined time (ntime) then flag and filter data set to just these images
    bad_view_pix <- dat %>% 
      mutate(problem_pix = 
               rep(rle(!!cond)$values & rle(!!cond)$lengths >= ntime, 
                   rle(!!cond)$lengths)) %>%
      filter(problem_pix) %>%
      arrange(NewLocationID, posix_date_time)
    
    return(bad_view_pix)
  }
  #'  Flag problematic images from each full time-trigger data set
  #'  ntime: 6 = 1 hr, 36 = 6 hrs, 72 = 12 hrs
  eoe_1hr_20s <- lapply(eoe20s_probs, sequential_probs, ntime = 6) %>% 
    #'  Merge motion & time trigger data sets into one dataframe - important for next step
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  eoe_1hr_20w <- lapply(eoe20w_probs, sequential_probs, ntime = 6) %>%
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  eoe_1hr_21s <- lapply(eoe21s_probs, sequential_probs, ntime = 6) %>%
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  
  wolf_1hr_19s <- lapply(wolf19s_probs, sequential_probs, ntime = 6) %>%
    #'  Merge motion & time trigger data sets into one dataframe
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  wolf_1hr_20s <- lapply(wolf20s_probs, sequential_probs, ntime = 6) %>%
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  wolf_1hr_21s <- lapply(wolf21s_probs, sequential_probs, ntime = 6) %>%
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID, posix_date_time)
  
  
  #'  Plot distribution of problem times - is there a temporal pattern?
  eoe_probtimes_20s <- eoe_1hr_20s %>%  #  Summer 2020
    mutate(Time = hms::as_hms(Time),
           Hours = hour(Time)) %>%
    ggplot(aes(Time)) + #Hours
    geom_freqpoly(binwidth = 1) +
    ggtitle("Number of images w/ OpState problems over 24 cycle \n(1hr Rule)")
  eoe_probtimes_20s
  ggsave("./Outputs/Figures/eoe_probtimes_20s_1hr.png", eoe_probtimes_20s, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  eoe_probtimes_20w <- eoe_1hr_20w %>%  #  Winter 2020-2021
    mutate(Time = hms::as_hms(Time),
           Hours = hour(Time)) %>%
    ggplot(aes(Time)) + #Hours
    geom_freqpoly(binwidth = 1) +
    ggtitle("Number of images w/ OpState problems over 24 cycle \n(1hr Rule)")
  eoe_probtimes_20w
  ggsave("./Outputs/Figures/eoe_probtimes_20w_1hr.png", eoe_probtimes_20w, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  
  eoe_probtimes_21s <- eoe_1hr_21s %>%  #  Summer 2021
    mutate(Time = hms::as_hms(Time),
           Hours = hour(Time)) %>%
    ggplot(aes(Time)) + #Hours
    geom_freqpoly(binwidth = 1) +
    ggtitle("Number of images w/ OpState problems over 24 cycle \n(1hr Rule)")
  eoe_probtimes_21s
  ggsave("./Outputs/Figures/eoe_probtimes_21s_1hr.png", eoe_probtimes_21s, units = "in", 
         width = 6, height = 6, dpi = 600, device = "png")
  

  #'  Pull out problem dates based on images when camera is obscured/misdirected
  #'  for 1+ hour on a given day
  prob_days <- function(dat) {
    prob_dates <- dat %>%
      arrange(NewLocationID, posix_date_time) %>%
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
    #'  Calculate length of each problem date range
    prob_dates <- prob_dates %>%
      group_by(NewLocationID, Burst) %>%
      mutate(ndays = as.numeric(difftime(Date, lag(Date), units = "days")),
             #'  Add 1 to each count so 1 day problems = 1, not 0, and start day
             #'  of each problem is included in the count
             ndays = ndays + 1) %>%   
      ungroup()
    return(prob_dates)
  }
  eoe_prob_dates_20s <- prob_days(eoe_1hr_20s) %>% 
    arrange(NewLocationID, Date, StartEnd) %>%
    #'  Drop this one observation - same camera labeled w/ 2 problems on same day 
    filter(NewLocationID != "GMU10A_U_67" | OpState != "completely obscured") %>%
    filter(NewLocationID != "GMU6_U_123" | OpState != "completely obscured") %>%
    filter(NewLocationID != "GMU6_U_99" | OpState != "completely obscured") %>% # STOP HERE if ntime > 6
    filter(NewLocationID != "GMU10A_P_110" | OpState != "completely obscured" | Date != "2020-10-11") %>%
    filter(NewLocationID != "GMU10A_P_110" | OpState != "completely obscured" | Date != "2020-10-14") %>%
    filter(NewLocationID != "GMU10A_P_110" | OpState != "completely obscured" | Date != "2020-10-16") %>%
    filter(NewLocationID != "GMU10A_P_110" | OpState != "severely misdirected" | Date != "2020-10-16" | StartEnd != "First") %>%
    mutate(ndays = ifelse(NewLocationID == "GMU10A_P_110" & OpState == "severely misdirected" & Date == "2020-10-16" & StartEnd == "Last", 3, ndays)) %>%
    filter(NewLocationID != "GMU10A_U_37" | OpState != "severely misdirected")
  eoe_prob_dates_20w <- prob_days(eoe_1hr_20w) %>%
    arrange(NewLocationID, Date, StartEnd) %>%
    #'  Drop these observation - same camera labeled w/ 2 problems on same day
    filter(NewLocationID != "GMU10A_U_23" | OpState != "completely obscured") %>%
    filter(NewLocationID != "GMU6_P_34" | OpState != "nightbad__dayok" | Date != "2020-12-15") %>%
    filter(NewLocationID != "GMU6_P_92" | OpState != "nightbad__dayok") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "completely obscured" | Date != "2020-11-01") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "completely obscured" | Date != "2020-11-05") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "completely obscured" | Date != "2021-01-02") %>%
    filter(NewLocationID != "GMU6_P_67" | OpState != "nightbad__dayok" | Date != "2021-01-02" | StartEnd != "Last") %>%
    mutate(StartEnd = ifelse(NewLocationID == "GMU6_P_67" & OpState == "nightbad__dayok" & Date == "2021-01-02", "Last", StartEnd)) %>%
    filter(NewLocationID != "GMU6_U_7" | OpState != "completely obscured") #%>%  # STOP HERE if ntime = 6, DON'T DO if ntime = 72
    # filter(NewLocationID != "GMU6_P_67" | OpState != "completely obscured" | Date != "2020-11-02") # Needed for ntime = 36 only
    # filter(NewLocationID != "GMU6_P_34" | Date != "2021-03-22" | StartEnd != "Last") %>% # Needed for ntime = 72
    # mutate(StartEnd = ifelse(NewLocationID == "GMU6_P_34" & Date == "2021-03-22", "Last", StartEnd)) %>% # Needed for ntime = 72
    # filter(NewLocationID != "GMU10A_U_78" | Date != "2021-02-06" | StartEnd != "Last") %>% # Needed for ntime = 72
    # mutate(StartEnd = ifelse(NewLocationID == "GMU10A_U_78" & Date == "2021-02-06", "Last", StartEnd)) # Needed for ntime = 72
  eoe_prob_dates_21s <- prob_days(eoe_1hr_21s) %>%
    arrange(NewLocationID, Date, StartEnd)
  
  wolf_prob_dates_19s <- prob_days(wolf_1hr_19s) %>%
    arrange(NewLocationID, Date, StartEnd) %>%
    #'  Drop this observation - same camera labeled w/ 2 problems on same day
    filter(NewLocationID != "GMU44_A_3110" | OpState != "completely obscured")
  wolf_prob_dates_20s <- prob_days(wolf_1hr_20s) %>%
    arrange(NewLocationID, Date, StartEnd)
  wolf_prob_dates_21s <- prob_days(wolf_1hr_21s) %>%
    arrange(NewLocationID, Date, StartEnd)
  
  #'  Add a "First" row to start of dataframe - necessary for wolf data sets
  #'  where first camera had only 1 problem day and lag in function above can't 
  #'  account for it so first row of dataframe is the "Last" problem date
  add_row1 <- function(dat) {
    first_row <- dat[1,] %>%
      mutate(StartEnd = ifelse(StartEnd == "Last", "First", StartEnd))
    dat <- bind_rows(first_row, dat) 
    return(dat)
  }
  wolf_prob_dates_20s <- add_row1(wolf_prob_dates_20s) %>%
    #'  Drop these observation - same camera labeled w/ 2 problems on same day
    filter(NewLocationID != "GMU10A_A_750" | OpState != "completely obscured")
  wolf_prob_dates_21s <- add_row1(wolf_prob_dates_21s)

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
      dplyr::select(-c(OpState, Burst, StartEnd, ndays)) %>%
      #'  Spread dataframe so problem dates are in a wide format (necessary if 
      #'  using camtrapR to create detection histories)
      spread(Problem, Date) %>%
      #'  Make sure problem columns are sequential
      relocate(c(`Problem1_from`, `Problem1_to`, `Problem2_from`, `Problem2_to`, `Problem3_from`, `Problem3_to`, 
                 `Problem4_from`, `Problem4_to`, `Problem5_from`, `Problem5_to`, `Problem6_from`, `Problem6_to`, 
                 `Problem7_from`, `Problem7_to`, `Problem8_from`, `Problem8_to`, `Problem9_from`, `Problem9_to`), .after = `NewLocationID`)
    
    return(wide_format)
  }
  eoe_wide_probs_20s <- problem_df(eoe_prob_dates_20s)   
  eoe_wide_probs_20w <- problem_df(eoe_prob_dates_20w)  
  eoe_wide_probs_21s <- problem_df(eoe_prob_dates_21s) 
  
  wolf_wide_probs_19s <- problem_df(wolf_prob_dates_19s) 
  wolf_wide_probs_20s <- problem_df(wolf_prob_dates_20s) 
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
  #'  Using setup/retrieval dates from only time-trigger data set for EoE cameras
  #'  because all cameras were set up to take both motion & time trigger images.
  #'  Time trigger data set is generally a better representation of start/stop 
  #'  and problem dates in the data.
  eoe_probcams_20s <- start_stop_dates(eoe20s_allT) %>%
    arrange(NewLocationID) %>% 
    full_join(eoe_wide_probs_20s, by = "NewLocationID")
  eoe_probcams_20w <- start_stop_dates(eoe20w_allT) %>%
    arrange(NewLocationID) %>% 
    full_join(eoe_wide_probs_20w, by = "NewLocationID")
  eoe_probcams_21s <- start_stop_dates(eoe21s_allT) %>%
    arrange(NewLocationID) %>% 
    full_join(eoe_wide_probs_21s, by = "NewLocationID")
  
  #'  Using both motion & time trigger data sets for wolf cameras b/c not all 
  #'  cameras were taking both types of image  
  wolf_probcams_19s <- lapply(wolf19s_all, start_stop_dates) %>% 
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID) %>%
    full_join(wolf_wide_probs_19s, by = "NewLocationID") %>%
    #'  Remove unnamed camera with no coordinate information
    filter(NewLocationID != "NA_O_NA")
  wolf_probcams_20s <- lapply(wolf20s_all, start_stop_dates) %>% 
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID) %>%
    full_join(wolf_wide_probs_20s, by = "NewLocationID")
  wolf_probcams_21s <- lapply(wolf21s_all, start_stop_dates) %>% 
    do.call(rbind.data.frame, .) %>%
    arrange(NewLocationID) %>%
    full_join(wolf_wide_probs_21s, by = "NewLocationID") %>%
    #'  Remove unnamed camera with no coordinate information
    filter(NewLocationID != "GMUNA_NA")
    
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
  
  
  ####  Summary stats on problem cameras  ####
  #'  -------------------------------------
  #'  Calculate summary stats on duration of typical problem days
  mean_problem_days <- function(dat) {
    #'  Overall mean, range, boxplot, median
    print("Mean number of problem days")
    print(mean(dat$ndays, na.rm = TRUE))
    print("Range of problem days")
    print(range(dat$ndays, na.rm = TRUE))
    boxplot(dat$ndays, na.rm = TRUE)
    
    print("Median number of problem days")
    print(median(dat$ndays, na.rm = TRUE))
    
    #'  A few cameras have MANY problems- weighted average to account for this
    mean_prob_days_wgtd <- dat %>%
      group_by(NewLocationID) %>%
      summarize(wgt = n()/2,
                mean_days = mean(ndays, na.rm = TRUE)) %>%
      ungroup()
    sum_wtd <- sum(mean_prob_days_wgtd$wgt)
    mean_prob_days_wgtd$wgt <- mean_prob_days_wgtd$wgt/sum_wtd
    mean_wgtd_prob_days <- weighted.mean(mean_prob_days_wgtd$mean_days, mean_prob_days_wgtd$wgt)
    print("Weighted mean number of problem days, weighted by # of problems per camera")
    print(mean_wgtd_prob_days)
    
    #'  Mean, range, boxplot, median of first problem per cameras (most cameras only had one)
    print("Mean number of problem days for first problem only")
    print(mean(dat$ndays[dat$Burst == 1], na.rm = TRUE))
    print("Range of problem days for first problem only")
    print(range(dat$ndays[dat$Burst == 1], na.rm = TRUE))
    boxplot(dat$ndays[dat$Burst == 1], na.rm = TRUE)
    print("Median number of problem days for first problem only")
    print(median(dat$ndays[dat$Burst == 1], na.rm = TRUE))
  }
  mean_problem_days(eoe_prob_dates_20s)
  mean_problem_days(eoe_prob_dates_20w)
  mean_problem_days(eoe_prob_dates_21s)
  
  mean_problem_days(wolf_prob_dates_19s)
  mean_problem_days(wolf_prob_dates_20s)
  mean_problem_days(wolf_prob_dates_21s)
  
  #'  What proportion of cameras had problems
  proportion_of_probcams <- function(dat) {
    #'  Total number of cameras
    ncams <- nrow(dat)
    #'  Number of cameras with 1 or more problems
    ncams_1prob <- length(na.omit(dat$Problem1_from))
    #'  Number of cameras with 2 or more problems
    ncams_2prob <- length(na.omit(dat$Problem2_from))
    #'  Number of cameras with 3 or more problems
    ncams_3prob <- length(na.omit(dat$Problem3_from))
    #'  Number of cameras with 10 or more problems
    ncams_10prob <- length(na.omit(dat$Problem10_from))
    nmbr_prob_cams <- c(ncams_1prob, ncams_2prob, ncams_3prob, ncams_10prob)
    proportion_prob_cams <- round(nmbr_prob_cams/ncams, 3)
    nmbr_probs <- c("1 or more", "2 or more", "3 or more", "4 or more")
    prop_prob_cams <- as.data.frame(cbind(nmbr_probs, nmbr_prob_cams, proportion_prob_cams)) 
    new_name <- paste0("% problem cams (n = ", ncams, ")")
    names(prop_prob_cams)[names(prop_prob_cams) == 'proportion_prob_cams'] <- new_name
    return(prop_prob_cams)
  }
  eoe20s_prop_probcams <- proportion_of_probcams(eoe_probcams_20s)
  eoe20w_prop_probcams <- proportion_of_probcams(eoe_probcams_20w)
  eoe21s_prop_probcams <- proportion_of_probcams(eoe_probcams_21s)
  
  wolf19s_prop_probcams <- proportion_of_probcams(wolf_probcams_19s)
  wolf20s_prop_probcams <- proportion_of_probcams(wolf_probcams_20s)
  wolf21s_prop_probcams <- proportion_of_probcams(wolf_probcams_21s)
  
  
  ####  Map problem cameras  ####
  #'  -----------------------
  #'  Read in shapefiles
  gmu <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
  eoe_gmus <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp")
  
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection(eoe_gmus)
  
  #'  Function to format camera location data for mapping
  format_probcam_locations <- function(dat) {
    #'  Label cameras as either having or not having a problem
    flag_problem_cams <- dat %>%
      mutate(ProblemCams = ifelse(is.na(Problem1_from), "No problems", "Problem Camera")) %>%
      dplyr::select(c(NewLocationID, Lat, Long, ProblemCams))
    
    #'  Make camera locations spatial
    cams <- st_as_sf(flag_problem_cams, coords = c("Long", "Lat"), crs = wgs84) %>%
      st_transform(cams, crs = sa_proj)
    return(cams)
  }
  eoe20s_probcam_locs <- format_probcam_locations(eoe_probcams_20s)
  eoe20w_probcam_locs <- format_probcam_locations(eoe_probcams_20w)
  eoe21s_probcam_locs <- format_probcam_locations(eoe_probcams_21s)
  
  wolf19s_probcam_locs <- format_probcam_locations(wolf_probcams_19s)
  wolf20s_probcam_locs <- format_probcam_locations(wolf_probcams_20s)
  wolf21s_probcam_locs <- format_probcam_locations(wolf_probcams_21s)
  

  #'  Plot Smr20 problem cameras
  EoE20s_problem_cam_map <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe20s_probcam_locs, aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "red")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ Operational Problems (1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20s_problem_cam_map
  # ggsave("./Outputs/Figures/Maps/EoE20s_problem_cams_1hr.png", EoE20s_problem_cam_map, units = "in", 
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  #'  Plot Wtr20 problem cameras
  EoE20w_problem_cam_map <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe20w_probcam_locs, aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "#edae49")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ Operational Problems (1h Rule)") +# UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20w_problem_cam_map
  # ggsave("./Outputs/Figures/Maps/EoE20w_problem_cams_1hr.png", EoE20w_problem_cam_map, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  #'  Plot Smr21 problem cameras
  EoE21s_problem_cam_map <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus, aes(fill = NAME)) +
    scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
    geom_sf(data = eoe21s_probcam_locs, aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "#darkred")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ Operational Problems (1h Rule)") +
    coord_sf(xlim = c(-13020000, -12700000), ylim = c(5800000, 6274865), expand = TRUE) +
    theme_bw()
  EoE21s_problem_cam_map
  ggsave("./Outputs/Figures/Maps/EoE21s_problem_cams_1hr.png", EoE21s_problem_cam_map, units = "in",
         width = 6, height = 6, dpi = 600, device = "png")
  
  
  #'  Map problem cameras where problem lasts longer than 1 week and 1 month
  wk_month_long_probs <- function(obs, cams) {
    #'  Snag camera coordinates from full observation data
    camlocs <- obs %>%
      dplyr::select(c(NewLocationID, Lat, Long)) %>%
      #'  Group by each camera site
      group_by(NewLocationID) %>%
      #'  Filter to just the first observation per camera
      slice(1) %>%
      ungroup() %>%
      filter(!is.na(NewLocationID))
    
    #'  Identify and format cameras that had a problem for >1 week long
    weeklong_probs <- cams %>%
      #'  Filter to only problems that lasted longer than 1 week
      filter(ndays >7) %>%
      #'  Retain only the first problem at each camera (doesn't matter if there
      #'  was more than one problem, just that at least one 1-wk problem occurred)
      group_by(NewLocationID) %>%
      slice(1) %>%
      ungroup() %>%
      #'  Merge with camera coordinate data
      full_join(camlocs, by = "NewLocationID") %>%
      #'  Flag each problem as having a problem or not
      mutate(ProblemCams = ifelse(is.na(ndays), "No problems", "Problem Camera")) %>%
      dplyr::select(c(NewLocationID, Lat, Long, ProblemCams)) %>%
      #'  Convert to a sf object for mapping
      st_as_sf(coords = c("Long", "Lat"), crs = wgs84) %>%
      st_transform(crs = sa_proj)
    
    #'  Identify and format cameras that had a problem for >1 month (30 days) long
    #'  Same deal as above but filtering to cameras with >30 day long problems
    monthlong_probs <- cams %>%
      filter(ndays >30) %>%
      full_join(camlocs, by = "NewLocationID") %>%
      mutate(ProblemCams = ifelse(is.na(ndays), "No problems", "Problem Camera")) %>%
      dplyr::select(c(NewLocationID, Lat, Long, ProblemCams)) %>%
      st_as_sf(coords = c("Long", "Lat"), crs = wgs84) %>%
      st_transform(crs = sa_proj)
    
    #'  List week and month-long problem cameras together and return this object
    prob_list <- list(weeklong_probs, monthlong_probs)
    return(prob_list)
  }
  eoe_prob_wkmo_20s <- wk_month_long_probs(obs = eoe20s_allT, cams = eoe_prob_dates_20s)
  eoe_prob_wkmo_20w <- wk_month_long_probs(obs = eoe20w_allT, cams = eoe_prob_dates_20w)
  eoe_prob_wkmo_21s <- wk_month_long_probs(obs = eoe21s_allT, cams = eoe_prob_dates_21s)
  
 
  #'  Plot cameras with >1week long problems and >1month long problems
  #'  Plot Smr20 problem cameras
  EoE20s_problem_cam_map_1wk <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_20s[[1]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "red")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Week Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20s_problem_cam_map_1wk
  # ggsave("./Outputs/Figures/Maps/EoE20s_problem_cams_1hr_1wk.png", EoE20s_problem_cam_map_1wk, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  EoE20s_problem_cam_map_1mo <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_20s[[2]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "red")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Month Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20s_problem_cam_map_1mo
  # ggsave("./Outputs/Figures/Maps/EoE20s_problem_cams_1hr_1mo.png", EoE20s_problem_cam_map_1mo, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  #'  Plot Wtr20 problem cameras
  EoE20w_problem_cam_map_1wk <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_20w[[1]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "#edae49")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Week Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20w_problem_cam_map_1wk
  # ggsave("./Outputs/Figures/Maps/EoE20w_problem_cams_1hr_1wk.png", EoE20w_problem_cam_map_1wk, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  EoE20w_problem_cam_map_1mo <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME != "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_20w[[2]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "#edae49")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Month Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12830000), ylim = c(5800000, 6050000), expand = TRUE) +
    theme_bw()
  EoE20w_problem_cam_map_1mo
  # ggsave("./Outputs/Figures/Maps/EoE20w_problem_cams_1hr_1mo.png", EoE20w_problem_cam_map_1mo, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  #'  Plot Smr21 problem cameras
  EoE21s_problem_cam_map_1wk <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus, aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_21s[[1]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "darkred")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Week Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12700000), ylim = c(5800000, 6274865), expand = TRUE) +
    theme_bw()
  EoE21s_problem_cam_map_1wk
  # ggsave("./Outputs/Figures/Maps/EoE21s_problem_cams_1hr_1wk.png", EoE21s_problem_cam_map_1wk, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  EoE21s_problem_cam_map_1mo <- ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus, aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC", "#66CC99")) +
    geom_sf(data = eoe_prob_wkmo_21s[[2]], aes(color = ProblemCams), shape = 16, size = 2) +
    scale_color_manual(values=c("black", "darkred")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("EoE Cameras w/ >1 Month Long Operational Problems \n(1hr Rule)") + # UPDATE THE TIME HERE
    coord_sf(xlim = c(-13020000, -12700000), ylim = c(5800000, 6274865), expand = TRUE) +
    theme_bw()
  EoE21s_problem_cam_map_1mo
  # ggsave("./Outputs/Figures/Maps/EoE21s_problem_cams_1hr_1mo.png", EoE21s_problem_cam_map_1mo, units = "in",
  #        width = 6, height = 6, dpi = 600, device = "png")
  
  
  
  #'  From here, look for date ranges where cameras just didn't take any pictures and 
  #'  add these problem date ranges to the above df
  #'  Don't forget to look at the Wolf data!
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  