  #'  ---------------------------------
  #'  Cleaning raw detection data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean camea deteciton data, make sure annual datasets are consistent, 
  #'  filter into manageable data sets, and connect to camera locations.
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
  
  #' #'  Full detection data sets
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
  # load("./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  # load("./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  # load("./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  #' 
  #' load("./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  #' load("./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  #'  Keeper data sets (all animal/human/vehicle images & noon timelapse images)
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
  eoe_seasons <- list("Smr20", "Wtr20", "Smr21")
  eoe_motion_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_motion_skinny, pred = "predator", SIMPLIFY = FALSE) 
  eoe_noon_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  #'  Double check it worked
  eoe21s_noon <- eoe_noon_list[[3]]
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  eoe_motion_list[[3]] <- dplyr::select(eoe_motion_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  eoe_noon_list[[3]] <- dplyr::select(eoe_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  
  eoe_motion_wtr20 <- eoe_deploy_info(dets = eoe20w_allM, season = "Wtr20", pred = "predator")
  
  #' #'  Run full eoe timelapse data sets through function
  #' eoe20s_allT <- eoe_deploy_info(season = "Smr20", dets = eoe20s_allT, pred = "predator")
  #' eoe20w_allT <- eoe_deploy_info(season = "Wtr20", dets = eoe20w_allT, pred = "predator")
  #' eoe21s_allT <- eoe_deploy_info(season = "Smr21", dets = eoe21s_allT, pred = "predator")
  #' save(eoe20s_allT, file = "./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  #' save(eoe20w_allT, file = "./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  #' save(eoe21s_allT, file = "./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  
  
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
             NewLocationID = paste0("GMU", Gmu, "_", Setup, "_", LocationID)) %>%
      dplyr::select(-c(acams, aocams)) %>%
      relocate(NewLocationID, .before = File) %>%
      relocate(Setup, .after = NewLocationID)
    return(dets)
  }
  wolf_seasons <- list("Smr19", "Smr20", "Smr21")
  wolf_motion_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_motion_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)  
  wolf_noon_list <- mapply(wolf_deploy_info, season = wolf_seasons, dets = wolf_time_skinny, abund = "Abundance", abund_occu = "Abund_Occu", SIMPLIFY = FALSE)  
  #'  Double check it worked
  wolf21s_noon <- wolf_noon_list[[3]]
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  wolf_motion_list[[3]] <- dplyr::select(wolf_motion_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  wolf_noon_list[[3]] <- dplyr::select(wolf_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  
  #' #'  Run full wolf timelapse data sets through function
  #' wolf19s_allT <- wolf_deploy_info(season = "Smr19", dets = wolf19s_allT, pred = "predator")
  #' wolf20s_allT <- wolf_deploy_info(season = "Smr20", dets = wolf20s_allT, pred = "predator")
  #' wolf21s_allT <- wolf_deploy_info(season = "Smr21", dets = wolf21s_allT, pred = "predator")
  #' save(wolf19s_allT, file = "./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  #' save(wolf20s_allT, file = "./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  #' save(wolf21s_allT, file = "./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  
  ####  Format & visualize date/time data  ####
  #  ---------------------------------------
  #'  Detection data was recorded in MST (America/Edmonton (UTC-07:00); tz="America/Edmonton")
  set_tzone <- function(dat) {
    tz(dat$posix_date_time) <- "America/Edmonton"
    print(attr(dat$posix_date_time, "tzone"))
    print(range(dat$posix_date_time, na.rm = T))
    return(dat)
  }
  eoe_motion_list <- lapply(eoe_motion_list, set_tzone)
  eoe_noon_list <- lapply(eoe_noon_list, set_tzone)
  wolf_motion_list <- lapply(wolf_motion_list, set_tzone)
  wolf_noon_list <- lapply(wolf_noon_list, set_tzone)
  
  eoe_motion_wtr20 <- set_tzone(eoe_motion_wtr20)

  #'  Visualize observations over time
  #'  https://r4ds.had.co.nz/dates-and-times.html
  #'  Note that when you use date-times in a numeric context (like in a histogram), 
  #'  1 means 1 second, so a binwidth of 86400 means one day. For dates, 1 means 1 day.
  #'  Noon timelapse images only- represents camera function, not animal activity/detections
  eoe_noon_list[[1]] %>%  #  Summer 2020
    filter(posix_date_time < ymd(20210101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) + # 86400 seconds = 1 day
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active EoE cameras, Summer 2020")
  
  eoe_noon_list[[2]] %>%  #  Winter 2020-2021
    filter(posix_date_time > ymd(20201101) & posix_date_time < ymd(20210531)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active EoE cameras, Winter 2020-2021")
  
  eoe_noon_list[[3]] %>%  #  Summer 2021
    filter(posix_date_time > ymd(20210501) & posix_date_time < ymd(20211101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active EoE cameras, Summer 2021")
  
  wolf_noon_list[[1]] %>%  # Summer 2019
    filter(posix_date_time < ymd(20200101)) %>% 
    filter(Time == "12:0:0") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) + # 86400 seconds = 1 day
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active wolf cameras, Summer 2019")
  
  wolf_noon_list[[2]] %>%  #  Summer 2020
    filter(posix_date_time < ymd(20210101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active wolf cameras, Summer 2020")
  
  wolf_noon_list[[3]] %>%  #  Summer 2021
    filter(posix_date_time < ymd(20220101)) %>% 
    filter(Time == "12:00:00") %>%
    ggplot(aes(posix_date_time)) + 
    geom_freqpoly(binwidth = 86400) +
    scale_x_datetime(
      name = "Dates",
      breaks = scales::date_breaks("month")) +
    ggtitle("Number of active wolf cameras, Summer 2021")
  
  
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
  
  eoe_mismatch <- mt_mismatch(motion_cams = eoe_motion_wtr20, noon_cams = eoe_noon_list[[2]])
  
  #'  Review motion triggered images
  eoe_m_smr20 <- eoe_motion_list[[1]]
  prob_eoe_m_smr20 <- eoe_m_smr20[eoe_m_smr20$NewLocationID == "GMU6_U_11" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_81" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_P_11" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_36" | 
                                    eoe_m_smr20$NewLocationID == "GMU6_U_129" | 
                                    eoe_m_smr20$NewLocationID == "GMU10A_U_160" | 
                                    eoe_m_smr20$NewLocationID == "GMU10A_U_8",]
  eoe_m_wtr20 <- eoe_motion_list[[2]]
  prob_eoe_m_smr20 <- eoe_m_wtr20[eoe_m_wtr20$NewLocationID == "GMU6_U_164",]
  eoe_m_smr21 <- eoe_motion_list[[3]]
  
  wolf_m_smr20 <- wolf_motion_list[[1]]
  wolf_m_wtr20 <- wolf_motion_list[[2]]
  wolf_m_smr21 <- wolf_motion_list[[3]]
  
  #'  Review noontime triggered images
  eoe_t_smr20 <- eoe_noon_list[[1]]
  prob_eoe_t_smr20 <- eoe_t_smr20[eoe_t_smr20$NewLocationID == "GMU6_U_69" | 
                                    eoe_t_smr20$NewLocationID == "GMU6_U_105" | 
                                    eoe_t_smr20$NewLocationID == "GMU10A_U_11",]
  eoe_t_wtr20 <- eoe_noon_list[[2]]
  eoe_t_smr21 <- eoe_noon_list[[3]]
  
  wolf_t_smr20 <- wolf_noon_list[[1]]
  wolf_t_wtr20 <- wolf_noon_list[[2]]
  wolf_t_smr21 <- wolf_noon_list[[3]]
  
  
  
  
  
  
  
  
  
  