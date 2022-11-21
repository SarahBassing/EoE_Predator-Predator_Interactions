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
  
  eoe_noon_list <- mapply(eoe_deploy_info, season = eoe_seasons, dets = eoe_time_skinny, pred = "predator", SIMPLIFY = FALSE)  
  #'  Fix NewLocationID info- this was recorded differently for Smr21 images in 
  #'  original data set so unnecessary info gets added in this function - need to
  #'  revert back to original information
  eoe_noon_list[[3]] <- dplyr::select(eoe_noon_list[[3]], -c(NewLocationID)) %>%
    mutate(NewLocationID = LocationID) %>%
    relocate(NewLocationID, .after = LocationID)
  
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
  
  eoe_motion_smr20 <- set_tzone(eoe_motion_smr20)
  eoe_motion_wtr20 <- set_tzone(eoe_motion_wtr20)
  eoe_motion_smr21 <- set_tzone(eoe_motion_smr21)
  
  eoe_time_smr20 <- set_tzone(eoe_time_smr20)
  eoe_time_wtr20 <- set_tzone(eoe_time_wtr20)
  eoe_time_smr21 <- set_tzone(eoe_time_smr21)
  
  wolf_motion_smr20 <- set_tzone(wolf_motion_smr20)
  wolf_motion_wtr20 <- set_tzone(wolf_motion_wtr20)
  wolf_motion_smr21 <- set_tzone(wolf_motion_smr21)
  
  wolf_time_smr20 <- set_tzone(wolf_time_smr20)
  wolf_time_wtr20 <- set_tzone(wolf_time_wtr20)
  wolf_time_smr21 <- set_tzone(wolf_time_smr21)

  
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
  #'  Pull out any camera with a noted misdirected/obscured viewshed
  #'  Looking for images labeled: "severely misdirected", "partially obscured", 
  #'  "completely obscured", "malfunction"
  #'  Images marked "minorly misdirected" are still OK to use per S.Thompson
  #'  "partially obscured": 25-75% obscured but probably can still see some animals
  problem_children <- function(dat) {
    #'  Retain NewLocationID info of cameras with known problems
    prob_cams <- as.data.frame(unique(dat$NewLocationID[dat$OpState == "severely misdirected" |
                                                          dat$OpState == "partially obscured" |
                                                          dat$OpState == "completely obscured" |
                                                          dat$OpState == "malfunction"]))
    colnames(prob_cams) <- "NewLocationID"
    #'  Retain all images from just problem cameras
    prob_pix <- semi_join(dat, prob_cams, by = "NewLocationID")
    #'  Double check I have the same number of cameras
    print(nrow(prob_cams))
    print(length(unique(prob_pix$NewLocationID)))
    return(prob_pix)
  }
  eoe_probs <- lapply(eoe_motion_list, problem_children)
  wolf_probs <- lapply(wolf_motion_list, problem_children)
  
  
  
  
  
  
  