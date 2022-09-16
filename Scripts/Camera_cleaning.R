  #'  ---------------------------------
  #'  Cleaning raw camera data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean deployment and deteciton data, make sure annual datasets 
  #'  are consistent, filter and merge.
  #'  ---------------------------------
  
  #'  Load libraries
  library(stringr)  
  library(tidyverse)
  
  #'  Camera deployment data
  load("./Data/IDFG camera data/eoe2020s_cameras.RData") 
  load("./Data/IDFG camera data/eoe2020w_cameras.RData") 
  load("./Data/IDFG camera data/eoe2021s_cameras.RData")
  
  load("./Data/IDFG camera data/wolf2019s_cameras.RData")
  load("./Data/IDFG camera data/wolf2020s_cameras.RData")
  load("./Data/IDFG camera data/wolf2021s_cameras.RData")
  
  ####  Format date/time data  ####
  #'  -------------------------
  #'  Identify time zones of deployment data (start/end dates)
  #'  Detetion data is in MST (America/Edmonton (UTC-07:00); tz="America/Edmonton")
  attr(cams_s20_eoe$Start_Date, "tzone") #"UTC"
  attr(cams_s19_wolf$Start_Date, "tzone") #NULL
  
  #'  Format data when time is included in Date_Deployed and Date_Retrieved
  format_time1 <- function(dat) {
    whatistime <- dat %>%
      mutate(
        #'  Format to POSIXct based on attributes above
        Date_Deployed = as.POSIXct(Date_Deployed, format="%m/%d/%Y %H:%M", tz="UTC"),
        Date_Retrieved = as.POSIXct(Date_Retrieved, format="%m/%d/%Y %H:%M", tz="UTC"),
        Date_Dep = as.Date(Date_Deployed, format = "%Y-%m-%d"),
        Date_Pull = as.Date(Date_Retrieved, format = "%Y-%m-%d")) %>%
      dplyr::select(-c(Date_Deployed, Date_Retrieved)) %>%
      relocate(Date_Dep, .before = Observers) %>%
      rename(Date_Deployed = Date_Dep) %>%
      relocate(Date_Pull, .before = Observers) %>%
      rename(Date_Retrieved = Date_Pull)
    return(whatistime)
  }
  cams_s20_eoe <- format_time1(cams_s20_eoe)
  cams_s21_eoe <- format_time1(cams_s21_eoe)
  
  #'  Format data when time is included in Date_Deployed and Date_Retrieved
  cams_w20_eoe <- mutate(cams_w20_eoe,
        Date_Dep = as.POSIXct(Date_Deployed, format="%m/%d/%Y", tz="UTC"),
        Date_Pull = as.POSIXct(Date_Retrieved, format="%m/%d/%Y", tz="UTC")) %>%
      dplyr::select(-c(Date_Deployed, Date_Retrieved)) %>%
      relocate(Date_Dep, .before = Observers) %>%
      rename(Date_Deployed = Date_Dep) %>%
      relocate(Date_Pull, .before = Observers) %>%
      rename(Date_Retrieved = Date_Pull)
  
  cams_s19_wolf <- mutate(cams_s19_wolf,
        Date_Deployed = as.Date(Date_Deployed, format = "%Y-%m-%d"),
        Date_Retrieved = as.Date(Date_Retrieved, format = "%Y-%m-%d"),
        Start_Date = as.Date(Start_Date, format = "%Y-%m-%d"),
        End_Date = as.Date(End_Date, format = "%Y-%m-%d"))
  
  cams_s20_wolf <- mutate(cams_s20_wolf,
        Date_Deployed = as.Date(Date_Deployed, format = "%m/%d/%Y"),
        Date_Retrieved = as.Date(Date_Retrieved, format = "%m/%d/%Y"),
        Start_dt = as.POSIXct(Start_Date, format = "%Y-%m-%d %H:%M:%S", tz="UTC"),
        End_dt = as.POSIXct(End_Date, format = "%Y-%m-%d %H:%M:%S", tz="UTC"),
        Start_dt = as.Date(Start_dt, format = "%Y-%m-%d"),
        End_dt = as.Date(End_dt, format = "%Y-%m-%d"))  %>%
      dplyr::select(-c(Start_Date, End_Date)) %>%
      relocate(Start_dt, .before = Observers) %>%
      rename(Start_Date = Start_dt) %>%
      relocate(End_dt, .before = Observers) %>%
      rename(End_Date = End_dt)
  
  cams_s21_wolf <- mutate(cams_s21_wolf,
        Date_Deployed = as.Date(Date_Deployed, format = "%m/%d/%Y"),
        Date_Retrieved = as.Date(Date_Retrieved, format = "%m/%d/%Y"),
        Start_dt = as.POSIXct(Start_Date, format = "%m/%d/%Y %H:%M", tz="UTC"),
        End_dt = as.POSIXct(End_Date, format = "%m/%d/%Y %H:%M", tz="UTC"),
        Start_dt = as.Date(Start_dt, format = "%Y-%m-%d"),
        End_dt = as.Date(End_dt, format = "%Y-%m-%d"))  %>%
      dplyr::select(-c(Start_Date, End_Date)) %>%
      relocate(Start_dt, .before = Observers) %>%
      rename(Start_Date = Start_dt) %>%
      relocate(End_dt, .before = Observers) %>%
      rename(End_Date = End_dt)
  
  
  cams_s21_wolf <- format_time2(cams_s21_wolf)
  
  ####  Consistent LocationID naming structure  ####
  #'  ------------------------------------------
  #'  Step 1: make sure GMUs with "A"s are always uppercase
  uppercase_time <- function(dat) {
    uppercase_A <- dat %>%
      mutate(Gmu = toupper(Gmu), 
             LocationID = toupper(LocationID))
    return(uppercase_A)
  }
  cams_s20_eoe <- uppercase_time(cams_s20_eoe) %>% mutate(Season = "Smr20")
  cams_w20_eoe <- uppercase_time(cams_w20_eoe) %>% mutate(Season = "Wtr20")
  cams_s21_eoe <- uppercase_time(cams_s21_eoe) %>% mutate(Season = "Smr21")
  cams_s19_wolf <- uppercase_time(cams_s19_wolf) %>% mutate(Season = "Smr19")
  cams_s20_wolf <- uppercase_time(cams_s20_wolf) %>% mutate(Season = "Smr20")
  cams_s21_wolf <- uppercase_time(cams_s21_wolf) %>% mutate(Season = "Smr21")
  #'  Step 2: differentiate Abundance and Occupancy Cameras for wolf data sets
  AOC <- function(dat) {
    cam_target <- dat %>%
      mutate(NewTarget = ifelse(Target == "Abundance", "A", "O"),
             NewTarget = ifelse(Target == "Abund_Occu", "B", NewTarget),
             NewLocationID = paste0(NewTarget, "_", LocationID)) %>%
      dplyr::select(-c(NewTarget, LocationID)) %>%
      rename(LocationID = NewLocationID) %>%
      relocate(LocationID, .after = "TargetSpecies")
    return(cam_target)
  }
  cams_s19_wolf <- AOC(cams_s19_wolf)
  cams_s20_wolf <- AOC(cams_s20_wolf)
  #'  Step 3: combine GMU, Setup, and Location information into a single LocationID
  new_LocationID <- function(dat) {
    camera_station <- dat %>%
      mutate(NewSetup = ifelse(Setup == "ungulate", "U", "P"),
             # NewSetup = ifelse(Setup == "WLF", "WLF", NewSetup), 
             NewLocationID = paste0("GMU", Gmu, "_", NewSetup, "_", LocationID)) %>%
      dplyr::select(-c(NewSetup, LocationID)) %>%
      rename(LocationID = NewLocationID) %>%
      relocate(LocationID, .after = "TargetSpecies")
    return(camera_station)
  }
  cams_s20_eoe <- new_LocationID(cams_s20_eoe)
  cams_w20_eoe <- new_LocationID(cams_w20_eoe)
  #'  Slight tweaks for wolf cams
  new_LocationID_wolf1 <- function(dat) {
    camera_station <- dat %>%
      mutate(NewLocationID = paste0(Gmu, "_", Setup, LocationID)) %>%
      dplyr::select(-LocationID) %>%
      rename(LocationID = NewLocationID) %>%
      relocate(LocationID, .after = "TargetSpecies")
    return(camera_station)
  }
  cams_s19_wolf <- new_LocationID_wolf1(cams_s19_wolf)
  new_LocationID_wolf2 <- function(dat) {
    camera_station <- dat %>%
      mutate(NewLocationID = paste0("GMU", Gmu, "_", Setup, LocationID)) %>%
      dplyr::select(-LocationID) %>%
      rename(LocationID = NewLocationID) %>%
      relocate(LocationID, .after = "TargetSpecies")
    return(camera_station)
  }
  cams_s20_wolf <- new_LocationID_wolf2(cams_s20_wolf)
  cams_s21_wolf <- new_LocationID_wolf2(cams_s21_wolf)
  
  #'  Merge all deployment data together
  eoe_cams <- rbind(cams_s20_eoe, cams_w20_eoe, cams_s21_eoe)
  length(unique(eoe_cams$LocationID))
  wolf_cams <- rbind(cams_s19_wolf, cams_s20_wolf, cams_s21_wolf)
  length(unique(wolf_cams$LocationID))
  
  ####  Explore deployment locations  ####
  #'  --------------------------------
  #'  Identify which locations had a camera at it for >1 season
  #'  CamID doesn't really matter so much as whether a camera was deployed at the
  #'  same LocationID with matching Lat/Long coordinates for >1 season 
  samesame <- function(dat1, dat2, season1, season2) {
    #'  Which cameras from season 1 stayed in the same place in season 2?
    same_locID1 <- dat1$LocationID %in% dat2$LocationID
    same_Lat1 <- dat1$Lat %in% dat2$Lat
    same_Long1 <- dat1$Long %in% dat2$Long
    same_camID1 <- dat1$CamID %in% dat2$CamID
    #'  Which cameras from season 2 were at the same location as the previous season?
    same_locID2 <- dat2$LocationID %in% dat1$LocationID
    same_Lat2 <- dat2$Lat %in% dat1$Lat
    same_Long2 <- dat2$Long %in% dat1$Long
    same_camID2 <- dat2$CamID %in% dat1$CamID
    
    #'  Add indicators to each data set to track consistencies across seasons
    #'  For sites monitored in season 1, is the same LocationID, etc. monitored in season 2?
    cams_t1to2 <- dat1
    cams_t1to2$same_locID <- same_locID1
    cams_t1to2$same_Lat <- same_Lat1
    cams_t1to2$same_Long <- same_Long1
    cams_t1to2$same_camID <- same_camID1
    cams_t1to2 <- dplyr::select(cams_t1to2, c(Region, Gmu, Setup, Target, LocationID, Lat, Long, CamID, same_locID, same_Lat, same_Long, same_camID)) %>%
      mutate(Monitor_season = season1) %>%
      relocate(Monitor_season, .after = CamID)
    #'  For sites montored in season 2, was the same LocationID, etc. monitored the previous season?
    cams_t2from1 <- dat2
    cams_t2from1$same_locID <- same_locID2
    cams_t2from1$same_Lat <- same_Lat2
    cams_t2from1$same_Long <- same_Long2
    cams_t2from1$same_camID <- same_camID2
    cams_t2from1 <- dplyr::select(cams_t2from1, c(Region, Gmu, Setup, Target, LocationID, Lat, Long, CamID, same_locID, same_Lat, same_Long, same_camID)) %>%
      mutate(Monitor_season = season2) %>%
      relocate(Monitor_season, .after = CamID)
    
    #'  Which LocationIDs, Lat, Long, & CamIDs appear in >1 season
    multi_season_cam_deploy <- rbind(cams_t1to2, cams_t2from1) %>%
      arrange(LocationID)
    return(multi_season_cam_deploy)
  }
  same_s20_w2021_eoe <- samesame(dat1 = cams_s20_eoe, dat2 = cams_w20_eoe, season1 = "Smr20", season2 = "Wtr2021")
  same_s20_s21_eoe <- samesame(dat1 = cams_s20_eoe, dat2 = cams_s21_eoe, season1 = "Smr20", season2 = "Smr21")
  same_s19_s20_wolf <- samesame(dat1 = cams_s19_wolf, dat2 = cams_s20_wolf, season1 = "Smr19", season2 = "Smr20")
  same_s19_s21_wolf <- samesame(dat1 = cams_s19_wolf, dat2 = cams_s21_wolf, season1 = "Smr19", season2 = "Smr21")
  same_s20_s21_wolf <- samesame(dat1 = cams_s20_wolf, dat2 = cams_s21_wolf, season1 = "Smr20", season2 = "Smr21")
  
  ####  Potential problem cameras  ####
  #'  -----------------------------
  #'  Cameras with same deployment & retrieval date
  same_startend_eoe <- filter(eoe_cams, Date_Deployed == Date_Retrieved)
  same_startend_wolf <- filter(wolf_cams, Date_Deployed == Date_Retrieved) #none
  
  #'  Cameras with deploy or retrieval dates outside focal monitoring season
  longdeploy_eoe 
  #'  Cameras with missing GMU information
  #'  Wolf locations that were only monitored for 1 summer
  #'  Look at covariates that weren't collected
  #'  What are some of these other fields?
  
  
  
  
  