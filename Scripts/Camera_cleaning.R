  #'  ---------------------------------
  #'  Cleaning raw camera data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean deployment and detection data, make sure annual data sets 
  #'  are consistent, filter and merge.
  #'  ---------------------------------
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(stringr)
  library(sf)
  library(raster)
  library(ggplot2)
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
      mutate(Target = ifelse(Target == "AbundOccu", "Abund_Occu", Target),
             NewTarget = ifelse(Target == "Abundance", "A", "O"),
             NewTarget = ifelse(Target == "Abund_Occu", "B", NewTarget),
             NewLocationID = paste0(NewTarget, "_", LocationID)) %>%
      dplyr::select(-c(NewTarget)) %>% #, LocationID
      # rename(LocationID = NewLocationID) %>%
      relocate(LocationID, .after = "TargetSpecies") %>%
      relocate(NewLocationID, .after = "LocationID")
    return(cam_target)
  }
  cams_s19_wolf <- AOC(cams_s19_wolf)
  cams_s20_wolf <- AOC(cams_s20_wolf)
  
  #'  Step 3: combine GMU, Setup, and Location information into a single LocationID
  new_LocationID <- function(dat) {
    camera_station <- dat %>%
      mutate(NewSetup = ifelse(Setup == "ungulate", "U", "P"),
             NewLocationID = paste0("GMU", Gmu, "_", NewSetup, "_", LocationID)) %>%
      dplyr::select(-c(NewSetup)) %>% 
      relocate(LocationID, .after = "TargetSpecies") %>%
      relocate(NewLocationID, .after = "LocationID")
    return(camera_station)
  }
  cams_s20_eoe <- new_LocationID(cams_s20_eoe)
  cams_w20_eoe <- new_LocationID(cams_w20_eoe)
  #'  Slight tweak for summer 2021 data - need to generate LocationID that matches
  #'  previous seasons
  cams_s21_eoe <- mutate(cams_s21_eoe, NewLocationID = LocationID, 
                         #'  Remove everything before underscore (do twice b/c 2 underscores)
                         LocationID = sub(".*_", "", LocationID),
                         LocationID = sub(".*_", "", LocationID)) %>%
    relocate(NewLocationID, .after = "LocationID")
  #'  Slight tweaks for wolf cams
  new_LocationID_wolf1 <- function(dat) {
    camera_station <- dat %>%
      mutate(NewLocationID = paste0(Gmu, "_", NewLocationID)) %>% #Setup, 
      relocate(LocationID, .after = "TargetSpecies") %>%
      relocate(NewLocationID, .after = "LocationID")
    return(camera_station)
  }
  cams_s19_wolf <- new_LocationID_wolf1(cams_s19_wolf)
  new_LocationID_wolf2 <- function(dat) {
    camera_station <- dat %>%
      mutate(NewLocationID = paste0("GMU", Gmu, "_", NewLocationID)) %>% #Setup, 
      relocate(LocationID, .after = "TargetSpecies") %>%
      relocate(NewLocationID, .after = "LocationID")
    return(camera_station)
  }
  cams_s20_wolf <- new_LocationID_wolf2(cams_s20_wolf)
  #'  summer 2021 wolf LocationID already has A/O info so just rename as NewLocationID
  cams_s21_wolf <- mutate(cams_s21_wolf, NewLocationID = LocationID)
  cams_s21_wolf <- new_LocationID_wolf2(cams_s21_wolf)
  
  #'  Merge all deployment data together
  eoe_cams <- rbind(cams_s20_eoe, cams_w20_eoe, cams_s21_eoe)
  length(unique(eoe_cams$LocationID))
  summary(eoe_cams)
  wolf_cams <- rbind(cams_s19_wolf, cams_s20_wolf, cams_s21_wolf)
  length(unique(wolf_cams$LocationID))
  summary(wolf_cams)
  
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
    cams_t1to2 <- dplyr::select(cams_t1to2, c(Region, Gmu, Setup, Target, LocationID, Lat, Long, CamID, CameraHeight_M, same_locID, same_Lat, same_Long, same_camID)) %>%
      mutate(Monitor_season = season1) %>%
      relocate(Monitor_season, .after = CamID)
    #'  For sites monitored in season 2, was the same LocationID, etc. monitored the previous season?
    cams_t2from1 <- dat2
    cams_t2from1$same_locID <- same_locID2
    cams_t2from1$same_Lat <- same_Lat2
    cams_t2from1$same_Long <- same_Long2
    cams_t2from1$same_camID <- same_camID2
    cams_t2from1 <- dplyr::select(cams_t2from1, c(Region, Gmu, Setup, Target, LocationID, Lat, Long, CamID, CameraHeight_M, same_locID, same_Lat, same_Long, same_camID)) %>%
      mutate(Monitor_season = season2) %>%
      relocate(Monitor_season, .after = CamID)
    
    #'  Which LocationIDs, Lat, Long, & CamIDs appear in >1 season
    multi_season_cam_deploy <- rbind(cams_t1to2, cams_t2from1) %>%
      arrange(LocationID)
    return(multi_season_cam_deploy)
  }
  same_s20_w2021_eoe <- samesame(dat1 = cams_s20_eoe, dat2 = cams_w20_eoe, season1 = "Smr20", season2 = "Wtr2021")
  same_w2021_s20_eoe <- samesame(dat1 = cams_w20_eoe, dat2 = cams_s20_eoe, season1 = "Wtr2021", season2 = "Smr20")
  same_s20_s21_eoe <- samesame(dat1 = cams_s20_eoe, dat2 = cams_s21_eoe, season1 = "Smr20", season2 = "Smr21")
  same_s19_s20_wolf <- samesame(dat1 = cams_s19_wolf, dat2 = cams_s20_wolf, season1 = "Smr19", season2 = "Smr20")
  same_s19_s21_wolf <- samesame(dat1 = cams_s19_wolf, dat2 = cams_s21_wolf, season1 = "Smr19", season2 = "Smr21")
  same_s20_s21_wolf <- samesame(dat1 = cams_s20_wolf, dat2 = cams_s21_wolf, season1 = "Smr20", season2 = "Smr21")
  
  #'  How many cameras did not run the next season or were new in the current season?
  #'  Summer EoE cameras that did not run in winter
  missing_w2021_eoe <- filter(same_s20_w2021_eoe, same_locID == FALSE) #202 smr20 cams did not run wtr2021
  #'  Winter EoE cameras that did not run in previous summer
  new_w2021_eoe <- filter(same_w2021_s20_eoe, Monitor_season == "Wtr2021" & same_locID == FALSE) #no new cameras went out for wtr2021
  #'  Summer EoE cameras that did not run the next summer
  missing_smr_eoe <- filter(same_s20_s21_eoe, same_locID == FALSE) # 254 smr21 cams but all in GMU1 where cams first went out smr21
  
  #'  Wolf locations that were only monitored for 1 summer
  #'  No consistent GMU or style of camera (Abundance vs Occupancy)
  deployed20 <- cams_s19_wolf$LocationID %in% cams_s20_wolf$LocationID
  deployed21 <- cams_s19_wolf$LocationID %in% cams_s21_wolf$LocationID
  lone_cams_s19_wolf <- cbind(cams_s19_wolf, deployed20, deployed21) %>%
    filter(deployed20 == FALSE) %>%
    filter(deployed21 == FALSE)                                           # 141 smr19 cams did not run again
  deployed19 <- cams_s20_wolf$LocationID %in% cams_s19_wolf$LocationID
  deployed21 <- cams_s20_wolf$LocationID %in% cams_s21_wolf$LocationID
  lone_cams_s20_wolf <- cbind(cams_s20_wolf, deployed19, deployed21) %>%
    filter(deployed19 == FALSE) %>%
    filter(deployed21 == FALSE)                                           # 44 smr20 cams did not run previously or again
  deployed19 <- cams_s21_wolf$LocationID %in% cams_s19_wolf$LocationID
  deployed20 <- cams_s21_wolf$LocationID %in% cams_s20_wolf$LocationID
  lone_cams_s21_wolf <- cbind(cams_s21_wolf, deployed19, deployed20) %>%
    filter(deployed19 == FALSE) %>%
    filter(deployed20 == FALSE)                                           # 61 smr20 cams did not run previously
  
  ####  Reorganize data for me  ####
  #'  ---------------------------
  #'  Function to reorganize camera data so it's easier to use
  organize_cols <- function(dat) {
    dat <- dat %>%
      #'  Drop unneeded columns
      dplyr::select(-c(Project, TargetSpecies, LocationIDSource, CameraBearing, Observers, Strata, GeneralLocation, M_pics, T_pics)) %>%
      #'  Add "NA" if cell left blank  
      mutate(across(c("CameraFacing", "MarkerType_A", "MarkerType_B", "DominantHabitatType", "Topography"), 
                    ~ifelse(.=="", NA, as.character(.)))) %>%
      #'  Fill in field for ungulate cameras, which monitored truly random locations
      mutate(CameraFacing = ifelse(is.na(CameraFacing) & Setup == "ungulate", "random", CameraFacing)) %>%
      #'  Move columns around so they're better organized
      relocate(Season, .after = CamID) %>%
      relocate(Date_Deployed, .after = Season) %>%
      relocate(Date_Retrieved, .after = Date_Deployed) %>%
      relocate(Start_Date, .after = Date_Retrieved) %>%
      relocate(End_Date, .after = Start_Date) %>%
      relocate(CameraHeight_M, .after = End_Date) %>%
      relocate(CameraFacing, .after = CameraHeight_M) %>%
      relocate(CamStatus_Arrival, .after = last_col()) %>%
      relocate(CamStatus_Leaving, .after = last_col())
    return(dat)
  }
  #'  Reformat data in long format
  eoe_cams_long <- arrange(eoe_cams, LocationID) %>%
    organize_cols(.)
  # write.csv(eoe_cams_long, file = "./Data/IDFG camera data/cams_eoe_long.csv")
  wolf_cams_long <- arrange(wolf_cams, LocationID) %>%
    organize_cols(.)
  # write.csv(wolf_cams_long, file = "./Data/IDFG camera data/cams_wolf_long.csv")
  
  #'  Reformat data in wide format
  #'  Resulting data frame is bananas b/c so many repeat columns that I don't 
  #'  know if they're important or not
  eoe_list <- list(cams_s20_eoe, cams_w20_eoe, cams_s21_eoe)
  eoe_trim <- lapply(eoe_list, organize_cols)
  eoe_cams_wide <- full_join(eoe_trim[[1]], eoe_trim[[2]], by = c("Region", "Gmu", "Setup", "Target", "LocationID", "AreaType")) %>%
    full_join(eoe_trim[[3]], by = c("Region", "Gmu", "Setup", "Target", "LocationID", "AreaType")) %>%
    dplyr::select(-c(DominantHabitatType.x, DominantHabitatType.y, Topography.x, Topography.y, CanopyCover.x, CanopyCover.y))
  #'  Make at least some of these repeat columns easier to keep track of
  #'  Final season's worth of data retains original column names
  names(eoe_cams_wide) <- gsub(".x", "_smr20", names(eoe_cams_wide), fixed = T)
  names(eoe_cams_wide) <- gsub(".y", "_wtr2021", names(eoe_cams_wide), fixed = T)
  
  wolf_list <- list(cams_s19_wolf, cams_s20_wolf, cams_s21_wolf)
  wolf_trim <- lapply(wolf_list, organize_cols)
  wolf_cams_wide <- full_join(wolf_trim[[1]], wolf_trim[[2]], by = c("Region", "Gmu", "Setup", "Target", "LocationID", "AreaType", "MarkerDistance_A_M", "MarkerType_A", "MarkerDistance_B_M", "MarkerType_B", "MarkerDistance_C_M", "MarkerType_C")) %>%
    full_join(wolf_trim[[3]], by = c("Region", "Gmu", "Setup", "Target", "LocationID", "AreaType", "MarkerDistance_A_M", "MarkerType_A", "MarkerDistance_B_M", "MarkerType_B", "MarkerDistance_C_M", "MarkerType_C")) %>%
    dplyr::select(-c(DominantHabitatType.x, DominantHabitatType.y, Topography.x, Topography.y, CanopyCover.x, CanopyCover.y))

  ####  Potential problem cameras  ####
  #'  -----------------------------
  #'  Cameras with same deployment & retrieval date
  same_startend_eoe <- filter(eoe_cams, Date_Deployed == Date_Retrieved | Start_Date == End_Date)
  same_startend_wolf <- filter(wolf_cams, Date_Deployed == Date_Retrieved | Start_Date == End_Date) 
  
  #'  Cameras with odd locations relative to GMUs
  prob_locs <- eoe_cams_long[eoe_cams_long$NewLocationID == "GMU6_P_109" | eoe_cams_long$NewLocationID == "GMU6_U_122" | eoe_cams_long$NewLocationID == "GMU6_U_109" | eoe_cams_long$NewLocationID == "GMU6_P_63" | eoe_cams_long$NewLocationID == "GMU6_U_63" | eoe_cams_long$NewLocationID == "GMU10A_U_101" | eoe_cams_long$NewLocationID == "UNKNOWN",]
  #write.csv(prob_locs, "./Data/IDFG camera data/problem_eoe_cam_locs.csv")
  
  
  #### NEED TO FINISH CLEANING WOLF CAMERAS TO IDENTIFY POTENTIAL PROBLEMS & MAKE SKINNY VERSION  ####
  
  
  #'  Reduced camera location data set to something more usable
  #'  For simplicity, retain each unique camera location with first lat/long coordinates
  #'  associated with the location and some site-specific sampling effort/habitat data 
  eoe_cams_skinny <- dplyr::select(eoe_cams_wide, c("Region", "Gmu", "Setup", "Target", "LocationID", "Lat_smr20", "Long_smr20", "CameraHeight_M_smr20", "CameraFacing_smr20", "DominantHabitatType", "Topography", "CanopyCover")) %>%
    #'  Drop data from GMU1 b/c these weren't deployed in 2020
    filter(Gmu != "1")
  names(eoe_cams_skinny)[names(eoe_cams_skinny) == "Lat_smr20"] <- "Lat"
  names(eoe_cams_skinny)[names(eoe_cams_skinny) == "Long_smr20"] <- "Long"
  names(eoe_cams_skinny)[names(eoe_cams_skinny) == "CameraHeight_M_smr20"] <- "CameraHeight_M"
  names(eoe_cams_skinny)[names(eoe_cams_skinny) == "CameraFacing_smr20"] <- "CameraFacing"
  eoe_cams_s21_skinny <- dplyr::select(eoe_cams_wide, c("Region", "Gmu", "Setup", "Target", "LocationID", "Lat", "Long", "CameraHeight_M", "CameraFacing", "DominantHabitatType", "Topography", "CanopyCover")) %>%
    #'  Retain only data from GMU1 to match data from cameras deployed in GMU6 & GMU10A in 2020
    filter(Gmu == "1")
  eoe_cams_skinny <- rbind(eoe_cams_skinny, eoe_cams_s21_skinny) %>%
    arrange(LocationID)
  # write.csv(eoe_cams_skinny, file = "./Data/IDFG camera data/cams_eoe_skinny.csv")
  
  
  ####  Visualize these locations  ####
  #'  ------------------------------
  #'  Read in shapefiles
  gmu <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
  eoe_gmus <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp")
  
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection(eoe_gmus)
  
  #'  Make camera locations spatial
  eoe_cams_long <- mutate(eoe_cams_long, Setup = factor(Setup, levels = c("ungulate", "predator")))
  cams <- st_as_sf(eoe_cams_long, coords = c("Long", "Lat"), crs = wgs84)
  cams_reproj <- st_transform(cams, crs = sa_proj)
  cams_reproj <- mutate(cams_reproj, Season = ifelse(Season == "Smr20", "Summer 2020", Season),
                        Season = ifelse(Season == "Wtr20", "Winter 2020-2021", Season),
                        Season = ifelse(Season == "Smr21", "Summer 2021", Season))
  
  wolf_cams_long <- mutate(wolf_cams_long, Target = factor(Target, levels = c("Abundance", "Occupancy", "Abund_Occu")))
  wcams <- st_as_sf(wolf_cams_long, coords = c("Long", "Lat"), crs = wgs84)
  wcams_reproj <- st_transform(wcams, crs = sa_proj)
  wcams_reproj <- mutate(wcams_reproj, Season = ifelse(Season == "Smr19", "Summer 2019", Season),
                        Season = ifelse(Season == "Smr20", "Summer 2020", Season),
                        Season = ifelse(Season == "Smr21", "Summer 2021", Season))
  
  #'  Quick plot to visualize camera locations within GMUs
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus, aes(fill = NAME)) +
    scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
    geom_sf(data = cams_reproj, aes(color = Season), shape = 16) +
    scale_color_manual(values=c("#d1495b", "#edae49", "#00798c")) +
    guides(fill=guide_legend(title="GMU")) +
    ggtitle("IDFG Cameras for EoE Project") +
    coord_sf(xlim = c(-13050000, -12700000), ylim = c(5700000, 6274865), expand = TRUE) +
    theme_bw() 
  #'  Several cameras appear to have incorrect coordinates - way outside focal GMUs
  
  #'  Problem cams!
  LocationID <- c("UNKNOWN", "GMU6_P_109", "GMU6_U_122", "GMU6_U_109", "GMU6_P_63", "GMU6_U_63", "GMU10A_U_101")
  Lat <- c(48.04700, 47.67241, 47.48708, 47.48655, 46.72188, 46.72188, 46.79938)
  Long <- c(-116.9418, -116.7394, -116.7189, -116.7187, -117.0134, -117.0134, -116.5604)
  Season <- c("Smr21", "Wtr20", "Wtr20", "Wtr20", "Wtr20", "Wtr20", "Wtr20")
  Probs <- rep("Problem locations", length(LocationID))
  prob_cams <- as.data.frame(cbind(LocationID, Lat, Long, Season, Probs)) 
  # write.csv(prob_cams, file = "./Data/IDFG camera data/problem_eoe_cam_locs.csv")
  prob_cams <- st_as_sf(prob_cams, coords = c("Long", "Lat"), crs = wgs84)
  probs_reproj <- st_transform(prob_cams, crs = sa_proj)
  
  #'  Map problem cams
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus, aes(fill = NAME)) +
    scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
    geom_sf(data = cams_reproj, aes(color = Season), shape = 16) +
    geom_sf(data = probs_reproj, aes(color = Probs), shape = 16, size = 3) +
    scale_color_manual(values=c("#B90E0A", "#d1495b", "#edae49", "#00798c")) +
    guides(fill = guide_legend(title = "GMU"), color = guide_legend(title = "Camera Setup")) +
    coord_sf(xlim = c(-13050000, -12700000), ylim = c(5700000, 6274865), expand = TRUE) +
    ggtitle("Problematic EoE Cameras?") +
    theme_bw() 
  
  #'  Zoom in on each GMU
  #'  GMU 1
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME == "1",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#CC6666")) +
    geom_sf(data = cams_reproj[cams_reproj$Gmu == "1",], aes(color = Setup), shape = 16, size = 2.5) +
    scale_color_manual(values=c("#000000", "#006A67"), labels = c("Ungulate", "Predator")) +
    guides(color = guide_legend(title = "Camera Setup"), fill = "none") +
    coord_sf(xlim = c(-13050000, -12900000), ylim = c(6100000, 6274865), expand = TRUE) +
    ggtitle("IDFG Cameras for EoE Project: GMU 1") +
    theme_bw() 
  #'  GMU 6
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME == "6",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#9999CC")) +
    geom_sf(data = cams_reproj[cams_reproj$Gmu == "6",], aes(color = Setup), shape = 16, size = 2.5) +
    scale_color_manual(values=c("#000000", "#ab5e00"), labels = c("Ungulate", "Predator")) +
    guides(color = guide_legend(title = "Camera Setup"), fill = "none") +
    coord_sf(xlim = c(-13000000, -12860000), ylim = c(5920000, 6050000), expand = TRUE) +    
    ggtitle("IDFG Cameras for EoE Project: GMU 6") +
    theme_bw()
  #'  GMU 10A
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = eoe_gmus[eoe_gmus$NAME == "10A",], aes(fill = NAME)) +
    scale_fill_manual(values=c("#66CC99")) +
    geom_sf(data = cams_reproj[cams_reproj$Gmu == "10A",], aes(color = Setup), shape = 16, size = 2.5) +
    scale_color_manual(values=c("#000000", "#ab5e00"), labels = c("Ungulate", "Predator")) +
    guides(color = guide_legend(title = "Camera Setup"), fill = "none") +
    coord_sf(xlim = c(-13000000, -12800000), ylim = c(5750000, 6000000), expand = TRUE) +
    ggtitle("IDFG Cameras for EoE Project: GMU 10A") +
    theme_bw()
  
  
  #'  Map wolf locations
  #'  By season of last deployment (most deployed all three seasons)
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = wcams_reproj, aes(color = Season), shape = 16) +
    guides(color = guide_legend(title = "Season last monitored")) +
    ggtitle("IDFG Cameras for Wolf Monitoring") +
    theme_bw() +
    theme(legend.position = c(0.77, 0.85))
  #'  By targeted inference based on sampling design
  ggplot() +
    geom_sf(data = gmu) +
    geom_sf(data = wcams_reproj, aes(color = Target), shape = 16) +
    scale_color_manual(values=c("#CC6666", "#66CC99", "#004D8C"), labels = c("Abundance", "Occupancy", "Abund/Occu")) +
    guides(color = guide_legend(title = "Target Inference")) +
    ggtitle("IDFG Cameras for Wolf Monitoring") +
    theme_bw() +
    theme(legend.position = c(0.80, 0.85))
  

  ####  Summary stats on camera deployment  ####
  #'  ---------------------------------------
  #'  Number of unique EoE sites (exact coords vary by year but site basically the same)
  n_eoe <- length(unique(eoe_cams_long$NewLocationID))
  #'  Number of predator vs ungulate setup cameras
  n_eoe_P <- length(unique(eoe_cams_long$NewLocationID[eoe_cams_long$Setup == "predator"]))
  n_eoe_U <- length(unique(eoe_cams_long$NewLocationID[eoe_cams_long$Setup == "ungulate"]))
  #'  Number of cameras per season and setup type
  yr_eoe <- eoe_cams_long %>%
    group_by(Setup, Season) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(Project = "EoE",
           Season = ifelse(Season == "Wtr20", "Wtr20-21", Season)) %>%
    relocate(Project, .before = Setup)
  #'  Number of cameras per season, GMU, and setup type
  yr_gmu_eoe <- eoe_cams_long %>%
    group_by(Setup, Season, Gmu) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(Project = "EoE",
           Season = ifelse(Season == "Wtr20", "Wtr20-21", Season)) %>%
    relocate(Project, .before = Setup)
  # write.csv(yr_gmu_eoe, file = "./Data/IDFG camera data/camera_EoE_summary.csv")

  #'  Number of unique wolf sites (exact coords vary by year but site basically the same)
  n_wolf <- length(unique(wolf_cams_long$NewLocationID))
  #'  Number of abundance vs occupancy cameras (both counts include cameras that 
  #'  were deployed for both abundance & occupancy)
  n_wolf_A <- length(unique(wolf_cams_long$NewLocationID[wolf_cams_long$Target == "Abundance" | wolf_cams_long$Target == "Abund_Occu"]))
  n_wolf_O <- length(unique(wolf_cams_long$NewLocationID[wolf_cams_long$Target == "Occupancy" | wolf_cams_long$Target == "Abund_Occu"]))
  #'  Number of cameras per season and setup type
  yr_wolf <- wolf_cams_long %>%
    group_by(Target, Season) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(Project = "Wolf") %>%
    relocate(Project, .before = Target) %>%
    rename(Setup = Target)
  
  #'  Merge summary tables together and save
  camera_summary <- rbind(yr_eoe, yr_wolf) %>%
    arrange(Project, Season)
  # write.csv(camera_summary, file = "./Data/IDFG camera data/camera_EoEWolf_summary.csv")
  
  
  
  
  
  
  
  
  
  
  