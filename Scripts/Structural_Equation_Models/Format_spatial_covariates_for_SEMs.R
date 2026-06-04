  #'  ----------------------------------
  #'  Format spatial covariates for SEM
  #'  Adapted from M. Parson code
  #'  June 2026
  #'  ----------------------------------
  #'  Format spatial covariates for SEM, including:
  #'    1) Winter severity: the cumulative number of days with >15” of snow and/or 
  #'       a minimum daily temperature < 0°F (DeLgiudice et al. 2002). This combines
  #'       SNODAS and Daymet data bases. Annual winter severity calculated at the 
  #'       GMU level to account for migration between summer and winter ranges. 
  #'    2) Road density - proxy for harvest exposure
  #'    3) Proportion of public land - proxy for harvest exposure
  #'    4) Species-specific harvest. Deer, elk, moose, lion, and bear harvest at
  #'       the GMU scale owing to level of reporting data. Wolf harvest at the 
  #'       cluster level. Coyote harvest data was unavailable.
  #'    5) Percent disturbed forest in the previous 20 years.  
  
  # Load packages
  library(terra)
  library(sf)
  library(tidyverse)
  library(tigris)
  library(mapview)
  
  #'  Read in cluster polygon shapefiles
  gmu1_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu1_05.11.26.shp") %>% mutate(GMU = "GMU1")
  gmu6_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu6_05.11.26.shp") %>% mutate(GMU = "GMU6")
  gmu10a_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu10a_05.11.26.shp") %>% mutate(GMU = "GMU10A")
  
  #'  Merge cluster polygons across GMUs
  cluster_poly <- bind_rows(gmu10a_poly, gmu1_poly, gmu6_poly) %>%
    rename(ClusterID = Clusters) %>%
    mutate(Cluster_unique = 1:nrow(.)) %>%
    relocate(Cluster_unique, .before = ClusterID) 
  mapview::mapview(cluster_poly, zcol = "Cluster_unique")
  
  #'  Read in EoE GMUs shapefile
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Load GMUs to get total area 
  idfg_gmus <- read_sf("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
  #'  Convert to a data frame and tidy up
  gmu_areas <- data.frame(GMU = idfg_gmus$NAME, area_km2 = as.numeric(st_area(idfg_gmus)/1000000)) %>% 
    filter(GMU %in% c("1","6","10A")) %>% 
    mutate(GMU = case_when(
      GMU == "1" ~ "GMU1",
      GMU == "6" ~ "GMU6",
      GMU == "10A" ~ "GMU10A"
    ))
  #'  Now make a spatial version
  gmu_areas_sf <- idfg_gmus %>%
    filter(NAME == "1" | NAME == "6" | NAME == "10A") %>%
    rename(GMU = NAME) %>%
    mutate(GMU = paste0("GMU", GMU)) %>%
    dplyr::select(GMU) %>%
    full_join(gmu_areas, by = "GMU") %>%
    relocate(area_km2, .after = "GMU")
  
  #'  ----------------------------------
  ####  Generate Winter Severity Index  ####
  #'  ----------------------------------
  #'  Read in snow depth data
  snow19 <- rast("./Shapefiles/IDFG spatial data 2023 additional/SNODAS/extracted_Depth_from_2019-01-01_to_2019-12-31.tiff")
  snow20 <- rast("./Shapefiles/IDFG spatial data 2023 additional/SNODAS/extracted_Depth_from_2020-01-01_to_2020-12-31.tiff")
  snow21 <- rast("./Shapefiles/IDFG spatial data 2023 additional/SNODAS/extracted_Depth_from_2021-01-01_to_2021-12-31.tiff")
  snow22 <- rast("./Shapefiles/IDFG spatial data 2023 additional/SNODAS/extracted_Depth_from_2022-01-01_to_2022-12-31.tiff")
  snow23 <- rast("./Shapefiles/IDFG spatial data 2023 additional/SNODAS/extracted_Depth_from_2023-01-01_to_2023-12-31.tiff")
  
  #'  Read in temperature data
  tmin19 <- rast("./Shapefiles/IDFG spatial data 2023 additional/Daymet/compiled/2019tmin.tiff")
  tmin20 <- rast("./Shapefiles/IDFG spatial data 2023 additional/Daymet/compiled/2020tmin.tiff")
  tmin21 <- rast("./Shapefiles/IDFG spatial data 2023 additional/Daymet/compiled/2021tmin.tiff")
  tmin22 <- rast("./Shapefiles/IDFG spatial data 2023 additional/Daymet/compiled/2022tmin.tiff")
  tmin23 <- rast("./Shapefiles/IDFG spatial data 2023 additional/Daymet/compiled/2023tmin.tiff")
  
  tmin19 <- terra::project(tmin19,snow19)
  tmin20 <- terra::project(tmin20,snow20)
  tmin21 <- terra::project(tmin21,snow21)
  tmin22 <- terra::project(tmin22,snow22)
  tmin23 <- terra::project(tmin23,snow23)
  
  #'  Compile snow data for each winter from October 1 through May 1
  win20_depth <- c(snow19[[274:365]],snow20[[1:122]]) # Snow depth from October 1 2019 - May 2020
  win21_depth <- c(snow20[[275:366]],snow21[[1:121]]) # Snow depth from October 1 2020 - May 2021
  win22_depth <- c(snow21[[274:365]],snow22[[1:121]]) # Snow depth from October 1 2021 - May 2022
  win23_depth <- c(snow22[[274:365]],snow23[[1:121]]) # Snow depth from October 1 2022 - May 2023
  win_depth <- list(win20_depth, win21_depth, win22_depth, win23_depth)
  
  #'  Compile min temperature data for each winter from October 1 through May 1
  win20_temp <- c(tmin19[[274:365]],tmin20[[1:121]]) # Winter temps from October 1 2019 - May 2020
  win21_temp <- c(tmin20[[274:365]],tmin21[[1:121]]) # Winter temps from October 1 2020 - May 2021
  win22_temp <- c(tmin21[[274:365]],tmin22[[1:121]]) # Winter temps from October 1 2021 - May 2022
  win23_temp <- c(tmin22[[274:365]],tmin23[[1:121]]) # Winter temps from October 1 2022 - May 2023
  win_temp <- list(win20_temp, win21_temp, win22_temp, win23_temp)
  
  #'  Create summed 1/0 rasters of whether snow depth is > 38cm and temp is < 17.7C
  win_depth_thresh <- function(depth_rast) {
    #'  Reclassify raster values as 0 if snow depth <= 380 mm, 1 if > 38 mm and sum 1's
    win_depth_threshold <- sum(classify(depth_rast, rcl = matrix(c(0,380,0,380,5000,1), byrow = T, nrow = 2)))
    return(win_depth_threshold)
  }
  win_depth_classified <- lapply(win_depth, win_depth_thresh)
  
  win_temp_thresh <- function(temp_rast) {
    #'  Reclassify raster values as 1 if temp <= -17.7C, 0 if > -17.7C and sum 1's
    win_temp_threshold <- sum(classify(temp_rast, rcl = matrix(c(-100,-17.7,1,-17.7,100,0), byrow = T, nrow = 2)))
    return(win_temp_threshold)
  }
  win_temp_classified <- lapply(win_temp, win_temp_thresh)
  
  #'  Calculate Winter Severity Index based on DelGuidice et al. (2002)
  #'  Consider changing - large effected by snow depth with minimal temperature influence
  win20_WSI <- win_depth_classified[[1]] + win_temp_classified[[1]]
  win21_WSI <- win_depth_classified[[2]] + win_temp_classified[[2]]
  win22_WSI <- win_depth_classified[[3]] + win_temp_classified[[3]]
  win23_WSI <- win_depth_classified[[4]] + win_temp_classified[[4]]
  
  #'  Visualize
  plot(win20_depth[[1]])
  plot(win20_temp[[1]])
  plot(win20_WSI)
  
  #'  Reproject GMUs to match WSI
  gmu_areas_sf_reporj <- st_transform(gmu_areas_sf, crs = crs(win20_WSI))
  
  #'  Calculate WSI at the GMU scale
  GMU_WSI <- data.frame(GMU = gmu_areas_sf,
                        WSI_20_gmu = as.numeric(zonal(win20_WSI, vect(gmu_areas_sf_reporj), na.rm = T)$sum),
                        WSI_21_gmu = as.numeric(zonal(win21_WSI, vect(gmu_areas_sf_reporj), na.rm = T)$sum),
                        WSI_22_gmu = as.numeric(zonal(win22_WSI, vect(gmu_areas_sf_reporj), na.rm = T)$sum),
                        WSI_23_gmu = as.numeric(zonal(win23_WSI, vect(gmu_areas_sf_reporj), na.rm = T)$sum)) %>%
    dplyr::select(-GMU.geometry) %>%
    rename("GMU" = "GMU.GMU")
  
  #'  -----------------------------------------------------
  ####  Percent disturbed forest within previous 20 years  #### 
  #'  -----------------------------------------------------
  #'  Read in Hansen forest disturbance layer
  hist_dist <- rast("./Shapefiles/IDFG spatial data 2023 additional/Hansen_GFC-2023-v1.11_lossyear_50N_120W.tif")
  hist_dist <- terra::project(hist_dist,snow19,method = "near")
  
  # 0 = undisturbed 
  # 1 = disturbed in 2001
  # ...
  # 23 = disturbed in 2023
  
  #'  Reclassify forest disturbance to area disturbed in the last 20 years for each year
  hdist_20 <- classify(hist_dist, rcl = matrix(c(0,0.9,0,0.95,20.5,1,20.6,25,0), byrow = T, nrow = 3)) # disturbed 2001 - 2020
  hdist_21 <- classify(hist_dist, rcl = matrix(c(0,1.9,0,1.95,21.5,1,21.6,25,0), byrow = T, nrow = 3)) # disturbed 2002 - 2021
  hdist_22 <- classify(hist_dist, rcl = matrix(c(0,2.9,0,2.95,22.5,1,22.6,25,0), byrow = T, nrow = 3)) # disturbed 2003 - 2022
  hdist_23 <- classify(hist_dist, rcl = matrix(c(0,3.9,0,3.95,23.5,1,23.6,25,0), byrow = T, nrow = 3)) # disturbed 2004 - 2023
  
  #'  Visualize
  plot(hdist_20)
  plot(hdist_23)
  
  #'  ----------------
  ####  Road Density  ####
  #'  ----------------
  #'  Load road density data
  roaddens <- rast("./Shapefiles/IDFG spatial data 2023 additional/AllRoads_2022_density_kmkm2.tif")
  
  #'  Get spatial extent of needed raster area
  xmin_rast <- xmin(win20_WSI)
  xmax_rast <- xmax(win20_WSI)
  ymin_rast <- ymin(win20_WSI)
  ymax_rast <- ymax(win20_WSI)
  
  #'  Create template raster that covers that spatial extent
  a <- rast(xmin = xmin_rast,
            xmax = xmax_rast,
            ymin = ymin_rast,
            ymax = ymax_rast,
            crs = win20_WSI)
  #'  Project the template raster to the CRS of the road density layer
  a2 <- project(a,crs(roaddens))
  
  #'  Crop road density to the smaller area, then project to same CRS as cluster polygons
  roaddens <- terra::crop(roaddens, a2)
  roaddens <- terra::project(roaddens, cluster_poly)
  plot(roaddens)
  
  #'  -----------------------------
  ####  Proportion of public land  ####
  #'  -----------------------------
  #'  Load land ownership data
  ownership <- read_sf("./Data/IDFG BGMR data/IDFG 2023 harvest data/IDOwnership/IDOwnership.shp")
  ownership <- ownership %>% 
    mutate(public = case_when(
      AGNCY_NAME %in% c("BLM","BOR","NWR","STATE","STATEFG","USFS") ~ 1,
      TRUE ~ 0)) %>%
    st_transform(., crs(roaddens))
  
  #'  Rasterize
  ownership_rast <- (rasterize(vect(ownership), field = "public", roaddens))
  
  #'  Visualize
  plot(ownership_rast)
  
  #'  ----------------
  ####  Harvest data  ####
  #'  ----------------
  #####  Deer and elk harvest  #####
  elk_cont <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Elk_controlled_2019-2024_1-6-10a.csv")%>% 
    dplyr::select(-CHunt) %>% 
    rename(Unit = Hunt.Area)
  elk_gen <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Elk_general_2019-2024_1-6-10a.csv")%>% 
    dplyr::select(-Method)
  deer_cont <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Deer_controlled_2019-2024_1-6-10a.csv") %>% 
    dplyr::select(-CHunt)%>% 
    rename(Unit = Hunt.Area)
  deer_gen <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Deer_general_2019-2024_1-6-10a.csv")%>% 
    dplyr::select(-Method)

  #'  Total elk harvest by year, unit, and harvest type
  elk_harv <- bind_rows(elk_cont,elk_gen) %>% 
    group_by(Year,Unit) %>% 
    summarize(Antlered_harvest = sum(Antlered, na.rm = T),
              Antlerless_harvest = sum(Antlerless, na.rm = T),
              Total_harvest = sum(Harvest, na.rm = T),
              Total_days = sum(Days, na.rm = T))
  #'  Total deer harvest by year, unit, and harvest type
  deer_harv <- bind_rows(deer_cont,deer_gen) %>% 
    group_by(Year,Unit) %>% 
    summarize(Antlered_harvest = sum(Antlered, na.rm = T),
              Antlerless_harvest = sum(Antlerless, na.rm = T),
              Total_harvest = sum(Harvest, na.rm = T),
              Total_days = sum(Days, na.rm = T))
  
  #'  Create lookup lists to help match harvest units to GMUs
  elk_lookup <- list(c("1","1-1","1-2"),
                     c("6","1-2"),
                     c("10A","10A-1","10A-1X"))
  names(elk_lookup) <- c("1","6","10A")
  deer_lookup <- list(c("1","1-1","1-1x"),
                      c("6","1-1"),
                      c("10A","10A-1X"))
  names(deer_lookup) <- c("1","6","10A")
  
  #'  Function to match harvest unit to GMU and calculate catch per unit effort
  x <- function(lookup, data){
    data %>% 
      filter(Unit %in% lookup) %>% 
      group_by(Year) %>% 
      summarise(Antlered_harvest = sum(Antlered_harvest),
                Antlerless_harvest = sum(Antlerless_harvest),
                Total_harvest = sum(Total_harvest),
                Total_effort = sum(Total_days)) %>%
      ungroup() %>%
      #'  Calculate harvest per hunter day
      mutate(Harvest_per_HunterDay = round(Total_harvest / Total_effort, 4))
  }
  #'  Format elk and deer data per GMU
  elkharv_list <- lapply(elk_lookup,x,elk_harv)
  elkharv_list[[1]]$GMU <- "GMU1"
  elkharv_list[[2]]$GMU <- "GMU6"
  elkharv_list[[3]]$GMU <- "GMU10A"
  deerharv_list <- lapply(deer_lookup,x,deer_harv)
  deerharv_list[[1]]$GMU <- "GMU1"
  deerharv_list[[2]]$GMU <- "GMU6"
  deerharv_list[[3]]$GMU <- "GMU10A"
  
  deer_harv <- bind_rows(deerharv_list)
  elk_harv <- bind_rows(elkharv_list)
  
  #####  Moose, Lion, Bear harvest  #####
  bear_harv <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Blackbear_hunt_2019-2024_1-6-10a.csv")
  moose_harv <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Moose_hunt_2019-2024_1-6-10a.csv")
  lion_harv <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/Lion_hunt_2019-2024_1-6-10a.csv")
  bear_moose_lion <- list(bear_harv, moose_harv, lion_harv)
  
  #'  Function to summarize harvest data
  spp_harv <- function(harv_data) {
    harvest <- harv_data %>%
      group_by(Year, Unit) %>%
      summarise(Total_harvest = sum(Harvest), 
                Total_effort = sum(Days)) %>%
      mutate(Harvest_per_HunterDay = round(Total_harvest / Total_effort, 4),
             GMU = case_when(
               Unit == '1' ~ "GMU1",
               Unit == "10A" ~ "GMU10A",
               Unit == "6" ~ "GMU6")) %>%
      ungroup()
  } 
  other_harv <- lapply(bear_moose_lion, spp_harv)
  bear_harv <- other_harv[[1]]
  moose_harv <- other_harv[[2]]
  lion_harv <- other_harv[[3]]
  
  #####  Join elk, deer, bear, moose, and lion data  #####
  #'  Pivot harvest data to wider format 
  elk_harv_wide <- elk_harv %>% 
    left_join(gmu_areas, by = "GMU") %>%
    rename(Elk_antlered = Antlered_harvest,
           Elk_antlerless = Antlerless_harvest,
           Elk_total = Total_harvest) %>% 
    mutate(Elk_total_per100km = round(Elk_total / area_km2*100, 3)) %>% 
    dplyr::select(Year, GMU, Elk_total_per100km) %>% 
    pivot_wider(names_from = Year, values_from = Elk_total_per100km, names_prefix = "Elk_per100km_") 
  deer_harv_wide <- deer_harv %>% 
    left_join(gmu_areas, by = "GMU") %>% 
    rename(Deer_antlered = Antlered_harvest,
           Deer_antlerless = Antlerless_harvest,
           Deer_total = Total_harvest) %>% 
    mutate(Deer_total_per100km = round(Deer_total / area_km2*100, 3)) %>% 
    dplyr::select(Year,GMU,Deer_total_per100km) %>% 
    pivot_wider(names_from = Year, values_from = Deer_total_per100km, names_prefix = "Deer_per100km_") 
  bear_harv_wide <- bear_harv %>% 
    left_join(gmu_areas, by = "GMU") %>% 
    rename(Bear_total = Total_harvest) %>%
    mutate(Bear_total_per100km = round(Bear_total / area_km2*100, 3)) %>% 
    dplyr::select(Year, GMU, Bear_total_per100km) %>% 
    pivot_wider(names_from = Year, values_from = Bear_total_per100km, names_prefix = "Bear_per100km_") 
  lion_harv_wide <- lion_harv %>% 
    left_join(gmu_areas, by = "GMU") %>% 
    rename(Lion_total = Total_harvest) %>%
    mutate(Lion_total_per100km = round(Lion_total / area_km2*100, 3)) %>% 
    dplyr::select(Year, GMU, Lion_total_per100km) %>% 
    pivot_wider(names_from = Year, values_from = Lion_total_per100km, names_prefix = "Lion_per100km_")
  moose_harv_wide <- moose_harv %>% 
    left_join(gmu_areas, by = "GMU") %>% 
    rename(Moose_total = Total_harvest) %>%
    mutate(Moose_total_per100km = round(Moose_total / area_km2*100, 3)) %>% 
    dplyr::select(Year, GMU, Moose_total_per100km) %>% 
    pivot_wider(names_from = Year, values_from = Moose_total_per100km, names_prefix = "Moose_per100km_")
  
  #####  Wolf harvest  #####
  #'  Load BGMR wolf harvest data
  BGMR_wolf_dat <- read.csv("Data/IDFG BGMR data/IDFG 2023 harvest data/BGMRWolfKills_2026-03-19.csv") 
   
  #'  Select needed columns and clean up kill locations
  BGMRdat <- BGMR_wolf_dat %>%
    dplyr::select(BGMR_GUID, BGMRId, unit, sex, age, killLocation, speciesCode, 
                  speciesCommonName, killDate, checkDate, gpsSystem, gpsDatum, 
                  gpsLat, gpsLong, harvestYear, dataAnalysisUnit, daysHunted) %>%
    mutate(killLocation = toupper(killLocation)) %>%
    mutate(killLocation_cleaned = case_when(
      grepl("JIM FORD",killLocation) ~ "JIM FORD CREEK",
      grepl("LASTA",killLocation) ~ "ATLASTA PEAK",
      grepl("BALDY",killLocation) ~ "BALDY BEAR",
      grepl("BERRY CREEK",killLocation) ~ "BERRY CREEK",
      grepl("BLACKWELL",killLocation) ~ "BLACKWELL DIVIDE",
      grepl("BOULDER CREEK",killLocation) ~ "BOULDER CREEK",
      grepl("BREAKFAST",killLocation) ~ "BREAKFAST CREEK",
      grepl("BRUSH LAKE",killLocation) ~ "BRUSH LAKE",
      grepl("CAMP 9|CAMP NINE",killLocation) ~ "CAMP NINE",
      grepl("CANAL CREEK",killLocation) ~ "CANAL CREEK",
      grepl("CEDAR CREEK",killLocation) ~ "CEDAR CREEK",
      grepl("CEMETERY|CEMETARY",killLocation) ~ "CEMETERY RIDGE",
      grepl("CLARK FORK",killLocation) ~ "CLARK FORK",
      grepl("DEER CREEJ|DEER CREEK",killLocation) ~ "DEER CREEK",
      grepl("EMERALD",killLocation) ~ "EMERALD CREEK",
      grepl("FISH HOOK|FISHHOOK",killLocation) ~ "FISH HOOK CREEK",
      grepl("GOAT",killLocation) ~ "GOAT MOUNTAIN",
      grepl("HELL ROARING|HELLROARING",killLocation) ~ "HELL ROARING CREEK",
      grepl("JERU|JERV",killLocation) ~ "JERU",
      grepl("KATKA|KATICA",killLocation) ~ "KATKA",
      grepl("MAGGIE CREEK",killLocation) ~ "MAGGIE CREEK",
      grepl("ST JOE|ST. JOE|JOE",killLocation) ~ "NORTH FORK ST. JOE",
      grepl("PHILIPS DRAW|PHILLIPS DRAW",killLocation) ~ "PHILLIPS DRAW",
      grepl("QUEEN",killLocation) ~ "QUEEN MOUNTAIN",
      grepl("SHANGHAI",killLocation) ~ "SHANGHAI",
      grepl("SWAMP CREEK",killLocation) ~ "SWAMP CREEK",
      grepl("TUNNEL",killLocation) ~ "TUNNEL CREEK",
      grepl("WASHINGTON",killLocation) ~ "WASHINGTON CREEK",
      grepl("WHISKEY|WHISKY",killLocation) ~ "WHISKEY CREEK",
      TRUE ~ killLocation)) %>%
    dplyr::select(BGMRId, unit, harvestYear, killDate, gpsLat, gpsLong, dataAnalysisUnit, 
                  killLocation, killLocation_cleaned, daysHunted)
  
  #'  Read in manually edited kill locations with approximate lat/lon kill locations
  wolf_morts_LL <- read.csv("./Data/IDFG BGMR data/IDFG 2023 harvest data/WolfKillLocations_LatLon.csv")
  #'  Filter out observations where lat/long kill location could not be approximated
  wolf_morts_LL <- wolf_morts_LL %>% 
    dplyr::select(BGMRId, unit, dataAnalysisUnit, killLocation, Lat, Long) %>% 
    filter(!is.na(Lat))
  
  #'  Join kills to the BGMR data
  BGMRdat <- BGMRdat %>% 
    left_join(wolf_morts_LL) %>% 
    filter(!is.na(Lat)) %>%
    arrange(unit, killDate) #harvestYear
  
  #'  Convert to sf object for plotting and intersecting
  BGMRdat <- st_as_sf(BGMRdat,coords = c("Long","Lat"), crs = 4326) # WGS 84
  
  #'  Plot to visualize points and check for GMU consistency
  plot(win20_WSI, ylim = c(48, 49.2))
  lines(cluster_poly, col = as.factor(cluster_poly$GMU), lwd = 2)
  points(BGMRdat, col = as.factor(BGMRdat$unit), cex = 2)
  
  #'  Intersect points and polygons to tally wolf harvest
  sf::sf_use_s2(FALSE)
  wolfharvest <- (st_intersection(BGMRdat, cluster_poly))
  
  #'  Wolf harvest is essentially year-round so need to reclassify harvest date
  #'  based on when it occurred relative to summer camera trapping, not calendar year.
  #'  Wolf harvest that occurs the previous summer, fall, winter, and spring are
  #'  expected to affect the current summer wolf population
  wolfharvest_adj <- wolfharvest %>%
    mutate(cameraYear = ifelse(killDate >= "2018-06-01" & killDate < "2019-06-01", 2019, 2018),
           cameraYear = ifelse(killDate >= "2019-06-01" & killDate < "2020-06-01", 2020, cameraYear),
           cameraYear = ifelse(killDate >= "2020-06-01" & killDate < "2021-06-01", 2021, cameraYear),
           cameraYear = ifelse(killDate >= "2021-06-01" & killDate < "2022-06-01", 2022, cameraYear),
           cameraYear = ifelse(killDate >= "2022-06-01" & killDate < "2023-06-01", 2023, cameraYear),
           cameraYear = ifelse(killDate >= "2023-06-01" & killDate < "2024-06-01", 2024, cameraYear),
           cameraYear = ifelse(killDate >= "2024-06-01" & killDate < "2025-06-01", 2025, cameraYear)) %>%
    relocate(cameraYear, .after = killDate) %>%
    arrange(unit, killDate)
  
  #'  Tally the number of wolves harvested in each cluster in each year (relative
  #'  to camera deployment)
  wolf_harv <- wolfharvest_adj %>% 
    st_drop_geometry() %>% 
    group_by(cameraYear, Cluster_unique) %>% #harvestYear
    summarize(n_harvest = n(),
              n_effort = sum(daysHunted, na.rm = TRUE),
              #'  Change n_effort = 0 to 1 so each reported harvest has at least
              #'  one day of effort (can't have less!). Clearly there is under-reporting
              #'  in the wolf harvest effort data.
              n_effort = ifelse(n_effort == 0, 1, n_effort)) %>% 
    mutate(cameraYear = paste0("wolfharvest_",cameraYear),
           Harvest_per_HunterDay = round(n_harvest / n_effort, 4)) 
  
  #'  Pivot data to create one column per year with total harvest per cluster
  #'  REMEMBER: wolfharvest_YEAR is actually the camera year it applies to - a one
  #'  year time lag is already built into these data based on which harvest events
  #'  were included in the total.
  harvestbyyear <- wolf_harv %>%
    dplyr::select(-c(n_effort, Harvest_per_HunterDay)) %>%
    pivot_wider(names_from = cameraYear,
                values_from = n_harvest) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    arrange(Cluster_unique)
  harvestbyyear
  #'  Pivot data to create one column per year with total harvest per hunter day per cluster
  harvestbyyear_cpu <- wolf_harv %>%
    dplyr::select(-c(n_effort, n_harvest)) %>%
    pivot_wider(names_from = cameraYear,
                values_from = Harvest_per_HunterDay) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    arrange(Cluster_unique)
  harvestbyyear_cpu
  #'  Problem with this approach is because hunting days was not reported by some
  #'  and so n_effort adjusted to 1 for those harvest events, CPU is artificially
  #'  high in those clusters compared to CPU in clusters where hunting days was
  #'  actually reported. This is likely a misleading metric of hunting intensity.
  
  #'  Join harvest data to cluster polygons
  cluster_poly_wolf <- cluster_poly %>% 
    left_join(harvestbyyear, by = "Cluster_unique") %>% 
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    arrange(Cluster_unique) %>%
    dplyr::select(-c(wolfharvest_2024, wolfharvest_2025))
  
  head(cluster_poly_wolf, 15)
  
  #'  ------------------------------------------
  ####  Extract / merge all spatial covariates  ####
  #'  ------------------------------------------
  #'  Calculate proportion disturbed forest, road density, and proportion of public
  #'  land for each cluster polygon  -----   THIS TAKES A HOT SECOND
  cluster_poly_covs <- cluster_poly_wolf %>% 
    mutate(#WSI20 = c(zonal(win20_WSI, vect(cluster_poly))$sum),
           # WSI21 = c(zonal(win21_WSI, vect(cluster_poly))$sum),
           # WSI22 = c(zonal(win22_WSI, vect(cluster_poly))$sum),
           # WSI23 = c(zonal(win23_WSI, vect(cluster_poly))$sum),
           distForest20 = c(zonal(hdist_20, vect(cluster_poly))$Layer_1),
           distForest21 = c(zonal(hdist_21, vect(cluster_poly))$Layer_1),
           distForest22 = c(zonal(hdist_22, vect(cluster_poly))$Layer_1),
           distForest23 = c(zonal(hdist_23, vect(cluster_poly))$Layer_1),
           roaddens = c(zonal(roaddens, vect(cluster_poly))$AllRoads_2022_density_kmkm2),
           proppub = c(zonal(ownership_rast, vect(cluster_poly), na.rm = T)$public))
  
  #'  Weight cluster-specific covaraites by cluster area and average within each GMU 
  cluster_poly_covs_df <- cluster_poly_covs %>% 
    st_drop_geometry() %>% 
    as.data.frame() %>% 
    group_by(Cluster_unique) %>% 
    mutate(total_area = sum(unique(area_km2))) %>% 
    # mutate(# weighted_WSI20 = round(area_km2 * WSI20, 2),
    #        # weighted_WSI21 = round(area_km2 * WSI21, 2),
    #        # weighted_WSI22 = round(area_km2 * WSI22, 2),
    #        # weighted_WSI23 = round(area_km2 * WSI23, 2),
    #        weighted_dist20 = round(area_km2 * dist20, 2),   # I don't understand the point of this since the next step leads to the original value
    #        weighted_dist21 = round(area_km2 * dist21, 2),
    #        weighted_dist22 = round(area_km2 * dist22, 2),
    #        weighted_dist23 = round(area_km2 * dist23, 2),
    #        weighted_roaddens = round(area_km2 * roaddens, 2),
    #        weighted_proppub = round(area_km2 * proppub, 2)) %>% 
    # ungroup() %>% 
    # dplyr::select(-c("area_km2","dist20","dist21","dist22","dist23")) %>%  #"WSI20","WSI21","WSI22","WSI23",
    # group_by(Cluster_unique,GMU,wolfharvest_2018,wolfharvest_2019,wolfharvest_2020,wolfharvest_2021,wolfharvest_2022,wolfharvest_2023,total_area) %>% 
    # summarise(# weighted_WSI20 = sum(weighted_WSI20) / total_area,
    #           # weighted_WSI21 = sum(weighted_WSI21) / total_area,
    #           # weighted_WSI22 = sum(weighted_WSI22) / total_area,
    #           # weighted_WSI23 = sum(weighted_WSI23) / total_area,
    #           weighted_dist20 = sum(weighted_dist20) / total_area,
    #           weighted_dist21 = sum(weighted_dist21) / total_area,
    #           weighted_dist22 = sum(weighted_dist22) / total_area,
    #           weighted_dist23 = sum(weighted_dist23) / total_area,
    #           weighted_roaddens = sum(weighted_roaddens) / total_area,
    #           weighted_proppub = sum(weighted_proppub) / total_area) %>% 
    # # slice(1) %>% 
    # ungroup() %>% 
    mutate(wolfharvest_2018_per100km = round(wolfharvest_2018 / total_area*100,3),
           wolfharvest_2019_per100km = round(wolfharvest_2019 / total_area*100,3),
           wolfharvest_2020_per100km = round(wolfharvest_2020 / total_area*100,3),
           wolfharvest_2021_per100km = round(wolfharvest_2021 / total_area*100,3),
           wolfharvest_2022_per100km = round(wolfharvest_2022 / total_area*100,3),
           wolfharvest_2023_per100km = round(wolfharvest_2023 / total_area*100,3))
  
  #'  Join all harvest data to larger covariate data frame and reduce to only necessary columns
  cluster_poly_covs_df <- cluster_poly_covs_df %>% 
    left_join(elk_harv_wide, by = "GMU") %>% 
    left_join(deer_harv_wide, by = "GMU") %>% 
    left_join(moose_harv_wide, by = "GMU") %>% 
    left_join(bear_harv_wide, by = "GMU") %>% 
    left_join(lion_harv_wide, by = "GMU") %>%
    dplyr::select(Cluster_unique, GMU, distForest20:Lion_per100km_2024) %>%
    left_join(GMU_WSI, by = "GMU")
  
  #' #'  Save
  #' write.csv(cluster_poly_covs_df, "./Data/SEM_covariates.csv", row.names = F)
  
  #'  --------------------------------------------------------
  ####  NOTE: These covariates do not need to be time lagged  ####
  #'  --------------------------------------------------------
  #'  wolfharvest_Year is the total wolves harvested from June 1 year-1 to May 30 year
  #'  and reflects the camera year it should be applied to. A one year time lag is 
  #'  already built into these data based on which harvest events were included.
  #'  WSIYear is the calculated winter severity index from October 1 year-1 to May 1 year. 
  #'  weighted_distyear is the calculated proportion of disturbed forest 20 years prior to year. 
  print("Do not time lag the following variables: wolfharvest, WSI, weighted_distyear")
  
