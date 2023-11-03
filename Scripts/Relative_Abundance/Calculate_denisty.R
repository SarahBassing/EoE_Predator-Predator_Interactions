  #'  ----------------------------------------
  #'  Calculate denisty with TIFC method
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  September 2023
  #'  ----------------------------------------
  #'  Calculate density using time-in-front-of-camera methods following analyses
  #'  described by Huggard et al. 2018, Warbington & Boyce 2020, and Becker et al. 2022.
  #'  Code modified from Becker et al. 2022 repository: https://github.com/mabecker89/tifc-method/tree/main
  #'  
  #'  Time-in-front-of-camera (TIFC) generated in Calculate_time-in-front-of-camera.R
  #'  Effective detection distance (EDD) appended to TIFC in Effective_detection_distance.R
  #'  Sequential probimgs generated in Detection_data_cleaning.R (Data_Formatting folder)
  #'  Sampling effort data generated in Detection_histories_for_occmod.R (MultiSpp_OccMod folder)
  #'  Leave probabilities taken directly from Becker et al. (2018)
  #'  ----------------------------------------
    
  #'  Load libraries
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  
  #'  -------------------------------
  ####  Load and format input data  ####
  #'  -------------------------------  
  #'  TIFC data
  load("./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time_avg_edd.RData") #eoe_all_fov-time_edd.RData
  # load("./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time.RData")
  # df_tt_full_20s <- read_csv("./Data/Relative abundance data/RAI Phase 2/eoe_all_20s_fov-time.csv") 
  # df_tt_full_21s <- read_csv("./Data/Relative abundance data/RAI Phase 2/eoe_all_21s_fov-time.csv")
  # df_tt_full_22s <- read_csv("./Data/Relative abundance data/RAI Phase 2/eoe_all_22s_fov-time.csv")
  
  #'  Filter to focal species
  focal_spp <- function(df_tt_full) {
    df_tt_full <- df_tt_full %>%
      # filter(common_name == "bear_black" | common_name == "bobcat" | common_name == "coyote" |
      #          common_name == "elk" | common_name == "moose" | common_name == "mountain_lion" |
      #          common_name == "muledeer" | common_name == "rabbit_hare" | 
      #          common_name == "whitetaileddeer" | common_name == "wolf")
      filter(common_name == "bear_black" | common_name == "bobcat" | common_name == "coyote" |
              common_name == "mountain_lion" | common_name == "wolf") %>%
      rename("detdist" = "avg_dist") %>%
      rename("detdist_v2" = "max_dist")
    return(df_tt_full)
  }
  df_tt_full_20s <- focal_spp(tt_list[[1]])
  df_tt_full_21s <- focal_spp(tt_list[[2]])
  df_tt_full_22s <- focal_spp(tt_list[[3]])
  
  #'  Load effective detection distance (EDD)
  #'  Currently using detection zone measured at camera by IDFG. Area (m2) measured
  #'  using walk test, flagging, or other methods to estimate detection zone
  #'  These values are not corrected for habitat type or species
  load("./Data/Relative abundance data/RAI Phase 2/area_m2.RData")
  area_m2_20s <- area_m2[[1]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  area_m2_21s <- area_m2[[2]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  area_m2_22s <- area_m2[[3]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  
  #'  Update species names for "distance groups" based on Becker et al. (2022) classifications
  spp_20s <- as.data.frame(unique(df_tt_full_20s$common_name)); names(spp_20s) <- "common_name"
  spp_21s <- as.data.frame(unique(df_tt_full_21s$common_name)); names(spp_21s) <- "common_name"
  spp_22s <- as.data.frame(unique(df_tt_full_22s$common_name)); names(spp_22s) <- "common_name"
  df_dist_groups_20s <- arrange(spp_20s, common_name)
  df_dist_groups_21s <- arrange(spp_21s, common_name)
  df_dist_groups_22s <- arrange(spp_22s, common_name)
  # dist_group <- c("Bear", "Lynx", "CoyoteFox","Elk", "BigForestUngulates", "BigForestCarnivores", 
  #                 "Muledeer", "SmallForest", "WTDeer", "BigForestCarnivores")
  dist_group <- c("Bear", "Lynx", "CoyoteFox","BigForestCarnivores", "BigForestCarnivores")
  df_dist_groups_20s$dist_group <- dist_group
  df_dist_groups_21s$dist_group <- dist_group
  df_dist_groups_22s$dist_group <- dist_group
  
  #'  Camera field of view angle (per Reconyx Hyperfire documentation)
  cam_fov_angle <- 42
  
  #'  --------------------------------------
  ####  Density: animal-time per area-time  ####
  #'  --------------------------------------
  
  #####  Step 1. Append EDD information  ####
  #'  -----------------------------------
  #' #'  Only need this if actually using detection distance!
  #' df_detdist <- dd %>%
  #'   as.data.frame() %>%
  #'   rownames_to_column(var = "location") %>%
  #'   gather(key = "SpGroupSeason", value = "detdist", WTDeer.summer:`Pronghorn.winter`) %>%
  #'   mutate(SpGroupSeason = str_replace_all(SpGroupSeason, "[//(//)// ]", "")) %>%
  #'   # Create two new columns: Detection Distance Group and Season, sep by "."
  #'   separate(SpGroupSeason, into = c("dist_group", "season")) %>%
  #'   mutate(dist_group = str_replace(dist_group, "wapiti", ""),
  #'          season = tolower(season))
  
  #'  Join all ingredients ('ing') required to calculate density
  all_ing <- function(df_tt_full, df_dist_groups, area_m2) {
    df_dens_ing <- df_tt_full %>%
      left_join(df_dist_groups, by = "common_name") %>%
      left_join(area_m2, by = "location") #%>%
    # left_join(df_detdist, by = c("location", "dist_group")) 
    
    return(df_dens_ing)
  }
  df_dens_ing_20s <- all_ing(df_tt_full_20s, df_dist_groups = df_dist_groups_20s, area_m2 = area_m2_20s)
  df_dens_ing_21s <- all_ing(df_tt_full_21s, df_dist_groups = df_dist_groups_21s, area_m2 = area_m2_21s)
  df_dens_ing_22s <- all_ing(df_tt_full_22s, df_dist_groups = df_dist_groups_22s, area_m2 = area_m2_22s)
  
  
  #####  Step 2. Calculate density  #####         
  #'  -----------------------------------
  #'  Density = animal-time / area-time = animals/area
  #'  D = sum(N x Tf) / (Af x To), where 
  #'  N = number of individuals, Tf = time in front of camera's field of view, so
  #'  N * Tf = total_duration (animal-time)
  #'  Af = area of camera's field-of-view (m), To = total camera operating time
  #'  Af = Area(m2) = (EDD^2 x pi x a) / 360, where a = camera angle of view, so
  #'  Af * To = sampling effort (area-time)
  calc_density <- function(df_dens_ing) {
    df_density <- df_dens_ing %>%
      mutate(effort = total_season_days * (detdist ^ 2 * pi * (cam_fov_angle / 360)) / 100,
             effort_walktest = total_season_days * (edd_area_m2 * pi * (cam_fov_angle / 360)) / 100,
             #'  Catch per unit effort
             cpue = total_duration / effort,
             cpue_walktest = total_duration / effort_walktest,
             #'  Catch per unit effort in km2 
             cpue_km2 = cpue / 60 / 60 / 24 * 10000,
             cpue_km2_walktest = cpue_walktest / 60 / 60 / 24 * 10000,
             #'  Catch per unit effort in 100km2
             cpue_100km2_og = cpue_km2 * 100,
             cpue_100km2_walktest = cpue_km2_walktest * 100,
             
             #'  Units that work for me:
             #'  Sampling effort = sampled area (m2) * total operation time (seconds) 
             effort_m2_sec = (total_season_days*60*60*24) * (detdist ^ 2 * pi * (cam_fov_angle / 360)),
             effort_m2_sec_v2 = (total_season_days*60*60*24) * (detdist_v2 ^ 2 * pi * (cam_fov_angle / 360)),
             #'  Catch per unit effort: animal seconds per 1 m2 seconds (will differ from cpue above)
             cpue_m2_sec = total_duration / effort_m2_sec,
             cpue_m2_sec_v2 = total_duration / effort_m2_sec_v2,
             #'  Catch per unit effort in km2 = animals per 1 km2 (should match cpue_km2 above)
             cpue_km2_sec = (total_duration / effort_m2_sec) * 1000000,
             cpue_km2_sec_v2 = (total_duration / effort_m2_sec_v2) * 1000000,
             #'  Catch per unit effort in 100km2 = animals per 100 km2
             cpue_100km2 = cpue_km2_sec * 100,
             cpue_100km2_v2 = cpue_km2_sec_v2 * 100) 
    
      return(df_density)
  }
  df_density_20s <- calc_density(df_dens_ing_20s); df_density_20s$season <- "Smr20"
  df_density_21s <- calc_density(df_dens_ing_21s); df_density_21s$season <- "Smr21"
  df_density_22s <- calc_density(df_dens_ing_22s); df_density_22s$season <- "Smr22"
  
  #'  List annual tifc relative density data
  tifc_density_list <- list(df_density_20s, df_density_21s, df_density_22s)
  
  #'  Bind into a single data frame
  df_density <- rbind(df_density_20s, df_density_21s, df_density_22s) %>%
    #'  Add GMU to dataset
    mutate(gmu = sub("_.*", "", location),
           setup = ifelse(grepl("P", location), "predator", "ungulate")) %>%
    relocate(gmu, .after = "season") %>%
    relocate(setup, .after = "gmu")
  
  #'  SAVE! 
  write_csv(df_density, "./Data/Relative abundance data/RAI Phase 2/eoe_all_years_density_avg_edd_predonly.csv")
  
  
  #####  Summarize density data across years, gmus, and setup  ####
  #'  ---------------------------------------------------------
  #'  Summary stats for density estimates
  density_stats <- df_density %>%
    group_by(gmu, common_name) %>%
    summarise(mean_density_km2 = round(mean(cpue_km2_sec, na.rm = TRUE), 2),
              se_density_km2 = round((sd(cpue_km2_sec, na.rm = TRUE)/sqrt(nrow(.))), 2),
              mean_density_100km2 = round(mean(cpue_100km2, na.rm = TRUE), 2),
              se_density_100km2 = round((sd(cpue_100km2, na.rm = TRUE)/sqrt(nrow(.))), 2)) %>%
    ungroup() %>%
    arrange(gmu, mean_density_km2)
  print(density_stats)
  write_csv(density_stats, "./Data/Relative abundance data/RAI Phase 2/tifc_density_stats_avg_edd_predonly.csv")
  
  
  #'  Summary stats by gmu, year, and camera set up
  more_density_stats <- df_density %>%
    group_by(gmu, common_name, setup, season) %>%
    summarise(mean_density_100km2 = round(mean(cpue_100km2, na.rm = TRUE), 2),
              se = round((sd(cpue_100km2, na.rm = TRUE)/sqrt(nrow(.))), 2)) %>%
    ungroup()
  predcam_density <- filter(more_density_stats, setup == "predator") %>% 
    dplyr::select(-setup) %>%
    rename("predator_cam_density_100km2" = "mean_density_100km2") %>%
    rename("predator_cam_se" = "se")
  ungcam_density <- filter(more_density_stats, setup == "ungulate") %>%
    dplyr::select(-setup) %>%
    rename("ungulate_cam_density_100km2" = "mean_density_100km2") %>%
    rename("ungulate_cam_se" = "se")
  density_by_setup <- full_join(predcam_density, ungcam_density, by = c("gmu", "common_name", "season"))
  
  gmu1_density_stats <- filter(density_by_setup, gmu == "GMU1") %>%
    arrange(season, predator_cam_density_100km2) %>%
    mutate(road_density_100km2 = paste0(predator_cam_density_100km2, " (", predator_cam_se, ")"),
           random_density_100km2 = paste0(ungulate_cam_density_100km2, " (", ungulate_cam_se, ")")) %>%
    dplyr::select(common_name, season, road_density_100km2, random_density_100km2)
  print(gmu1_density_stats)
  gmu6_density_stats <- filter(density_by_setup, gmu == "GMU6") %>%
    arrange(season, predator_cam_density_100km2) %>%
    mutate(road_density_100km2 = paste0(predator_cam_density_100km2, " (", predator_cam_se, ")"),
           random_density_100km2 = paste0(ungulate_cam_density_100km2, " (", ungulate_cam_se, ")")) %>%
    dplyr::select(common_name, season, road_density_100km2, random_density_100km2)
  print(gmu6_density_stats)
  gmu10A_density_stats <- filter(density_by_setup, gmu == "GMU10A") %>%
    arrange(season, predator_cam_density_100km2) %>%
    mutate(road_density_100km2 = paste0(predator_cam_density_100km2, " (", predator_cam_se, ")"),
           random_density_100km2 = paste0(ungulate_cam_density_100km2, " (", ungulate_cam_se, ")")) %>%
    dplyr::select(common_name, season, road_density_100km2, random_density_100km2)
  print(gmu10A_density_stats)
  
  gmu_density_stats <- rbind(gmu1_density_stats, gmu6_density_stats, gmu10A_density_stats)
  write_csv(gmu_density_stats, "./Data/Relative abundance data/RAI Phase 2/tifc_density_stats_by_gmu_predonly.csv")
  
  
  #####  Visualize TIFC density data  #####
  #'  --------------------------------
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
  
  #'  Append TIFC relative density data to spatial data
  spatial_tifc <- function(tifc, spp, cams) {
    #'  Filter tifc data to single species
    single_spp_tifc <- tifc %>%
      filter(common_name == spp) %>%
      #'  Rename camera location column to match spatial data
      rename("NwLctID" = "location")
    
    #'  Join spatial data with tifc data
    tifc_shp <- full_join(cams, single_spp_tifc, by = "NwLctID") %>%
      filter(!is.na(common_name)) %>%
      dplyr::select(c("NwLctID", "common_name", "season", "Gmu", "cpue_m2_sec", 
                      "cpue_km2_sec", "cpue_100km2", "geometry")) %>%
      #'  Change camera location column back to something less awkward
      rename("location" = "NwLctID") %>%
      #'  Remove duplicate observations (happens in the TIFC stage when multiple 
      #'  walktest viewsheds are reported for same site)
      distinct()
    
    return(tifc_shp)
  }
  spatial_tifc_bear <- mapply(tifc = tifc_density_list, spatial_tifc, spp = "bear_black", cams = cam_list, SIMPLIFY = FALSE)
  spatial_tifc_bob <- mapply(tifc = tifc_density_list, spatial_tifc, spp = "bobcat", cams = cam_list, SIMPLIFY = FALSE)
  spatial_tifc_coy <- mapply(tifc = tifc_density_list, spatial_tifc, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  spatial_tifc_lion <- mapply(tifc = tifc_density_list, spatial_tifc, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  spatial_tifc_wolf <- mapply(tifc = tifc_density_list, spatial_tifc, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  Function to map TIFC relative density index
  map_tifc <- function(sf_tifc, spp) {
    #'  Filter spatial tifc data by study area
    sf_tifc_gmu1 <- sf_tifc[sf_tifc$Gmu == "1",]
    sf_tifc_gmu6 <- sf_tifc[sf_tifc$Gmu == "6",]
    sf_tifc_gmu10a <- sf_tifc[sf_tifc$Gmu == "10A",]
    
    #'  Snag year data correspond to
    yr <- sf_tifc %>% dplyr::select(season) %>%
      slice(1L) %>%
      mutate(year = ifelse(season == "Smr20", "2020", season),
             year = ifelse(year == "Smr21", "2021", year),
             year = ifelse(year == "Smr22", "2022", year)) 
    yr <- as.data.frame(yr)
    yr <- yr$year
    
    #'  GMU 1 plot
    gmu1_tifc <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "1",], fill = NA) +
      geom_sf(data = sf_tifc_gmu1, aes(size = cpue_km2_sec), shape  = 21, 
              col = "darkred", fill = "darkred", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "CPUE km2/sec", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 6 plot
    gmu6_tifc <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "6",], fill = NA) +
      geom_sf(data = sf_tifc_gmu6, aes(size = cpue_km2_sec), shape  = 21, 
              col = "darkgreen", fill = "darkgreen", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "CPUE km2/sec", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 10A plot
    gmu10a_tifc <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
      geom_sf(data = sf_tifc_gmu10a, aes(size = cpue_km2_sec), shape = 21, 
              col = "darkblue", fill = "darkblue", alpha = 3/10) +
      scale_size_continuous(range = c(0,12)) +
      labs(size = "CPUE km2/sec", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  Plot each map
    print(gmu1_tifc); print(gmu6_tifc); print(gmu10a_tifc)
    
    #'  List GMU TIFC maps together
    gmu_maps <- list(gmu1_tifc, gmu6_tifc, gmu10a_tifc)
    
    return(gmu_maps)
  }
  tifc_maps_bear <- lapply(spatial_tifc_bear, map_tifc, spp = "Black bear")
  tifc_maps_bob <- lapply(spatial_tifc_bob, map_tifc, spp = "Bobcat")
  tifc_maps_coy <- lapply(spatial_tifc_coy, map_tifc, spp = "Coyote")
  tifc_maps_lion <- lapply(spatial_tifc_lion, map_tifc, spp = "Mountain lion")
  tifc_maps_wolf <- lapply(spatial_tifc_wolf, map_tifc, spp = "wolf")
  
  #'  Plot same species and study area back to back across years
  gmu_by_yr_plots <- function(fig, spp) {
    #'  Note: list order is [[i]][[j]] 
    #'  where i = 2020, 2021, or 2022 and j = GMU1, GMU6, or GMU10a 
    #'  GMU 1 plots
    gmu1_patch <- fig[[1]][[1]] + fig[[2]][[1]] + fig[[3]][[1]] +
      plot_annotation(paste("GMU 1", spp, "relative density index"))
    
    #'  GMU 6 plots
    gmu6_patch <- fig[[1]][[2]] + fig[[2]][[2]] + fig[[3]][[2]] +
      plot_annotation(paste("GMU 6", spp, "relative density index"))
    
    #'  GMU 10A plots
    gmu10a_patch <- fig[[1]][[3]] + fig[[2]][[3]] + fig[[3]][[3]] +
      plot_annotation(paste("GMU 10A", spp, "relative density index"))
    
    #'  Print figure panels
    print(gmu1_patch); print(gmu6_patch); print(gmu10a_patch)
    
    #'  List
    gmu_patchwork_list <- list(gmu1_patch, gmu6_patch, gmu10a_patch)
    
    return(gmu_patchwork_list)
  }
  tifc_gmu_bear <- gmu_by_yr_plots(tifc_maps_bear, spp = "black bear")
  tifc_gmu_bob <- gmu_by_yr_plots(tifc_maps_bob, spp = "bobcat")
  tifc_gmu_coy <- gmu_by_yr_plots(tifc_maps_coy, spp = "coyote")
  tifc_gmu_lion <- gmu_by_yr_plots(tifc_maps_lion, spp = "mountain lion")
  tifc_gmu_wolf <- gmu_by_yr_plots(tifc_maps_wolf, spp = "wolf")
  
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu1_bear.tiff", tifc_gmu_bear[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu6_bear.tiff", tifc_gmu_bear[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu10A_bear.tiff", tifc_gmu_bear[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu1_bob.tiff", tifc_gmu_bob[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu6_bob.tiff", tifc_gmu_bob[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu10A_bob.tiff", tifc_gmu_bob[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu1_coy.tiff", tifc_gmu_coy[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu6_coy.tiff", tifc_gmu_coy[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu10A_coy.tiff", tifc_gmu_coy[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu1_lion.tiff", tifc_gmu_lion[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu6_lion.tiff", tifc_gmu_lion[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu10A_lion.tiff", tifc_gmu_lion[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu1_wolf.tiff", tifc_gmu_wolf[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu6_wolf.tiff", tifc_gmu_wolf[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/TIFC/Figures/tifc_gmu10A_wolf.tiff", tifc_gmu_wolf[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  
  
  #'  Function to plot two species together in same year vs across years?
  
  
  
  #'  Next challenges...
  #'  1. Bootstrap estimates so we have standard errors for camera-level density
  #'  2. Do we trust EDD based on camera setups?
  #'  3. If yes, probably need EDD for prey too so effort is consistent across species
  #'  4. Which density metric to use (cpue, cpue km2, cpue 100km2)? -- going with cpue 100km2
  #'  5. Detection rate vs TIFC as response variable (both are going to be indices of abund/density)?
  #'  6. should I be using mean EDD or max EDD when calculating effort????
  
  