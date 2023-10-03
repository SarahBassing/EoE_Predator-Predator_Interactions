  #'  ----------------------------------------
  #'  Calculate denisty with TIFC method
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  September 2023
  #'  ----------------------------------------
  #'  Calculate density using time-in-front-of-camera methods following analyses
  #'  described by Huggard et al. 2018, Warbington & Boyce 2020, and Becker et al. 2022.
  #'  Code modified from Becker et al. 2018 repository: https://github.com/mabecker89/tifc-method/tree/main
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
  load("./Data/Relative abundance data/RAI Phase 2/eoe_all_fov-time_edd.RData")
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
      rename("detdist" = "avg_dist")
    return(df_tt_full)
  }
  df_tt_full_20s <- focal_spp(tt_list[[1]])
  df_tt_full_21s <- focal_spp(tt_list[[2]])
  df_tt_full_22s <- focal_spp(tt_list[[3]])
  
  #' #'  Load effective detection distance (EDD) 
  #' #'  Currently using detection zone measured at camera by IDFG. Area (m2) measured 
  #' #'  using walk test, flagging, or other methods to estimate detection zone
  #' #'  These values are not corrected for habitat type or species
  #' load("./Data/Relative abundance data/RAI Phase 2/area_m2.RData")
  #' area_m2_20s <- area_m2[[1]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  #' area_m2_21s <- area_m2[[2]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  #' area_m2_22s <- area_m2[[3]] %>% rename("location" = "NewLocationID") %>% rename("edd_area_m2" = "Area_M2")
  
  #'  Update species names for "distance groups" based on Becker et al. (2018) classifications
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
  
  #'  ---------------------
  ####  Calculate density  ####
  #'  ---------------------
  
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
      left_join(df_dist_groups, by = "common_name") #%>%
      # left_join(area_m2, by = "location") #%>%
    # left_join(df_detdist, by = c("location", "dist_group")) %>%
    # # Remove random species (mostly birds) 
    # filter(!is.na(detdist))
    return(df_dens_ing)
  }
  df_dens_ing_20s <- all_ing(df_tt_full_20s, df_dist_groups = df_dist_groups_20s, area_m2 = area_m2_20s)
  df_dens_ing_21s <- all_ing(df_tt_full_21s, df_dist_groups = df_dist_groups_21s, area_m2 = area_m2_21s)
  df_dens_ing_22s <- all_ing(df_tt_full_22s, df_dist_groups = df_dist_groups_22s, area_m2 = area_m2_22s)
  
  
  #####  Step 2. Calculate density  #####         # should I be using mean EDD or max EDD when calculating effort????
  #'  -----------------------------------
  calc_density <- function(df_dens_ing) {
    df_density <- df_dens_ing %>%
      mutate(effort = total_season_days * (detdist ^ 2 * pi * (cam_fov_angle / 360)) / 100,
      # mutate(effort = total_season_days * (edd_area_m2 * pi * (cam_fov_angle / 360)) / 100,
             #'  Catch per unit effort
             cpue = total_duration / effort,
             #'  Catch per unit effort in km2
             cpue_km2 = cpue / 60 / 60 / 24 * 10000,
             cpue_100km2 = cpue_km2 * 100) %>%
      return(df_density)
  }
  df_density_20s <- calc_density(df_dens_ing_20s); df_density_20s$season <- "Smr20"
  df_density_21s <- calc_density(df_dens_ing_21s); df_density_21s$season <- "Smr21"
  df_density_22s <- calc_density(df_dens_ing_22s); df_density_22s$season <- "Smr22"
  
  #'  Bind into a single data frame
  df_density <- rbind(df_density_20s, df_density_21s, df_density_22s) %>%
    #'  Add GMU to dataset
    mutate(gmu = sub("_.*", "", location),
           setup = ifelse(grepl("P", location), "predator", "ungulate")) %>%
    relocate(gmu, .after = "season") %>%
    relocate(setup, .after = "gmu")
  
  
  #'  SAVE! 
  write_csv(df_density, "./Data/Relative abundance data/RAI Phase 2/eoe_all_years_density_edd_predonly.csv")
  
  #'  Summary stats for density estimates
  density_stats <- df_density %>%
    group_by(gmu, common_name) %>%
    summarise(mean_density_km2 = round(mean(cpue_km2, na.rm = TRUE), 2),
              se_density_km2 = round((sd(cpue_km2, na.rm = TRUE)/sqrt(nrow(.))), 2),
              mean_density_100km2 = round(mean(cpue_100km2, na.rm = TRUE), 2),
              se_density_100km2 = round((sd(cpue_100km2, na.rm = TRUE)/sqrt(nrow(.))), 2)) %>%
    ungroup() %>%
    arrange(gmu, mean_density_km2)
  print(density_stats)
  
  write_csv(density_stats, "./Data/Relative abundance data/RAI Phase 2/tifc_density_stats_edd_predonly.csv")

  #'  Next challenges...
  #'  1. Bootstrap estimates so we have standard errors for camera-level density
  #'  2. Do we trust EDD based on camera setups?
  #'  3. If yes, probably need EDD for prey too so effort is consistent across species
  #'  4. Which density metric to use (cpue, cpue km2, cpue 100km2)?
  #'  5. How to remove investigation images for bears across entire dataset?
  #'  6. Detection rate vs TIFC as response variable (both are going to be indices of abund/density)?
  #'  ALSO- should I be using mean EDD or max EDD when calculating effort????
  
  