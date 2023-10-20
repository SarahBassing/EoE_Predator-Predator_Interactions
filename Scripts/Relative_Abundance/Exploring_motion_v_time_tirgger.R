  #'  ----------------------------------
  #'  Time vs Motion trigger detections
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  October 2023
  #'  ----------------------------------
  #'  Compare detections of common and rare predators from time vs motion trigger
  #'  images to assess how much detection probability affects index of density
  #'  for each species relative to other species. Are our relative density estimates
  #'  for each predator proportionally consistent across image times or are do 
  #'  motion trigger images and imperfect detection have the potential to affect
  #'  relative relationshps between species?
  #'  ----------------------------------
  
  #'  Load libraries
  library(readr)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  -------------------------------
  ####  Load and format input data  ####
  #'  -------------------------------
  #'  Load motion trigger images
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_2023-09-26.RData") 
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
  #'  Filter motion trigger images to time period of interest and remove problem images
  thin_detections <- function(dets, seqprobs, start_date, end_date) {
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    
    clean_dets <- skinny_dets %>%
      filter(Species != "none") %>%
      filter(Vehicle != "TRUE") %>%
      filter(OpState != "maintenance") %>%
      #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
      mutate(Count = ifelse(Count == 0, 1, Count),
             Count = ifelse(is.na(Count), 1, Count),
             Date = as.Date(Date, format = "%Y-%m-%d")) %>% #format = "%d-%b-%Y"
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      arrange(NewLocationID, posix_date_time) %>%
      dplyr::select(c("NewLocationID", "Setup", "posix_date_time", "Species", "Count"))

    return(clean_dets)
  }
  mt_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15") 
  mt_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15")
  mt_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-07-01", end_date = "2022-09-15")
  
  #' #'  Load time trigger images
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allT_NewLocationID.RData")
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allT_NewLocationID.RData")
  #' load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allT_NewLocationID.RData")
  #' 
  #' #'  Filter time trigger images to time period of interest and remove empty images
  #' thin_detections <- function(dets, start_date, end_date) {
  #'   clean_dets <- dets %>%
  #'     filter(Species != "none") %>%
  #'     filter(Vehicle != "TRUE") %>%
  #'     filter(OpState != "maintenance") %>%
  #'     #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
  #'     mutate(Count = ifelse(Count == 0, 1, Count),
  #'            Count = ifelse(is.na(Count), 1, Count),
  #'            Date = as.Date(Date, format = "%Y-%m-%d")) %>% #format = "%d-%b-%Y"
  #'     #'  Filter to images to desired date range
  #'     filter(Date >= start_date & Date <= end_date) %>%
  #'     #'  Remove observations that can't be linked to a camera with coordinates
  #'     filter(!is.na(NewLocationID)) %>%
  #'     arrange(NewLocationID, posix_date_time) %>%
  #'     dplyr::select(c("NewLocationID", "Setup", "posix_date_time", "Species", "Count"))
  #' 
  #'   return(clean_dets)
  #' }
  #' tt_20s <- thin_detections(eoe20s_allT, start_date = "2020-07-01", end_date = "2020-09-15")
  #' tt_21s <- thin_detections(eoe21s_allT, start_date = "2021-07-01", end_date = "2021-09-15")
  #' tt_22s <- thin_detections(eoe22s_allT, start_date = "2022-07-01", end_date = "2022-09-15")
  #' 
  #' save(tt_20s, file = "./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_20s.RData")
  #' save(tt_21s, file = "./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_21s.RData")
  #' save(tt_22s, file = "./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_22s.RData")
  
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_20s.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_21s.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/tt_22s.RData")
  
  #'  Sampling effort data (number of days cameras operation)
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData") 
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe22s.RData") 
  
  #'  Calculate total number of images per camera for each species and trigger type
  n_pix <- function(dets, days_operable) {
    #'  Sum all images pre species at each camera across entire season
    all_imgs <- dets %>%
      group_by(NewLocationID, Setup, Species) %>%
      summarize(n_imgs = n()) %>%
      ungroup()
    
    #'  Add sampling effort to data frame
    effort <- dplyr::select(days_operable, c(NewLocationID, ndays))
    dets <- left_join(all_imgs, effort, by = "NewLocationID") %>%
      #'  Retain data from cameras where camera was operable for 30 days or more
      filter(ndays >= 30)
    return(dets)
  }
  nimgs_20s_m <- n_pix(mt_20s, days_operable = effort_20s)
  nimgs_21s_m <- n_pix(mt_21s, days_operable = effort_21s)
  nimgs_22s_m <- n_pix(mt_22s, days_operable = effort_22s)
  nimgs_20s_t <- n_pix(tt_20s, days_operable = effort_20s)
  nimgs_21s_t <- n_pix(tt_21s, days_operable = effort_21s)
  nimgs_22s_t <- n_pix(tt_22s, days_operable = effort_22s)
  
  
  #'  Compare detection counts between different trigger types
  compare_counts <- function(mt, tt, season, not_setup, setup) {
    #'  Bind tallied detections from motion and time trigger images
    all_triggers <-  left_join(mt, tt, by = c("NewLocationID", "Setup", "Species", "ndays")) %>%
      dplyr::select(-ndays) %>%
      #'  Filter to camera setup of interest (using != so I can feed NULL to it if desired)
      filter(Setup != not_setup) %>%
      #'  Reduce to species of interest and remove sites with all NAs
      # filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
      #          Species == "elk" | Species == "human" | Species == "moose" |
      #          Species == "mountain_lion" | Species == "muledeer" | Species == "rabbit_hare" |
      #          Species == "whitetaileddeer" | Species == "wolf" | Species == "cattle_cow") %>%
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "mountain_lion" | Species == "wolf") %>%
      #'  Replace NAs with 0 since one of the two trigger types failed to detect 
      #'  that species when the other trigger type did
      mutate(n_imgs.x = ifelse(is.na(n_imgs.x), 0, n_imgs.x),
             n_imgs.y = ifelse(is.na(n_imgs.y), 0, n_imgs.y),
             Season = season) %>%
      rename("n_imgs_mt" = "n_imgs.x") %>%
      rename("n_imgs_tt" = "n_imgs.y") %>%
      relocate(Season, .after = "NewLocationID")
    
    #'  Sum number of images per species and trigger type
    sum_imgs <- all_triggers %>%
      group_by(Species) %>%
      summarise(n_mt = sum(n_imgs_mt),
                n_tt = sum(n_imgs_tt)) %>%
      ungroup() 
    
    #'  Calculate proportion of images for each species relative to black bear 
    #'  images for each trigger type
    bear_mt <- sum_imgs$n_mt[sum_imgs$Species == "bear_black"]
    bear_tt <- sum_imgs$n_tt[sum_imgs$Species == "bear_black"]
    bob_mt <- sum_imgs$n_mt[sum_imgs$Species == "bobcat"]
    bob_tt <- sum_imgs$n_tt[sum_imgs$Species == "bobcat"]
    proportion_imgs <- sum_imgs %>%
      mutate(relative2bear_mt = round(n_mt/bear_mt, 3),
             relative2bear_tt = round(n_tt/bear_tt, 3),
             relative2bob_mt = round(n_mt/bob_mt, 3),
             relative2bob_tt = round(n_tt/bob_tt, 3),
             #'  If sum of bobcat images = 0 then can't calculate a proportion so change to NA
             relative2bob_tt = ifelse(is.infinite(relative2bob_tt), NA, relative2bob_tt),
             relative2bob_tt = ifelse(relative2bob_tt == "NaN", NA, relative2bob_tt),
             Season = season,
             Setup = setup) %>%
      relocate(Season, .after = "Species") %>%
      relocate(Setup, .after = "Season")
    print(proportion_imgs)
    
    return_dat <- list(all_triggers, proportion_imgs)
    
    return(return_dat)
  }
  mt_tt_20s_ung <- compare_counts(nimgs_20s_m, nimgs_20s_t, season = "Smr20", not_setup = "P", setup = "U")
  mt_tt_20s_pred <- compare_counts(nimgs_20s_m, nimgs_20s_t, season = "Smr20", not_setup = "U", setup = "P") 
  # mt_tt_20s <- compare_counts(nimgs_20s_m, nimgs_20s_t, season = "Smr20", not_setup = "NULL", setup = "Both") 
  
  mt_tt_21s_ung <- compare_counts(nimgs_21s_m, nimgs_21s_t, season = "Smr21", not_setup = "P", setup = "U") 
  mt_tt_21s_pred <- compare_counts(nimgs_21s_m, nimgs_21s_t, season = "Smr21", not_setup = "U", setup = "P") 
  # mt_tt_21s <- compare_counts(nimgs_21s_m, nimgs_21s_t, season = "Smr21", not_setup = "NULL", setup = "Both") 
  
  mt_tt_22s_ung <- compare_counts(nimgs_22s_m, nimgs_22s_t, season = "Smr22", not_setup = "P", setup = "U") 
  mt_tt_22s_pred <- compare_counts(nimgs_22s_m, nimgs_22s_t, season = "Smr22", not_setup = "U", setup = "P") 
  # mt_tt_22s <- compare_counts(nimgs_22s_m, nimgs_22s_t, season = "Smr22", not_setup = "NULL", setup = "Both") 
  
  #'  Bind summary info together
  mt_tt_summary <- rbind(mt_tt_20s_pred[[2]], mt_tt_20s_ung[[2]], mt_tt_21s_pred[[2]], 
                         mt_tt_21s_ung[[2]], mt_tt_22s_pred[[2]], mt_tt_22s_ung[[2]])
    
    
  #'  Compare relationships between predators for MT, TT, and TIFC data
  #'  Read in TIFC data
  tifc <- read_csv("./Data/Relative abundance data/RAI Phase 2/eoe_all_years_density_avg_edd_predonly.csv")
  
  #'  Summary stats
  density_stats <- tifc %>%
    group_by(common_name, setup, season) %>%
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
  density_by_setup <- full_join(predcam_density, ungcam_density, by = c("common_name", "season"))
  print(density_by_setup)
  
  
  #'  Setup data sets for regression analyses
  mt_tt_nimgs <- rbind(mt_tt_20s_pred[[1]], mt_tt_20s_ung[[1]], mt_tt_21s_pred[[1]], 
                       mt_tt_21s_ung[[1]], mt_tt_22s_pred[[1]], mt_tt_22s_ung[[1]]) 
  mt_nimgs <- mt_tt_nimgs %>%
    dplyr::select(-n_imgs_tt) %>%
    mutate(obs = seq(1:nrow(.))) %>%
    #'  Make sure each site has an observation for each species - fill in with 0 if not detected
    spread(Species, n_imgs_mt, fill = 0) %>%
    dplyr::select(-obs) %>%
    #'  Convert back to long format
    gather(Species, n_imgs_mt, -c(NewLocationID, Season, Setup)) %>%
    group_by(NewLocationID, Species, Season, Setup) %>%
    #'  Drop all the extra observations
    slice_max(n_imgs_mt) %>%
    ungroup() %>%
    #'  Remove duplicate rows where n_imgs = 0 for a given species and site
    distinct() %>%
    arrange(NewLocationID, Season, Setup, Species) %>%
    rename("n_imgs" = "n_imgs_mt")
  
  #'  Same but with time trigger images
  tt_nimgs <- mt_tt_nimgs %>%
    dplyr::select(-n_imgs_mt) %>%
    mutate(obs = seq(1:nrow(.))) %>%
    spread(Species, n_imgs_tt, fill = 0) %>%
    dplyr::select(-obs) %>%
    gather(Species, n_imgs_tt, -c(NewLocationID, Season, Setup)) %>%
    group_by(NewLocationID, Species, Season, Setup) %>%
    slice_max(n_imgs_tt) %>%
    ungroup() %>%
    distinct() %>%
    arrange(NewLocationID, Season, Setup, Species) %>%
    rename("n_imgs" = "n_imgs_tt")
  
  #'  Regressions to assess if there are at least consistent relationships btwn
  #'  bobcat and wolf metrics, even if the counts/densities relative to one another
  #'  differ depending on camera set up, trigger type, and count metric
  #'  Ideally tt and tifc relationship are similar...
  mt_wolf_bob <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "bobcat"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_bob) # (+) NS effect
  tt_wolf_bob <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "bobcat"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_bob) # (-) NS effect
  tifc_wolf_bob <- lm(tifc$cpue_km2_sec[tifc$common_name == "bobcat"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"])
  summary(tifc_wolf_bob) # (+) NS effect
  
  mt_wolf_bear <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "bear_black"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_bear) # (+) NS effect
  tt_wolf_bear <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "bear_black"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_bear) # (-) NS effect
  tifc_wolf_bear <- lm(tifc$cpue_km2_sec[tifc$common_name == "bear_black"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"])
  summary(tifc_wolf_bear) # (+) NS effect
  
  mt_wolf_lion <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "mountain_lion"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_lion) # (+) NS effect
  tt_wolf_lion <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "mountain_lion"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_lion) # (-) NS effect
  tifc_wolf_lion <- lm(tifc$cpue_km2_sec[tifc$common_name == "mountain_lion"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"])
  summary(tifc_wolf_lion) # (+) NS effect
  
  mt_wolf_lion_setup <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "mountain_lion"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"] + mt_nimgs$Setup[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_lion_setup) # wolf (1) NS effect, U cam (-) SIGNIFICANT effect
  tt_wolf_lion_setup <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "mountain_lion"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"] + tt_nimgs$Setup[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_lion_setup) # wolf (`+`) NS effect, U cam (-) significant effect
  tifc_wolf_lion_setup <- lm(tifc$cpue_km2_sec[tifc$common_name == "mountain_lion"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"] + tifc$setup[tifc$common_name == "wolf"])
  summary(tifc_wolf_lion_setup) # wolf (1) NS effect, U cam (-) SIGNIFICANT effect
  
  mt_wolf_coy <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "coyote"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_coy) # (+) Signif effect
  tt_wolf_coy <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "coyote"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_coy) # (+) NS effect
  tifc_wolf_coy <- lm(tifc$cpue_km2_sec[tifc$common_name == "coyote"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"])
  summary(tifc_wolf_coy) # (+) Signif effect
  
  mt_wolf_coy_setup <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "coyote"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"] + mt_nimgs$Setup[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_coy_setup) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  tt_wolf_coy_setup <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "coyote"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"] + tt_nimgs$Setup[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_coy_setup) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  tifc_wolf_coy_setup <- lm(tifc$cpue_km2_sec[tifc$common_name == "coyote"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"] + tifc$setup[tifc$common_name == "wolf"])
  summary(tifc_wolf_coy_setup) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  
  mt_wolf_coy_setup_season <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "coyote"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "wolf"] + mt_nimgs$Setup[mt_nimgs$Species == "wolf"] + mt_nimgs$Season[mt_nimgs$Species == "wolf"])
  summary(mt_wolf_coy_setup_season) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  tt_wolf_coy_setup_season <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "coyote"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "wolf"] + tt_nimgs$Setup[tt_nimgs$Species == "wolf"] + tt_nimgs$Season[tt_nimgs$Species == "wolf"])
  summary(tt_wolf_coy_setup_season) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  tifc_wolf_coy_setup_season <- lm(tifc$cpue_km2_sec[tifc$common_name == "coyote"] ~ tifc$cpue_km2_sec[tifc$common_name == "wolf"] + tifc$setup[tifc$common_name == "wolf"] + tifc$season[tifc$common_name == "wolf"])
  summary(tifc_wolf_coy_setup_season) # wolf (+) NS effect, U cam (-) SIGNIFICANT effect
  
  mt_coy_bob <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "bobcat"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "coyote"])
  summary(mt_coy_bob) # (+) Signif effect
  tt_coy_bob <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "bobcat"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "coyote"])
  summary(tt_coy_bob) # (+) Signif effect
  tifc_coy_bob <- lm(tifc$cpue_km2_sec[tifc$common_name == "bobcat"] ~ tifc$cpue_km2_sec[tifc$common_name == "coyote"])
  summary(tifc_coy_bob) # (+) Signif effect
  
  mt_coy_bob_setup <- lm(mt_nimgs$n_imgs[mt_nimgs$Species == "bobcat"] ~ mt_nimgs$n_imgs[mt_nimgs$Species == "coyote"] + mt_nimgs$Setup[mt_nimgs$Species == "coyote"])
  summary(mt_coy_bob_setup) # coy (+) Signif effect, U cam (-) SIGNIFICANT effect
  tt_coy_bob_setup <- lm(tt_nimgs$n_imgs[tt_nimgs$Species == "bobcat"] ~ tt_nimgs$n_imgs[tt_nimgs$Species == "coyote"] + tt_nimgs$Setup[tt_nimgs$Species == "coyote"])
  summary(tt_coy_bob_setup) # coy (+) NS effect, U cam (-) SIGNIFICANT effect
  tifc_coy_bob_setup <- lm(tifc$cpue_km2_sec[tifc$common_name == "bobcat"] ~ tifc$cpue_km2_sec[tifc$common_name == "coyote"] + tifc$setup[tifc$common_name == "coyote"])
  summary(tifc_coy_bob_setup) # coy (+) Signif effect, U cam (-) SIGNIFICANT effect
  
  
  
  