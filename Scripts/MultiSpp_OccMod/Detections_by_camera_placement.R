  #'  ----------------------------------
  #'  Spatial variation in detection
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  January 2023
  #'  ----------------------------------
  #'  Explore extent to which detections correlate with camera placement
  #'  ----------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(tidyverse)
  library(ggplot2)

  #'  Load detection histories
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20w_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")

  #'  Load camera station, problem cams, and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Target", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) #' Remove duplicate camera where height changed slightly
  
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  format_cam_station <- function(cams, season) {
    cams <- dplyr::select(cams, c("NewLocationID", "Lat", "Long")) %>%
      mutate(GMU = str_extract(NewLocationID, "[^_]+"),
             Season = season) %>%
      left_join(cams_eoe_long[cams_eoe_long$Season == season,], by = c("NewLocationID", "Season")) %>%
      #'  Drop handful of duplicated rows (not sure why this happens with left_join)
      unique() %>%
      #'  Reduce categories to random, road, trail
      mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | 
                                     CameraFacing == "decommision" | CameraFacing == "decommission", "road", CameraFacing),
             CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing),
             CameraFacing = ifelse(is.na(CameraFacing), "trail", CameraFacing)) %>% # predator cameras so either trail or road
      arrange(NewLocationID)
    return(cams)
  }
  cams_eoe20s <- format_cam_station(eoe_probcams_20s, season = "Smr20")
  cams_eoe20w <- format_cam_station(eoe_probcams_20w, season = "Wtr20")
  cams_eoe21s <- format_cam_station(eoe_probcams_21s, season = "Smr21") 
  
  #'  Load relative abundance indices
  load("./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Wtr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  
  #'  ----------------------------------------------------
  ####  Summary of detection events by camera deployment  ####
  #'  ----------------------------------------------------
  #'  Function to tally detections across primary period based on camera deployment
  #'  and plot distribution of detection/non-detections based on deployment type/location
  detection_by_placement <- function(dets, cams, sppseason) {
    #'  Sum total detections per camera and indicate if 1 or more detection occurred
    dets <- as.data.frame(dets) %>%
      rownames_to_column(var = "NewLocationID") %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(ndets = rowSums(.[2:11], na.rm = T),
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels = c("random", "road", "trail"))) %>%
      rowwise() %>%
      mutate(det = max(c_across(2:11), na.rm = T)) %>%
      dplyr::select(c(Setup, CameraFacing, det, ndets)) %>%
      filter(det != "-Inf")
    
    #'  Calculate percentage of data by camera setup and location
    print(sppseason)
    det_by_setup <- dplyr::select(dets, c(Setup, det)) %>%                                    
      group_by(Setup, det) %>%
      summarise(n_dets = table(det)) %>% ungroup() %>%
      group_by(det) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(det_by_setup)
    print(ggplot(dets, aes(x = det, fill = Setup)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      labs(fill="") +
      ggtitle(paste(sppseason, "detections by camera deployment")))
    
    det_bytrailtype <- dplyr::select(dets, c(CameraFacing, det)) %>%                                    
      group_by(CameraFacing, det) %>%
      summarise(n_dets = table(det)) %>% ungroup() %>%
      group_by(det) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(CameraFacing) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(det_bytrailtype)
    print(ggplot(dets, aes(x = det, fill = CameraFacing)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
      labs(fill="") +
      ggtitle(paste(sppseason, "detections by camera deployment")))
    
    ndet_by_setup <- dplyr::select(dets, c(Setup, ndets)) %>%                                    
      group_by(Setup, ndets) %>%
      summarise(n_dets = table(ndets)) %>% ungroup() %>%
      group_by(ndets) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(ndet_by_setup)
    print(ggplot(dets, aes(x = ndets, fill = Setup)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080")) +
      labs(fill="") +
      ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    ndet_bytrailtype <- dplyr::select(dets, c(CameraFacing, ndets)) %>%                                    
      group_by(CameraFacing, ndets) %>%
      summarise(n_dets = table(ndets)) %>% ungroup() %>%
      group_by(ndets) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(CameraFacing) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(ndet_bytrailtype)
    print(ggplot(dets, aes(x = ndets, fill = CameraFacing)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
      scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
      labs(fill="") +
      ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    return(dets)
  }
  wolf20s <- detection_by_placement(DH_eoe20s_predators[[5]][[1]], cams = cams_eoe20s, sppseason = "Wolf Smr20")
  wolf21s <- detection_by_placement(DH_eoe21s_predators[[5]][[1]], cams = cams_eoe21s, sppseason = "Wolf Smr21")
  coy20s <- detection_by_placement(DH_eoe20s_predators[[3]][[1]], cams = cams_eoe20s, sppseason = "Coyote Smr20")
  coy21s <- detection_by_placement(DH_eoe21s_predators[[3]][[1]], cams = cams_eoe21s, sppseason = "Coyote Smr21")
  lion20s <- detection_by_placement(DH_eoe20s_predators[[4]][[1]], cams = cams_eoe20s, sppseason = "Lion Smr20")
  lion21s <- detection_by_placement(DH_eoe21s_predators[[4]][[1]], cams = cams_eoe21s, sppseason = "Lion Smr21")
  bear20s <- detection_by_placement(DH_eoe20s_predators[[1]][[1]], cams = cams_eoe20s, sppseason = "Bear Smr20")
  bear21s <- detection_by_placement(DH_eoe21s_predators[[1]][[1]], cams = cams_eoe21s, sppseason = "Bear Smr21")
  bob20s <- detection_by_placement(DH_eoe20s_predators[[2]][[1]], cams = cams_eoe20s, sppseason = "Bobcat Smr20")
  bob21s <- detection_by_placement(DH_eoe21s_predators[[2]][[1]], cams = cams_eoe21s, sppseason = "Bobcat Smr21")
  
  
  #### table human detections, then make plots  ###
  load("./Data/Relative abundance data/EoE_RelativeN_30minElapsed_SamplingOcc.RData")

  human_detection_by_placement <- function(dets, cams, sppseason) {
    humans <- dets %>%
      filter(Species == "human") %>%
      group_by(NewLocationID) %>%
      mutate(ndets = sum(n_dets)) %>%
      slice(1L) %>%
      ungroup() %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(det = ifelse(is.na(ndets), 0, 1),
             ndets = ifelse(is.na(ndets), 0, ndets),
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels =c("random", "road", "trail"))) %>%
      dplyr::select(c(NewLocationID, Setup, CameraFacing, ndets, det)) 
    
    #'  Table data by camera set up and trail type
    print(sppseason)
    det_by_setup <- dplyr::select(humans, c(Setup, det)) %>%                                    
      group_by(Setup, det) %>%
      summarise(n_dets = table(det)) %>% ungroup() %>%
      group_by(det) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(det_by_setup)
    print(ggplot(humans, aes(x = det, fill = Setup)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
            scale_fill_manual(values=c("#69b3a2", "#404080")) +
            labs(fill="") +
            ggtitle(paste(sppseason, "detections by camera deployment")))
    
    det_bytrailtype <- dplyr::select(humans, c(CameraFacing, det)) %>%                                    
      group_by(CameraFacing, det) %>%
      summarise(n_dets = table(det)) %>% ungroup() %>%
      group_by(det) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(CameraFacing) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(det_bytrailtype)
    print(ggplot(humans, aes(x = det, fill = CameraFacing)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
            scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
            labs(fill="") +
            ggtitle(paste(sppseason, "detections by camera deployment")))
    
    ndet_by_setup <- dplyr::select(humans, c(Setup, ndets)) %>%                                    
      group_by(Setup, ndets) %>%
      summarise(n_dets = table(ndets)) %>% ungroup() %>%
      group_by(ndets) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(ndet_by_setup)
    print(ggplot(humans, aes(x = ndets, fill = Setup)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 5) +
            scale_fill_manual(values=c("#69b3a2", "#404080")) +
            labs(fill="") +
            ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    ndet_bytrailtype <- dplyr::select(humans, c(CameraFacing, ndets)) %>%                                    
      group_by(CameraFacing, ndets) %>%
      summarise(n_dets = table(ndets)) %>% ungroup() %>%
      group_by(ndets) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(CameraFacing) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(ndet_bytrailtype)
    print(ggplot(humans, aes(x = ndets, fill = CameraFacing)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 5) +
            scale_fill_manual(values=c("#69b3a2", "#404080", "#F8C8DC")) +
            labs(fill="") +
            ggtitle(paste("Sum of", sppseason, "detections by camera deployment")))
    
    return(humans)
  }
  humans20s <- human_detection_by_placement(eoe_30min_sampocc_list[[1]], cams = cams_eoe20s, sppseason = "Human Smr20")
  humans21s <- human_detection_by_placement(eoe_30min_sampocc_list[[3]], cams = cams_eoe21s, sppseason = "Human Smr21")
  
  
  #'  Table species pairings
  paired_detections <- function(detsA, detsB, cams, spp1, spp2) {
    #'  Count number of detections of spp 1 & indicate if ever detected at site
    spp1_dets <- as.data.frame(detsA) %>%
      rownames_to_column(var = "NewLocationID") %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(spp1 = spp1,
             ndetsA = rowSums(.[2:11], na.rm = T),
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels = c("random", "road", "trail"))) %>%
      rowwise() %>%
      mutate(detA = max(c_across(2:11), na.rm = T)) %>%
      dplyr::select(c(NewLocationID, Setup, CameraFacing, detA, ndetsA)) %>%
      filter(detA != "-Inf")
    #'  Count number of detections of spp2 & indicate if ever detected at site
    spp2_dets <- as.data.frame(detsB) %>%
      rownames_to_column(var = "NewLocationID") %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(spp1 = spp2,
             ndetsB = rowSums(.[2:11], na.rm = T),
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels = c("random", "road", "trail"))) %>%
      rowwise() %>%
      mutate(detB = max(c_across(2:11), na.rm = T)) %>%
      dplyr::select(c(NewLocationID, Setup, CameraFacing, detB, ndetsB)) %>%
      filter(detB != "-Inf")
    #'  Merge into single data frame & note if both species detected at same site
    spppair_dets <- full_join(spp1_dets, spp2_dets, by = c("NewLocationID", "Setup", "CameraFacing")) %>%
      mutate(detAB = sum(detA, detB, na.rm = T))
    #'  Create column with categories for each species, both, or none detected at a site
    det2spp_by_setup <- as.data.frame(spppair_dets) %>%
      dplyr::select(c(Setup, detA, detB, detAB)) %>%
      mutate(spp_dets = ifelse(detAB == 2, "both", detAB),
             spp_dets = ifelse(detA == 1 & detB == 0, spp1, spp_dets),
             spp_dets = ifelse(detA == 0 & detB == 1, spp2, spp_dets),
             spp_dets = ifelse(detAB == 0, "none", spp_dets)) %>%
      dplyr::select(Setup, spp_dets) %>%
      group_by(Setup, spp_dets) %>%
      #'  Count number of detections in each grouping
      summarize(n = n()) %>%
      #'  Ignore sites where neither species were detected
      filter(spp_dets != "none") %>%
      #'  Calculate percentage of detections of individual species and both species
      #'  across sites where detections occurred
      mutate(prec_dets = n / sum(n),
             prec_dets = round(prec_dets, 3)) %>%
      ungroup()
    print(det2spp_by_setup)
  }
  wolfbear_20s <- paired_detections(detsA = DH_eoe20s_predators[[5]][[1]], detsB = DH_eoe20s_predators[[1]][[1]], 
                                    cams = cams_eoe20s, spp1 = "wolf", spp2 = "bear")
  wolfbear_21s <- paired_detections(detsA = DH_eoe21s_predators[[5]][[1]], detsB = DH_eoe21s_predators[[1]][[1]], 
                                    cams = cams_eoe21s, spp1 = "wolf", spp2 = "bear")
  wolfbob_20s <- paired_detections(detsA = DH_eoe20s_predators[[5]][[1]], detsB = DH_eoe20s_predators[[2]][[1]], 
                                   cams = cams_eoe20s, spp1 = "wolf", spp2 = "bobcat")
  wolfbob_21s <- paired_detections(detsA = DH_eoe21s_predators[[5]][[1]], detsB = DH_eoe21s_predators[[2]][[1]], 
                                   cams = cams_eoe21s, spp1 = "wolf", spp2 = "bobcat")
  wolfcoy_20s <- paired_detections(detsA = DH_eoe20s_predators[[5]][[1]], detsB = DH_eoe20s_predators[[3]][[1]], 
                                   cams = cams_eoe20s, spp1 = "wolf", spp2 = "coyote")
  wolfcoy_21s <- paired_detections(detsA = DH_eoe21s_predators[[5]][[1]], detsB = DH_eoe21s_predators[[3]][[1]], 
                                   cams = cams_eoe21s, spp1 = "wolf", spp2 = "coyote")
  wolflion_20s <- paired_detections(detsA = DH_eoe20s_predators[[5]][[1]], detsB = DH_eoe20s_predators[[4]][[1]], 
                                    cams = cams_eoe20s, spp1 = "wolf", spp2 = "lion")
  wolflion_21s <- paired_detections(detsA = DH_eoe21s_predators[[5]][[1]], detsB = DH_eoe21s_predators[[4]][[1]], 
                                    cams = cams_eoe21s, spp1 = "wolf", spp2 = "lion")
  lionbear_20s <- paired_detections(detsA = DH_eoe20s_predators[[4]][[1]], detsB = DH_eoe20s_predators[[1]][[1]], 
                                    cams = cams_eoe20s, spp1 = "lion", spp2 = "bear")
  lionbear_21s <- paired_detections(detsA = DH_eoe21s_predators[[4]][[1]], detsB = DH_eoe21s_predators[[1]][[1]], 
                                    cams = cams_eoe21s, spp1 = "lion", spp2 = "bear")
  lionbob_20s <- paired_detections(detsA = DH_eoe20s_predators[[4]][[1]], detsB = DH_eoe20s_predators[[2]][[1]], 
                                   cams = cams_eoe20s, spp1 = "lion", spp2 = "bobcat")
  lionbob_21s <- paired_detections(detsA = DH_eoe21s_predators[[4]][[1]], detsB = DH_eoe21s_predators[[2]][[1]], 
                                   cams = cams_eoe21s, spp1 = "lion", spp2 = "bobcat")
  lioncoy_20s <- paired_detections(detsA = DH_eoe20s_predators[[4]][[1]], detsB = DH_eoe20s_predators[[3]][[1]], 
                                    cams = cams_eoe20s, spp1 = "lion", spp2 = "coyote")
  lioncoy_21s <- paired_detections(detsA = DH_eoe21s_predators[[4]][[1]], detsB = DH_eoe21s_predators[[3]][[1]], 
                                    cams = cams_eoe21s, spp1 = "lion", spp2 = "coyote")
  bearbob_20s <- paired_detections(detsA = DH_eoe20s_predators[[1]][[1]], detsB = DH_eoe20s_predators[[2]][[1]], 
                                    cams = cams_eoe20s, spp1 = "bear", spp2 = "bobcat")
  bearbob_21s <- paired_detections(detsA = DH_eoe21s_predators[[1]][[1]], detsB = DH_eoe21s_predators[[2]][[1]], 
                                    cams = cams_eoe21s, spp1 = "bear", spp2 = "bobcat")
  bearcoy_20s <- paired_detections(detsA = DH_eoe20s_predators[[1]][[1]], detsB = DH_eoe20s_predators[[3]][[1]], 
                                   cams = cams_eoe20s, spp1 = "bear", spp2 = "coyote")
  bearcoy_21s <- paired_detections(detsA = DH_eoe21s_predators[[1]][[1]], detsB = DH_eoe21s_predators[[3]][[1]], 
                                   cams = cams_eoe21s, spp1 = "bear", spp2 = "coyote")
  bobcoy_20s <- paired_detections(detsA = DH_eoe20s_predators[[2]][[1]], detsB = DH_eoe20s_predators[[3]][[1]], 
                                   cams = cams_eoe20s, spp1 = "bobcat", spp2 = "coyote")
  bobcoy_21s <- paired_detections(detsA = DH_eoe21s_predators[[2]][[1]], detsB = DH_eoe21s_predators[[3]][[1]], 
                                   cams = cams_eoe21s, spp1 = "bobcat", spp2 = "coyote")
    
  
  ####  Relative abundance by camera setup  ####
  prey_detection_by_placement <- function(ra, spp, cams, sppseason) {
    #'  Select necessary covariate data
    ra <- dplyr::select(ra, c("NewLocationID", spp))
    #'  Sum relative abundance per camera and indicate if 1 or more detection occurred
    relative_n <- as.data.frame(ra) %>%
      full_join(cams, by = "NewLocationID") %>%
      mutate(n = .[,2],
             Setup = factor(Setup, levels = c("ungulate", "predator")),
             CameraFacing = factor(CameraFacing, levels = c("random", "road", "trail"))) %>%
      rowwise() %>%
      mutate(det = ifelse(n > 0, 1, 0)) %>%
      dplyr::select(c(Setup, CameraFacing, det, n)) 
    
    #'  Calculate percentage of data by camera setup
    #'  Binary detection of prey at any point in time per season
    print(sppseason)
    det_by_setup <- dplyr::select(relative_n, c(Setup, det)) %>%                                    
      group_by(Setup, det) %>%
      summarise(n_dets = table(det)) %>% ungroup() %>%
      group_by(det) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(det_by_setup)
    print(ggplot(relative_n, aes(x = det, fill = Setup)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
            scale_fill_manual(values=c("#69b3a2", "#404080")) +
            labs(fill="") +
            ggtitle(paste(sppseason, "detections by camera deployment")))
    
    #'  Relative abundance of prey at each camera
    n_by_setup <- dplyr::select(relative_n, c(Setup, n)) %>%                                    
      group_by(Setup, n) %>%
      summarise(n_dets = table(n)) %>% ungroup() %>%
      group_by(n) %>%
      mutate(`percent of 0s or 1s` = n_dets/sum(n_dets)) %>% ungroup() %>%
      group_by(Setup) %>%
      mutate(`percent of U vs P dets` = n_dets/sum(n_dets)) %>% ungroup()
    print(n_by_setup)
    print(ggplot(relative_n, aes(x = n, fill = Setup)) +
            geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 2) +
            scale_fill_manual(values=c("#69b3a2", "#404080")) +
            labs(fill="") +
            ggtitle(paste(sppseason, "relative abundance by camera deployment")))
  }
  elk_20s <- prey_detection_by_placement(eoe_covs_20s, spp = "elk", cams = cams_eoe20s, sppseason = "Elk Summer 2020")
  elk_21s <- prey_detection_by_placement(eoe_covs_21s, spp = "elk", cams = cams_eoe21s, sppseason = "Elk Summer 2021")
  moose_20s <- prey_detection_by_placement(eoe_covs_20s, spp = "moose", cams = cams_eoe20s, sppseason = "Moose Summer 2020")
  moose_21s <- prey_detection_by_placement(eoe_covs_21s, spp = "moose", cams = cams_eoe21s, sppseason = "Moose Summer 2021")
  md_20s <- prey_detection_by_placement(eoe_covs_20s, spp = "muledeer", cams = cams_eoe20s, sppseason = "Mule Deer Summer 2020")
  md_21s <- prey_detection_by_placement(eoe_covs_21s, spp = "muledeer", cams = cams_eoe21s, sppseason = "Mule Deer Summer 2021")
  wtd_20s <- prey_detection_by_placement(eoe_covs_20s, spp = "whitetaileddeer", cams = cams_eoe20s, sppseason = "White-tailed Deer Summer 2020")
  wtd_21s <- prey_detection_by_placement(eoe_covs_21s, spp = "whitetaileddeer", cams = cams_eoe21s, sppseason = "White-tailed Deer Summer 2021")
  bunnies_20s <- prey_detection_by_placement(eoe_covs_20s, spp = "lagomorphs", cams = cams_eoe20s, sppseason = "Lagomorphs Summer 2020")
  bunnies_21s <- prey_detection_by_placement(eoe_covs_21s, spp = "lagomorphs", cams = cams_eoe21s, sppseason = "Lagomorphs Summer 2021")
  
  
  