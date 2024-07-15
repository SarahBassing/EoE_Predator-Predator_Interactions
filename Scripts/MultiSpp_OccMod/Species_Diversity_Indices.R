  #'  ------------------------------
  #'  Diversity indices
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ------------------------------
  #'  Code to calculate Species Richness (S) and Shannon's Diversity Index (H) 
  #'  for prey species detected at each camera site. Prey species include wild
  #'  (elk, moose, mule deer, white-tailed deer) and domestic (cattle) ungulates.
  #'  -----------------------------
  
  #'  Load packages
  library(tidyverse)
  
  #'  Relative abundance index data for prey species
  #'  Using Hour of Detection as RA index b/c highly correlated with other 
  #'  definitions of independent detection events and consistent with Ausband et al. 2023
  load("./Data/Relative abundance data/EoE_RelativeN_30minElapsed.RData") 
  load("./Data/Relative abundance data/EoE_RelativeN_HrOfDetection.RData")
  
  #'  Sampling effort data (number of days cameras operation)
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe20w.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/SamplingEffort_eoe21s.RData")
  
  #'  Format relative abundance indices
  reformat_relativeN_data <- function(dat) {
    RelativeN <- dat %>%
      #'  Create one column per species
      spread(Species, n_dets) %>%
      rowwise() %>%
      #'  Aggrigate relative abundance data into functional groups
      mutate(human_plus = sum(human, cat_domestic, dog_domestic, horse, na.rm = TRUE),
             ungulate = sum(elk, moose, muledeer, whitetaileddeer, na.rm = TRUE),
             big_deer = sum(elk, moose, na.rm = TRUE),
             small_deer = sum(muledeer, whitetaileddeer, na.rm = TRUE)) %>%
      #'  Rename and drop columns
      rename(livestock = cattle_cow) %>%
      rename(lagomorphs = rabbit_hare) %>%
      dplyr::select(-c(cat_domestic, dog_domestic, horse)) %>%
      #'  Change all NAs introduced during spread to 0's (i.e., non-detection)
      replace(is.na(.), 0) %>%
      #'  Relocate columns so more intuitive order
      relocate(human_plus, .after = human) %>%
      relocate(lagomorphs, .before = moose) %>%
      relocate(livestock, .before = moose) %>%
      relocate(ungulate, .after = whitetaileddeer) %>%
      relocate(big_deer, .after = ungulate) %>%
      relocate(small_deer, .after = big_deer)
    return(RelativeN)
  }
  # RA_Smr20_df <- reformat_relativeN_data(eoe_30min_list[[1]]) 
  # RA_Wtr20_df <- reformat_relativeN_data(eoe_30min_list[[2]])
  # RA_Smr21_df <- reformat_relativeN_data(eoe_30min_list[[3]])
  RA_Smr20_df <- reformat_relativeN_data(eoe_dethr_list[[1]]) 
  RA_Wtr20_df <- reformat_relativeN_data(eoe_dethr_list[[2]])
  RA_Smr21_df <- reformat_relativeN_data(eoe_dethr_list[[3]])
  
  #'  Adjust relative abundance indices by sampling effort to create a DAILY DETECTION RATE
  #'  Divide total number of hours with at least one detection (summed across season)
  #'  by the number of days camera was operating, i.e., the average number of hours
  #'  per day at least one animal was detected
  weighted_RA <- function(RA, effort) {
    effort <- dplyr::select(effort, c("NewLocationID", "ndays", "nhrs"))
    ra_scaled_by_nhrs <- RA %>%
      full_join(effort, by = "NewLocationID") %>%
      # mutate(elk_perday = round(elk/ndays, 3) * 100,
      #        human_perday = round(human/ndays, 3) * 100,
      #        human_perday = round(human_plus/ndays, 3) * 100,
      #        human_motorized_perday = round(human_motorized/ndays, 3) * 100,
      #        lagomorphs_perday = round(lagomorphs/ndays, 3) * 100,
      #        livestock_perday = round(livestock/ndays, 3) * 100,
      #        moose_perday = round(moose/ndays, 3) * 100,
      #        muledeer_perday = round(muledeer/ndays, 3) * 100,
      #        whitetaileddeer_perday = round(whitetaileddeer/ndays, 3) * 100,
      #        ungulate_perday = round(ungulate/ndays, 3) * 100,
      #        big_deer_perday = round(big_deer/ndays, 3) * 100,
      #        small_deer_perday = round(small_deer/ndays, 3) * 100) %>%
      mutate(elk_perday = round(elk/ndays, 3),
             human_perday = round(human/ndays, 3),
             human_perday = round(human_plus/ndays, 3),
             human_motorized_perday = round(human_motorized/ndays, 3),
             lagomorphs_perday = round(lagomorphs/ndays, 3),
             livestock_perday = round(livestock/ndays, 3),
             moose_perday = round(moose/ndays, 3),
             muledeer_perday = round(muledeer/ndays, 3),
             whitetaileddeer_perday = round(whitetaileddeer/ndays, 3),
             ungulate_perday = round(ungulate/ndays, 3),
             big_deer_perday = round(big_deer/ndays, 3),
             small_deer_perday = round(small_deer/ndays, 3)) %>%
      dplyr::select(-c("ndays", "nhrs"))
    return(ra_scaled_by_nhrs)
  }
  RA_Smr20_df <- weighted_RA(RA_Smr20_df, effort = effort_20s)
  RA_Wtr20_df <- weighted_RA(RA_Wtr20_df, effort = effort_20w)
  RA_Smr21_df <- weighted_RA(RA_Smr21_df, effort = effort_21s)
  
  
  #'  ---------------------------------------
  ####  Calculate Species Diversity Indices  ####
  #'  ---------------------------------------
  #'  Function to calculate species richness, Shannon's diversity index, and flag
  #'  the most frequently detected species at each camera
  species_diversity <- function(RA) {
    #'  Species Richness (S)
    #'  --------------------
    #'  Sum number of unique species detected at each camera
    SR <- RA %>% 
      dplyr::select(c("NewLocationID", "elk", "lagomorphs", "livestock", "moose", "muledeer", "whitetaileddeer")) %>%
      mutate(elk_det = ifelse(elk > 0, 1, 0),
             lagomorph_det = ifelse(lagomorphs > 0, 1, 0),
             livestock_det = ifelse(livestock > 0, 1, 0),
             moose_det = ifelse(moose > 0, 1, 0),
             muledeer_det = ifelse(muledeer > 0, 1, 0),
             whitetaileddeer_det = ifelse(whitetaileddeer > 0, 1, 0),
             SR = sum(elk_det, lagomorph_det, livestock_det, moose_det, muledeer_det, whitetaileddeer_det))
    
    #'  Shannon's diversity index (H)
    #'  -----------------------------
    #'  Considers species richness and evenness (abundance of each species)
    #'  https://www.programmingr.com/shannon-diversity-index-the-diversity-function-in-r/
    # Shannon <- as.data.frame(RA) %>% 
    #   dplyr::select(c("NewLocationID", "elk_perday", "lagomorphs_perday", "livestock_perday",
    #                   "moose_perday", "muledeer_perday", "whitetaileddeer_perday")) %>%
    #   filter(!is.na(elk_perday))
    #'  Alternatively, use un-weighted RA index
    #'  FYI: H values are almost identical to those of weighted RA index
    Shannon <- as.data.frame(RA) %>%
      dplyr::select(c("NewLocationID", "elk", "lagomorphs", "livestock", "moose", "muledeer", "whitetaileddeer")) %>%
      filter(!is.na(elk))
    
    Shannon_noLago <- as.data.frame(RA) %>%
      dplyr::select(c("NewLocationID", "elk", "livestock", "moose", "muledeer", "whitetaileddeer")) %>%
      filter(!is.na(elk))
  
    #'  Loop through each camera site to calculate H
    H <- c(NA)
    for(i in 1:nrow(Shannon)) {
      #'  Relative abundance of each species
      n <- c(Shannon[i,2], Shannon[i,3], Shannon[i,4], Shannon[i,5], Shannon[i,6], Shannon[i,7]) 
      #'  Remove species that were not detected (RA = 0)
      n <- n[n != 0]
      #'  Calculate proportion of community each species represents
      N <- sum(n)
      p <- n/N
      #'  Calculate Shannon's diversity index (H)
      H[i] <- -sum(p * log(p))
    }
    Shannon <- cbind(Shannon, H) %>%
      mutate(H = round(H, 5))
    
    H_noLago <- c(NA)
    for(i in 1:nrow(Shannon_noLago)) {
      #'  Relative abundance of each species
      n <- c(Shannon_noLago[i,2], Shannon_noLago[i,3], Shannon_noLago[i,4], Shannon_noLago[i,5], Shannon_noLago[i,6]) 
      #'  Remove species that were not detected (RA = 0)
      n <- n[n != 0]
      #'  Calculate proportion of community each species represents
      N <- sum(n)
      p <- n/N
      #'  Calculate Shannon's diversity index (H)
      H_noLago[i] <- -sum(p * log(p))
    }
    Shannon_noLago <- cbind(Shannon_noLago, H_noLago) %>%
      mutate(H_noLago = round(H_noLago, 5)) 
    
    #'  List wild ungulate species detected most frequently at each camera
    dominantSpp <- Shannon %>% 
      # dplyr::select(-c(NewLocationID, lagomorphs_perday, livestock_perday, H)) %>%
      dplyr::select(-c(NewLocationID, lagomorphs, livestock, H)) %>%
      rowwise() %>%
      mutate(dominantprey = names(.)[which.max(c_across(everything()))],
             # dominantprey = ifelse(dominantprey == "moose_perday", "other", dominantprey),
             # dominantprey = ifelse(dominantprey == "muledeer_perday", "other", dominantprey)) %>%
             dominantprey = ifelse(dominantprey == "moose", "other", dominantprey),
             dominantprey = ifelse(dominantprey == "muledeer", "other", dominantprey)) %>%
      dplyr::select(dominantprey) %>%
      cbind(Shannon) %>%
      cbind(Shannon_noLago$H_noLago) %>%
      rename("H_noLago" = "Shannon_noLago$H_noLago")
    
    #'  Merge SR, H, and most frequently detected wild ungulate species with raw and 
    #'  weighted relative abundance indices
    Spp_diversity <- full_join(SR, dominantSpp, by = c("NewLocationID")) %>%
      relocate(SR, .after = NewLocationID) %>%
      relocate(H, .after = SR) %>%
      relocate(H_noLago, .after = H) %>%
      relocate(dominantprey, .after = H_noLago) %>%
      # mutate(dominantprey = gsub("_perday", "", dominantprey)) %>%
      dplyr::select(c(NewLocationID, SR, H, H_noLago, dominantprey)) %>%
      full_join(RA, by = "NewLocationID")
    
    return(Spp_diversity)
  }
  spp_diversity_Smr20 <- species_diversity(RA_Smr20_df)
  spp_diversity_Wtr20 <- species_diversity(RA_Wtr20_df)
  spp_diversity_Smr21 <- species_diversity(RA_Smr21_df)
  
  