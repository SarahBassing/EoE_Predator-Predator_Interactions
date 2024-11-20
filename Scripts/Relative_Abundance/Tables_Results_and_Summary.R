  #'  ----------------------------------
  #'  Royle-Nichols model result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2023
  #'  ----------------------------------
  #'  Script to summarize results from species and year specific Royle-Nichols
  #'  abundance models.
  #'  ----------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(stringr)
  library(tidyverse)
    
  #'  Load model outputs (make sure only one run per model is in this folder)
  filenames <- list.files("./Outputs/Relative_Abundance/RN_model/JAGS_out", pattern="*.RData", full.names=TRUE)
  lapply(filenames, load, environment())
  
  #'  List models based on species and year
  yr_list <- list("2020", "2021", "2022")
  bear_list <- list(RN_bear_20s, RN_bear_21s, RN_bear_22s)
  bob_list <- list(RN_bob_20s, RN_bob_21s, RN_bob_22s)
  coy_list <- list(RN_coy_20s, RN_coy_21s, RN_coy_22s)
  lion_list <- list(RN_lion_20s, RN_lion_21s, RN_lion_22s)
  wolf_list <- list(RN_wolf_20s, RN_wolf_21s, RN_wolf_22s)
  elk_list <- list(RN_elk_20s, RN_elk_21s, RN_elk_22s)
  moose_list <- list(RN_moose_20s, RN_moose_21s, RN_moose_22s)
  wtd_list <- list(RN_wtd_20s, RN_wtd_21s, RN_wtd_22s)
  lago_list <- list(RN_lago_20s, RN_lago_21s, RN_lago_22s)
  
  #'  Define number of significant digits to round to
  rounddig <- 3
  
  #'  ----------------------------------------------
  ####  Average parameter estimates in wide format  ####
  #'  ----------------------------------------------
  #'  Model outputs (mean lambda, p, and r) per species
  params <- function(mod, spp, yr) {
    #'  Grab parameter estimates and SD per site
    lambda <- round(mod$mean$mu.lambda, 2)
    lambda.ll <- round(mod$q2.5$mu.lambda, rounddig)
    lambda.ul <- round(mod$q97.5$mu.lambda, rounddig)
    r <- round(mod$mean$mu.r, rounddig)
    r.ll <- round(mod$q2.5$mu.r, rounddig)
    r.ul <- round(mod$q97.5$mu.r, rounddig)
    p <- round(mod$mean$mean.p, rounddig)
    p.ll <- round(mod$q2.5$mean.p, rounddig)
    p.ul <- round(mod$q97.5$mean.p, rounddig)
    psi <- round(mod$mean$mean.psi, rounddig)
    psi.ll <- round(mod$q2.5$mean.psi, rounddig)
    psi.ul <- round(mod$q97.5$mean.psi, rounddig)
    #'  Merge and format into single data frame
    out <- cbind(lambda, lambda.ll, lambda.ul, r, r.ll, r.ul, p, p.ll, p.ul, psi, psi.ll, psi.ul)
    param_est <- as.data.frame(out) %>%
      mutate(Species = spp,
             Year = yr,
             # lambda = as.numeric(lambda),
             lambda.cri = paste0(lambda, " (", lambda.ll, " - ", lambda.ul, ")"),
             # r = as.numeric(r),
             r.cri = paste0(r, " (", r.ll, " - ", r.ul, ")"),
             # p = as.numeric(p),
             p.cri = paste0(p, " (", p.ll, " - ", p.ul, ")"),
             # psi = as.numeric(psi),
             psi.cri = paste0(psi, " (", psi.ll, " - ", psi.ul, ")")) %>%
      dplyr::select(c(Species, Year, lambda.cri, r.cri, p.cri, psi.cri))
    return(param_est)
  }
  study_yr <- list("2020", "2021", "2022")
  bear_mean_params <- mapply(params, bear_list, spp = "Black bear", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  bob_mean_params <- mapply(params, bob_list, spp = "Bobcat", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  coy_mean_params <- mapply(params, coy_list, spp = "Coyote", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  lion_mean_params <- mapply(params, lion_list, spp = "Mountain lion", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  wolf_mean_params <- mapply(params, wolf_list, spp = "Wolf", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  elk_mean_params <- mapply(params, elk_list, spp = "Elk", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  lago_mean_params <- mapply(params, lago_list, spp = "Lagomorph", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  moose_mean_params <- mapply(params, moose_list, spp = "Moose", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  wtd_mean_params <- mapply(params, wtd_list, spp = "White-tailed deer", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  mean_params <- bind_rows(bear_mean_params, bob_mean_params, coy_mean_params, lion_mean_params, wolf_mean_params, elk_mean_params, lago_mean_params, moose_mean_params, wtd_mean_params)
  write_csv(mean_params, "./Outputs/Relative_Abundance/RN_model/Tables/RN_mean_lambda_table.csv")
  
  #'  --------------------------------------------
  ####  All coefficient estiamtes in long format  ####
  #'  --------------------------------------------
  mod_out_summary_table <- function(mod, spp, yr) {
    #'  Retain & reformat model coefficients and mean derived parameters (r, lambda, psi)
    out <- as.data.frame(mod$summary) %>%
      rownames_to_column(., "Parameter") %>%
      transmute(
        Species = spp,
        Year = yr,
        Parameter = Parameter,
        Mean = round(mean, rounddig),
        Lower_CRI = round(`2.5%`, rounddig),
        Upper_CRI = round(`97.5%`, rounddig),
        Mean = as.numeric(Mean),
        Lower_CRI = as.numeric(Lower_CRI),
        Upper_CRI = as.numeric(Upper_CRI),
        `95% CRI` = paste0(" ", Lower_CRI, " - ", Upper_CRI),
        overlap0 = ifelse(overlap0 == 0, FALSE, TRUE)) %>%
      dplyr::select(-c(Lower_CRI, Upper_CRI)) %>%
      #'  Give parameters more meaningful names
      mutate(
        Parameter = ifelse(Parameter == "beta0", "Intercept [GMU10A]", Parameter),
        # Parameter = ifelse(Parameter == "beta1", "Percent forest", Parameter),
        # Parameter = ifelse(Parameter == "beta2", "Elevation", Parameter),
        # Parameter = ifelse(Parameter == "beta3", "Elevation^2", Parameter),
        Parameter = ifelse(Parameter == "beta4[2]", "GMU6", Parameter),              #beta4[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "beta4[3]", "GMU1", Parameter),
        Parameter = ifelse(Parameter == "alpha0", "Intercept [Random]", Parameter),
        Parameter = ifelse(Parameter == "alpha2[2]", "Setup [Road]", Parameter),     #alpha2[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "mu.lambda", "mean lambda", Parameter),
        Parameter = ifelse(Parameter == "mu.r", "mean Pr(ind. detection)", Parameter),
        Parameter = ifelse(Parameter == "mean.p", "mean Pr(detection)", Parameter),
        Parameter = ifelse(Parameter == "mean.psi", "mean Pr(occupancy)", Parameter),
        overlap0 = ifelse(str_detect(Parameter, "Intercept "), NA, overlap0))  %>%
      #'  Save only coefficients and mean derived parameters
      filter(Parameter == "Intercept [GMU10A]" | Parameter == "GMU6" | Parameter == "GMU1" | 
               # Parameter == "Percent forest" | Parameter == "Elevation" | Parameter == "Elevation^2" | 
               Parameter == "Intercept [Random]" | Parameter == "Setup [Road]" | 
               Parameter == "mean lambda" | Parameter == "mean Pr(ind. detection)" | 
               Parameter == "mean Pr(detection)" | Parameter == "mean Pr(occupancy)")
    
    print(out)
    return(out)
  }
  bear_tbl <- mapply(mod_out_summary_table, mod = bear_list, spp = "Black bear", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  bob_tbl <- mapply(mod_out_summary_table, mod = bob_list, spp = "Bobcat", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  coy_tbl <- mapply(mod_out_summary_table, mod = coy_list, spp = "Coyote", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  lion_tbl <- mapply(mod_out_summary_table, mod = lion_list, spp = "Mountain lion", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  wolf_tbl <- mapply(mod_out_summary_table, mod = wolf_list, spp = "Wolf", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  elk_tbl <- mapply(mod_out_summary_table, mod = elk_list, spp = "Elk", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  moose_tbl <- mapply(mod_out_summary_table, mod = moose_list, spp = "Moose", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  wtd_tbl <- mapply(mod_out_summary_table, mod = wtd_list, spp = "White-tailed deer", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  lago_tbl <- mapply(mod_out_summary_table, mod = lago_list, spp = "Lagomorphs", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  
  RN_result_tbl <- rbind(bear_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl) #bob_tbl, lago_tbl
  RN_result_tbl_list <- list(bear_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl) #bob_tbl, lago_tbl
  
  #'  Save
  write_csv(RN_result_tbl, file = "./Outputs/Relative_Abundance/RN_model/Tables/RN_result_tbl.csv")
  save(RN_result_tbl_list, file = "./Outputs/Relative_Abundance/RN_model/Tables/RN_result_tbl_list.R")
  
  #'  ------------------------
  ####  Data summary tables  ####
  #'  ------------------------
  #'  Unique detections
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_20s.RData") 
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_21s.RData") 
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_22s.RData") 
  #'  Seasonal detection histories
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  #'  Name lists so life is easier
  spp_order <- c("bear_black", "bobcat", "coyote", "mountain_lion", "wolf", "elk", "moose", "muledeer", "whitetaileddeer", "rabbit_hare")
  names(DH_npp20s_RNmod) <- spp_order
  names(DH_npp21s_RNmod) <- spp_order
  names(DH_npp22s_RNmod) <- spp_order
  #' Sampling effort
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp20s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp21s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp22s_RNmod.RData")
  
  #'  ----------------------------------
  #####  Summary stats for publication  #####
  #'  ----------------------------------
  #'  Number of cameras used in analyses and length of primary period
  dim(DH_npp20s_RNmod[[1]])
  dim(DH_npp21s_RNmod[[1]])
  dim(DH_npp22s_RNmod[[1]])
  
  #'  Average number of operable camera days
  (ndays20 <- mean(effort_20s_RNmod$ndays)); sd(effort_20s_RNmod$ndays)/sqrt(length(effort_20s_RNmod$ndays))
  (ndays21 <- mean(effort_21s_RNmod$ndays)); sd(effort_21s_RNmod$ndays)/sqrt(length(effort_21s_RNmod$ndays))
  (ndays22 <- mean(effort_22s_RNmod$ndays)); sd(effort_22s_RNmod$ndays)/sqrt(length(effort_22s_RNmod$ndays))
  ndays <- c(ndays20, ndays21, ndays22)
  mean(ndays)
  sd(ndays)/sqrt(length(ndays))
  
  #'  Average number of detection events per species of interest
  dets_20s <- mutate(npp20s_det_events, Year = "2020")
  dets_21s <- mutate(npp21s_det_events, Year = "2021")
  dets_22s <- mutate(npp22s_det_events, Year = "2022")
  dets <- bind_rows(dets_20s, dets_21s, dets_22s) %>%
    filter(Species == "bear_black" | Species == "coyote" | Species == "mountain_lion" | Species == "wolf" | 
             Species == "elk" | Species == "moose" | Species == "whitetaileddeer") %>%
    group_by(Species, Year) %>%
    summarise(ndets = n()) %>%
    ungroup() %>%
    mutate(Species = ifelse(Species == "bear_black", "Black bear", Species),
           Species = ifelse(Species == "coyote", "Coyote", Species), 
           Species = ifelse(Species == "mountain_lion", "Mountain lion", Species), 
           Species = ifelse(Species == "wolf", "Wolf", Species),
           Species = ifelse(Species == "elk", "Elk", Species), 
           Species = ifelse(Species == "moose", "Moose", Species), 
           Species = ifelse(Species == "whitetaileddeer", "White-tailed deer", Species))
    
  #'  -----------------------------
  #####  Summarize detection data  ####
  #'  ------------------------------
  #'  Function to summarize species-specific detection data
  all_detections <- function(dh1, dh2, dh3, spp) {
    #'  Add year to each data set
    dh1$Year <- "2020"
    dh2$Year <- "2021"
    dh3$Year <- "2022"
    #'  Bind data and summarize detections
    dh <- bind_rows(dh1, dh2, dh3) %>%
      rownames_to_column(., "NewLocationID") %>%
      relocate(Year, .before = o1) %>%
      #'  Sum detection events across sampling occasions
      mutate(total_dets = dplyr::select(., o1:o107) %>% rowSums(na.rm = TRUE),
             #'  Binary whether a detection ever occurred
             dets_YN = ifelse(total_dets > 1, 1, total_dets)) 
    #'  Summarize detection data per year
    det_summary <- dh %>%
      group_by(Year) %>%
      #'  Number of detection events
      summarize(n_dets = sum(total_dets),
                #'  Number of cameras with 1+ detection
                n_cam_dets = sum(dets_YN),
                #'  Number of unique camera locations
                n_cams = n(),
                #'  Proportion of cameras that had 1+ detection
                prop_cam_dets = round(n_cam_dets/n_cams, 2)) %>%
      ungroup() %>%
      #'  Add species name to table
      mutate(Species = spp) %>%
      relocate(Species, .before = Year)
    colnames(det_summary) <- c("Species", "Year", "Total detection events", "Total cameras with detections", 
                               "Number of operating cameras", "Proportion cameras with detections")
    return(det_summary)
  }
  DH_bear <- all_detections(DH_npp20s_RNmod$bear_black, DH_npp21s_RNmod$bear_black, DH_npp22s_RNmod$bear_black, spp = "Black bear")
  DH_bob <- all_detections(DH_npp20s_RNmod$bobcat, DH_npp21s_RNmod$bobcat, DH_npp22s_RNmod$bobcat, spp = "Bobcat")
  DH_coy <- all_detections(DH_npp20s_RNmod$coyote, DH_npp21s_RNmod$coyote, DH_npp22s_RNmod$coyote, spp = "Coyote")
  DH_lion <- all_detections(DH_npp20s_RNmod$mountain_lion, DH_npp21s_RNmod$mountain_lion, DH_npp22s_RNmod$mountain_lion, spp = "Mountain lion")
  DH_wolf <- all_detections(DH_npp20s_RNmod$wolf, DH_npp21s_RNmod$wolf, DH_npp22s_RNmod$wolf, spp = "Wolf")
  DH_elk <- all_detections(DH_npp20s_RNmod$elk, DH_npp21s_RNmod$elk, DH_npp22s_RNmod$elk, spp = "Elk")
  DH_moose <- all_detections(DH_npp20s_RNmod$moose, DH_npp21s_RNmod$moose, DH_npp22s_RNmod$moose, spp = "Moose")
  DH_wtd <- all_detections(DH_npp20s_RNmod$whitetaileddeer, DH_npp21s_RNmod$whitetaileddeer, DH_npp22s_RNmod$whitetaileddeer, spp = "White-tailed deer")
  
  #'  Final detection history summary table
  DH_summary <- rbind(DH_bear, DH_coy, DH_lion, DH_wolf, DH_elk, DH_moose, DH_wtd) %>% #DH_bob
    #'  Swap total detection events (based on number of days with detections) for total raw detection events
    full_join(dets, by = c("Species", "Year")) %>%
    dplyr::select(-c(`Total detection events`, `Number of operating cameras`)) %>%
    relocate(ndets, .after = Year) %>%
    rename("Total detection events" = "ndets")
  
  #'  Calculate mean number of detection events per species per year
  (mean_dets <- mean(DH_summary$`Total detection events`))
  (se_dets <- sd(DH_summary$`Total detection events`)/sqrt(nrow(DH_summary)))
  
  DH_pred <- DH_summary[DH_summary$Species == "Black bear" | DH_summary$Species == "Coyote" | DH_summary$Species == "Mountain lion" | DH_summary$Species == "Wolf",]
  (mean_dets_pred <- mean(DH_pred$`Total detection events`))
  (se_dets_pred <- sd(DH_pred$`Total detection events`)/sqrt(nrow(DH_pred)))
  DH_prey <- DH_summary[DH_summary$Species == "Elk" | DH_summary$Species == "Moose" | DH_summary$Species == "White-tailed Deer",]
  (mean_dets_prey <- mean(DH_prey$`Total detection events`))
  (se_dets_prey <- sd(DH_prey$`Total detection events`)/sqrt(nrow(DH_prey)))
  
  
  #'  Save!
  write.csv(DH_summary, "./Outputs/Relative_Abundance/RN_model/Tables/Summary_table_detections.csv")
  
  
  
  
  