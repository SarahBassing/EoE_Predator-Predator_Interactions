  #'  -------------------------------------
  #'  Result tables and figures
  #'  Northern Idaho Predator-Prey Project
  #'  Sarah B. Bassing
  #'  December 2024
  #'  -------------------------------------
  #'  Create result tables and figures based on the best supported model for each species and month
  #'  -------------------------------------

  #'  Load libraries
  # library(jagsUI)
  library(stringr)
  library(tidyverse)
    
  #'  Load top models
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/elk_july_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/elk_aug_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/wtd_july_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/wtd_aug_topmods.RData")
  
  #'  Add names to lists 
  elk_july_modnames <- c("RN_elk_july_null", "RN_elk_july_selected")
  # elk_july_modnames <- c("RN_elk_july_null", "RN_elk_july_selected", "RN_elk_july_selected.propSelected")
  names(elk_july_topmods) <- elk_july_modnames
  elk_aug_modnames <- c("RN_elk_aug_null", "RN_elk_aug_predicted.propSelected")
  # elk_aug_modnames <- c("RN_elk_aug_null", "RN_elk_aug_cvHQ", "RN_elk_aug_cvTbio", "RN_elk_aug_predicted.propSelected")
  names(elk_aug_topmods) <- elk_aug_modnames
  wtd_july_modnames <- c("RN_wtd_july_global1")
  names(wtd_july_topmods) <- wtd_july_modnames
  wtd_aug_modnames <- c("RN_wtd_aug_global1", "RN_wtd_aug_predicted.propSelected")
  # wtd_aug_modnames <- c("RN_wtd_aug_predicted.propSelected", "RN_wtd_aug_global1")
  names(wtd_aug_topmods) <- wtd_aug_modnames
  
  #' #'  Load saved detection histories (list order: July [[1]], Aug [[2]])
  #' load("./Data/Side_quests/Hilger/DH_elk_RNmod.RData")
  #' load("./Data/Side_quests/Hilger/DH_wtd_RNmod.RData")
  
  #'  Load saved rownames for annual detection histories
  load("./Data/Side_quests/Hilger/stacked_rownames_elk.RData")
  load("./Data/Side_quests/Hilger/stacked_rownames_wtd.RData")
  
  #'  Load saved covariate data (stacked to mirror detection histories)
  load("./Data/Side_quests/Hilger/stacked_stations_elk_july_aug.RData")
  load("./Data/Side_quests/Hilger/stacked_stations_wtd_july_aug.RData")
  
  #'  Grab estimated N per site from each top model
  estimated_N <- function(mod, locs, covs, spp) {
    #'  Grab estimated N and SD per site
    RN.n <- mod$mean$N
    RN.sd <- mod$sd$N
    #' #'  Grab camera location
    #' locs <- locs
    #'  Merge and format into single data frame with corresponding N & SD per site
    out <- cbind(locs, RN.n, RN.sd)
    RN_est <- as.data.frame(out) %>%
      mutate(NewLocationID = NewLocationID,
             Season = Season,
             Species = spp,
             RN.n = as.numeric(RN.n),
             RN.sd = as.numeric(RN.sd)) %>% 
      right_join(covs, by = c("NewLocationID", "Season")) %>%
      dplyr::select(c(NewLocationID, Species, Season, GMU, Setup, RN.n, RN.sd)) %>%
      mutate(GMU = ifelse(GMU == 1, "GMU10A", GMU),
             GMU = ifelse(GMU == 2, "GMU6", GMU),
             GMU = ifelse(GMU == 3, "GMU1", GMU),
             Setup = ifelse(Setup == 1, "U", "P"))
    return(RN_est)
  }
  rn_elk_july_out <- lapply(elk_july_topmods, estimated_N, locs = stacked_rownames_elk, covs = station_elk_list[[1]], spp = "elk")
  rn_elk_aug_out <- lapply(elk_aug_topmods, estimated_N, locs = stacked_rownames_elk, covs = station_elk_list[[2]], spp = "elk")
  rn_wtd_july_out <- lapply(wtd_july_topmods, estimated_N, locs = stacked_rownames_wtd, covs = station_wtd_list[[1]], spp = "whitetailed_deer")
  rn_wtd_aug_out <- lapply(wtd_aug_topmods, estimated_N, locs = stacked_rownames_wtd, covs = station_wtd_list[[2]], spp = "whitetailed_deer")
  
  #'  Site-specific N from top model
  rn_elk_july_out_top <- rn_elk_july_out[[1]]
  rn_elk_aug_out_top <- rn_elk_aug_out[[1]]
  rn_wtd_july_out_top <- rn_wtd_july_out[[1]]
  rn_wtd_aug_out_top <- rn_wtd_aug_out[[1]]
  
  #'  Save site-specific indices of relative abundances - use for mapping
  write_csv(rn_elk_july_out_top, "./Outputs/Hilger_RNmodel/Tables/RN_abundance_elk_july_top_mods.csv")
  write_csv(rn_elk_aug_out_top, "./Outputs/Hilger_RNmodel/Tables/RN_abundance_elk_aug_top_mods.csv")
  write_csv(rn_wtd_july_out_top, "./Outputs/Hilger_RNmodel/Tables/RN_abundance_wtd_july_top_mods.csv")
  write_csv(rn_wtd_aug_out_top, "./Outputs/Hilger_RNmodel/Tables/RN_abundance_wtd_aug_top_mods.csv")
  
  #'  Define number of significant digits to round to
  rounddig <- 3
  
  #'  ----------------------------------------------
  ####  Average parameter estimates in wide format  ####
  #'  ----------------------------------------------
  #'  Model outputs (mean lambda, p, and r) per species
  params <- function(mod, spp, modname) {
    #'  Grab parameter estimates and SD per site
    lambda <- round(mod$mean$mu.lambda, 2)               
    lambda.ll <- round(mod$q2.5$mu.lambda, rounddig)
    lambda.ul <- round(mod$q97.5$mu.lambda, rounddig)     ### Could grab lambdaYR estimates as well  ###
    r <- round(mod$mean$mu.r, rounddig)
    r.ll <- round(mod$q2.5$mu.r, rounddig)
    r.ul <- round(mod$q97.5$mu.r, rounddig)
    p <- round(mod$mean$mean.p, rounddig)
    p.ll <- round(mod$q2.5$mean.p, rounddig)
    p.ul <- round(mod$q97.5$mean.p, rounddig)
    #'  Merge and format into single data frame
    out <- cbind(lambda, lambda.ll, lambda.ul, r, r.ll, r.ul, p, p.ll, p.ul)
    param_est <- as.data.frame(out) %>%
      mutate(Species = spp,
             Model = modname,
             # lambda = as.numeric(lambda),
             lambda.cri = paste0(lambda, " (", lambda.ll, " - ", lambda.ul, ")"),
             # r = as.numeric(r),
             r.cri = paste0(r, " (", r.ll, " - ", r.ul, ")"),
             # p = as.numeric(p),
             p.cri = paste0(p, " (", p.ll, " - ", p.ul, ")")) %>%
      dplyr::select(c(Species, Model, lambda.cri, r.cri, p.cri))
    return(param_est)
  }
  elk_july_mean_params <- mapply(params, elk_july_topmods, spp = "Elk, July", modname = elk_july_modnames, SIMPLIFY = FALSE) %>% bind_rows(.)
  elk_aug_mean_params <- mapply(params, elk_aug_topmods, spp = "Elk, August", modname = elk_aug_modnames, SIMPLIFY = FALSE) %>% bind_rows(.)
  wtd_july_mean_params <- mapply(params, wtd_july_topmods, spp = "White-tailed deer, July", modname = wtd_july_modnames, SIMPLIFY = FALSE) %>% bind_rows(.)
  wtd_aug_mean_params <- mapply(params, wtd_aug_topmods, spp = "White-tailed deer, August", modname = wtd_aug_modnames, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  mean_params <- bind_rows(elk_july_mean_params, elk_aug_mean_params, wtd_july_mean_params, wtd_aug_mean_params)
  write_csv(mean_params, "./Outputs/Hilger_RNmodel/Tables/RN_mean_lambda_table.csv")
  
  #'  --------------------------------------------
  ####  All coefficient estimates in long format  ####
  #'  --------------------------------------------
  mod_out_summary_table <- function(mod, spp, mo) {
    #'  Retain & reformat model coefficients and mean derived parameters (r, lambda, psi)
    out <- as.data.frame(mod$summary) %>%
      rownames_to_column(., "Parameter") %>%
      transmute(
        Species = spp,
        Month = mo,
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
        Submodel = ifelse(Parameter == "alpha0" | Parameter == "a.setup[1]" | Parameter == "a.setup[2]", "Detection", "Abundance"),
        Submodel = ifelse(Parameter == "mu.lambda" | Parameter == "mu.r" | Parameter == "mean.p", NA, Submodel),
        Parameter = ifelse(Parameter == "beta0", "Intercept [Year: 2020]", Parameter),
        Parameter = ifelse(Parameter == "b.year[2]", "Year: 2021", Parameter),
        Parameter = ifelse(Parameter == "b.year[3]", "Year: 2022", Parameter),
        Parameter = ifelse(Parameter == "b.selected", "Total selected species", Parameter),
        Parameter = ifelse(Parameter == "b.predicted", "Total predicted species", Parameter),
        Parameter = ifelse(Parameter == "b.prop.selected", "Proportion selected species", Parameter),
        Parameter = ifelse(Parameter == "b.meanTbio", "Mean total biomass", Parameter),
        Parameter = ifelse(Parameter == "b.maxTbio", "Maximum total biomass", Parameter),
        Parameter = ifelse(Parameter == "b.cvTbio", "CV total biomass", Parameter),
        Parameter = ifelse(Parameter == "b.meanHQ", "Mean high quality biomass", Parameter),
        Parameter = ifelse(Parameter == "b.maxHQ", "Maximum high quality biomass", Parameter),
        Parameter = ifelse(Parameter == "b.cvHQ", "CV high quality biomass", Parameter),
        Parameter = ifelse(Parameter == "alpha0", "Intercept [Random camera]", Parameter),
        Parameter = ifelse(Parameter == "a.setup[2]", "Setup [Road camera]", Parameter),     #a.setup[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "mu.lambda", "Mean Lambda", Parameter),
        Parameter = ifelse(Parameter == "mu.r", "Mean Pr(Ind. detection)", Parameter),
        Parameter = ifelse(Parameter == "mean.p", "Mean Pr(Detection)", Parameter),
        overlap0 = ifelse(str_detect(Parameter, "Intercept "), NA, overlap0))  %>%
      relocate(Submodel, .before = Parameter) %>%
      #'  Save only coefficients and mean derived parameters
      filter(!grepl("N", Parameter) & !grepl("loglike.new", Parameter) & !grepl("lambda", Parameter) &
               !grepl("deviance", Parameter) & !grepl("rSetup", Parameter) & !grepl("b.year", Parameter) & 
               !grepl("a.setup", Parameter)) %>%
      filter(Parameter != "Mean Lambda" & Parameter != "Mean Pr(Ind. detection)" & Parameter != "Mean Pr(Detection)") %>%
      rename("Coefficient" = "Parameter") %>%
      arrange(Submodel)
    
    print(out)
    return(out)
  }
  #'  Extract coefficient estimates from top model (or 2nd best if NULL was top)
  elk_july_coef_tbl <- mod_out_summary_table(elk_july_topmods[[2]], spp = "Elk", mo = "July") 
  elk_aug_coef_tbl <- mod_out_summary_table(elk_aug_topmods[[2]], spp = "Elk", mo = "Aug") 
  wtd_july_coef_tbl <- mod_out_summary_table(wtd_july_topmods[[1]], spp = "White-tailed deer", mo = "July") 
  wtd_aug_coef_tbl <- mod_out_summary_table(wtd_aug_topmods[[1]], spp = "White-tailed deer", mo = "Aug") 
  
  RN_result_tbl <- rbind(elk_july_coef_tbl, elk_aug_coef_tbl, wtd_july_coef_tbl, wtd_aug_coef_tbl) 
  
  #'  Save
  write_csv(RN_result_tbl, file = "./Outputs/Hilger_RNmodel/Tables/RN_result_tbl.csv")
  
  #'  ------------------------
  ####  Data summary tables  ####
  #'  ------------------------
  #'  Unique detection events
  load("./Data/Side_quests/Hilger/DetectionEvents_eoe20s.RData")
  load("./Data/Side_quests/Hilger/DetectionEvents_eoe21s.RData")
  load("./Data/Side_quests/Hilger/DetectionEvents_eoe22s.RData")
  #'  Seasonal detection histories
  load("./Data/Side_quests/Hilger/DH_eoe20s_RNmod.RData") 
  load("./Data/Side_quests/Hilger/DH_eoe21s_RNmod.RData") 
  load("./Data/Side_quests/Hilger/DH_eoe22s_RNmod.RData") 
  #'  Split by month and name lists by species and month
  split_DH <- function(dh) {
    july <- dh[,1:31]
    august <- dh[,32:62]
    month_list <- list(july, august)
    names(month_list) <- c("July", "August")
    return(month_list)
  }
  #'  Monthly list order: Elk[[1]], White-tailed Deer[[2]]; July[[1]], August[[2]]
  spp_order <- c("Elk", "WTD")
  DH_eoe20s_monthly <- lapply(DH_eoe20s_RNmod, split_DH); names(DH_eoe20s_monthly) <- spp_order
  DH_eoe21s_monthly <- lapply(DH_eoe21s_RNmod, split_DH); names(DH_eoe21s_monthly) <- spp_order
  DH_eoe22s_monthly <- lapply(DH_eoe22s_RNmod, split_DH); names(DH_eoe22s_monthly) <- spp_order
  #' Sampling effort
  load("./Data/Side_quests/Hilger/Effort_eoe20s_RNmod.RData")
  load("./Data/Side_quests/Hilger/Effort_eoe21s_RNmod.RData")
  load("./Data/Side_quests/Hilger/Effort_eoe22s_RNmod.RData")
  
  #'  ----------------------------------
  #####  Summary stats for publication  #####
  #'  ----------------------------------
  #'  Number of cameras used in analyses and length of primary period
  dim(DH_eoe20s_RNmod[[1]])
  dim(DH_eoe21s_RNmod[[1]])
  dim(DH_eoe22s_RNmod[[1]])
  
  #'  Average number of operable camera days
  (July_ndays20 <- mean(Effort_eoe20s_RNmod$July_nDays)); sd(Effort_eoe20s_RNmod$July_nDays)/sqrt(length(Effort_eoe20s_RNmod$July_nDays))
  (Aug_ndays20 <- mean(Effort_eoe20s_RNmod$August_nDays)); sd(Effort_eoe20s_RNmod$August_nDays)/sqrt(length(Effort_eoe20s_RNmod$August_nDays))
  (July_ndays21 <- mean(Effort_eoe21s_RNmod$July_nDays)); sd(Effort_eoe21s_RNmod$July_nDays)/sqrt(length(Effort_eoe21s_RNmod$July_nDays))
  (Aug_ndays21 <- mean(Effort_eoe21s_RNmod$August_nDays)); sd(Effort_eoe21s_RNmod$August_nDays)/sqrt(length(Effort_eoe21s_RNmod$August_nDays))
  (July_ndays22 <- mean(Effort_eoe22s_RNmod$July_nDays)); sd(Effort_eoe22s_RNmod$July_nDays)/sqrt(length(Effort_eoe22s_RNmod$July_nDays))
  (Aug_ndays22 <- mean(Effort_eoe22s_RNmod$August_nDays)); sd(Effort_eoe22s_RNmod$August_nDays)/sqrt(length(Effort_eoe22s_RNmod$August_nDays))
  July_nDays <- c(July_ndays20, July_ndays21, July_ndays22)
  Aug_nDays <- c(Aug_ndays20, Aug_ndays21, Aug_ndays22)
  mean(July_nDays); sd(July_nDays)/sqrt(length(July_nDays))
  mean(Aug_nDays); sd(Aug_nDays)/sqrt(length(Aug_nDays))
  
  #'  Total number of unique detection events per species of interest
  dets_20s <- mutate(eoe20s_det_events, Year = "2020")
  dets_21s <- mutate(eoe21s_det_events, Year = "2021")
  dets_22s <- mutate(eoe22s_det_events, Year = "2022")
  dets <- bind_rows(dets_20s, dets_21s, dets_22s) %>%
    group_by(Species, Year) %>%
    summarise(ndets = n()) %>%
    ungroup() %>%
    mutate(Species = ifelse(Species == "elk", "Elk", Species), 
           Species = ifelse(Species == "whitetaileddeer", "White-tailed deer", Species))
  
  #'  -----------------------------
  #####  Summarize detection data  ####
  #'  ------------------------------
  #'  Function to summarize species-specific detection data (from detection histories)
  all_detections <- function(dh1, dh2, dh3, locs, spp) {    ###  COULD SPLIT THIS OUT BY MONTH AS WELL  ###
    #'  Bind data and summarize detection events
    dh <- bind_rows(dh1, dh2, dh3) %>%
      #'  Sum detection events across sampling occasions
      mutate(total_dets = rowSums(., na.rm = TRUE),
             #'  Binary whether a detection ever occurred
             dets_YN = ifelse(total_dets > 1, 1, total_dets)) %>%
      bind_cols(locs)
    #'  Summarize detection data per year
    det_summary <- dh %>%
      group_by(Season) %>%
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
      mutate(Species = spp,
             Season = ifelse(Season == "Smr20", "2020", Season),
             Season = ifelse(Season == "Smr21", "2021", Season),
             Season = ifelse(Season == "Smr22", "2022", Season)) %>%
      dplyr::select(Species, Season, n_dets, n_cam_dets, n_cams, prop_cam_dets)
    colnames(det_summary) <- c("Species", "Year", "Total detection events", "Total cameras with detections", 
                               "Number of operating cameras", "Proportion cameras with detections")
    return(det_summary)
  }
  DH_elk <- all_detections(DH_eoe20s_RNmod[[1]], DH_eoe21s_RNmod[[1]], DH_eoe22s_RNmod[[1]], locs = stacked_rownames_elk, spp = "Elk")
  DH_wtd <- all_detections(DH_eoe20s_RNmod[[2]], DH_eoe21s_RNmod[[2]], DH_eoe22s_RNmod[[2]], locs = stacked_rownames_wtd, spp = "White-tailed deer")
  
  #'  Final detection history summary table
  DH_summary <- rbind(DH_elk, DH_wtd) %>%
    #'  Swap total detection events (based on number of days with detections) for total raw detection events
    full_join(dets, by = c("Species", "Year")) %>%
    dplyr::select(-c(`Total detection events`, `Number of operating cameras`)) %>%
    relocate(ndets, .after = Year) %>%
    rename("Total detection events" = "ndets")
  
  #'  Calculate mean number of detection events per species per year
  DH_Elk <- DH_summary[DH_summary$Species == "Elk",]
  (mean_dets_elk <- mean(DH_Elk$`Total detection events`))
  (se_dets_elk <- sd(DH_Elk$`Total detection events`)/sqrt(nrow(DH_Elk)))
  DH_WTD <- DH_summary[DH_summary$Species == "White-tailed deer",]
  (mean_dets_wtd <- mean(DH_WTD$`Total detection events`))
  (se_dets_wtd <- sd(DH_WTD$`Total detection events`)/sqrt(nrow(DH_WTD)))
  
  
  #'  Save!
  write.csv(DH_summary, "./Outputs/Hilger_RNmodel/Tables/Summary_table_detections.csv")
  
  
  #'  -------------------------
  ####  Partial effects plots  ####
  #'  -------------------------
  #'  Grab covariates in top models (List order: July [[1]], August [[2]])
  elk_july_selected <- station_elk_list[[1]]$total_selected_species
  elk_aug_predicted <- station_elk_list[[2]]$total_predicted_species
  elk_aug_propSelected <- station_elk_list[[2]]$prop_selected_species_weighted
  wtd_july_cvHQ <- station_wtd_list[[1]]$cv_HQ
  wtd_july_selected <- station_wtd_list[[1]]$total_selected_species
  wtd_july_propSelected <- station_wtd_list[[1]]$prop_selected_species_weighted
  wtd_aug_meanTbio <- station_wtd_list[[2]]$mean_Tbio_kg_ha
  wtd_aug_cvHQ <- station_wtd_list[[2]]$cv_HQ
  wtd_aug_selected <- station_wtd_list[[2]]$total_selected_species
  wtd_aug_propSelected <- station_wtd_list[[2]]$prop_selected_species_weighted
  
  #'  Define number of new observations to predict across
  npoints <- 100
  
  #'  Create and z-transform new data based on covariates in top models
  newcovs <- function(focal_cov, npoints) {
    #'  Create range of covariate values to predict across
    print(r <- range(focal_cov))
    print(mean.cov <- mean(focal_cov))
    print(sd.cov <- sd(focal_cov))
    cov.pred.orig <- seq(r[1], r[2], length.out = npoints)
    #'  Scale new data to be consistent with data used in model
    cov.new <- (cov.pred.orig - mean.cov) / sd.cov
    return(cov.new)
  }
  elk_july_z.covs <- newcovs(elk_july_selected, npoints = npoints) 
  elk_aug_covs <- list(elk_aug_predicted, elk_aug_propSelected)
  elk_aug_z.covs <- lapply(elk_aug_covs, newcovs, npoints = npoints)
  wtd_july_covs <- list(wtd_july_cvHQ, wtd_july_selected, wtd_july_propSelected)
  wtd_july_z.covs <- lapply(wtd_july_covs, newcovs, npoints = npoints)
  wtd_aug_covs <- list(wtd_aug_meanTbio, wtd_aug_cvHQ, wtd_aug_selected, wtd_aug_propSelected)
  wtd_aug_z.covs <- lapply(wtd_aug_covs, newcovs, npoints = npoints)
  
  #'  Create a matrix for each *focal covariate* (newdata, z-transformed) in the order they appear in 
  #'  the top species & season-specific regression where Intercept = 1, Year = 3 (representing 2022), 
  #'  and all other covariates are held at their mean = 0
  #'  Top Elk July model: intercept, year, *total selected* (3 columns)
  elk_july_z.selected_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), elk_july_z.covs), nrow = npoints, ncol = 3)
  #'  Top Elk Aug model: intercept, year, *total predicted*, *proportion selected* (4 columns)
  elk_aug_z.predicted_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), elk_aug_z.covs[[1]], rep(0,npoints)), nrow = npoints, ncol = 4)
  elk_aug_z.propSelected_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), rep(0,npoints), elk_aug_z.covs[[2]]), nrow = npoints, ncol = 4)
  #'  Top WTD July model: intercept, year, mean Tbio, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 columns)
  wtd_july_z.cvHQ_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), rep(0,npoints), rep(0,npoints), wtd_july_z.covs[[1]], rep(0,npoints), rep(0,npoints)), nrow = npoints, ncol = 7)
  wtd_july_z.selected_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), wtd_july_z.covs[[2]], rep(0,npoints)), nrow = npoints, ncol = 7)
  wtd_july_z.propSelected_matrix <- matrix(data = c(rep(1,npoints), rep(3,npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), wtd_july_z.covs[[3]]), nrow = npoints, ncol = 7)
  #'  Top WTD Aug model: intercept, year, *mean Tbio*, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 columns)
  wtd_aug_z.meanTbio_matrix <- matrix(data = c(rep(1,npoints), rep(3, npoints), wtd_aug_z.covs[[1]], rep(0,npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints)), nrow = npoints, ncol = 7)
  wtd_aug_z.cvHQ_matrix <- matrix(data = c(rep(1,npoints), rep(3, npoints), rep(0,npoints), rep(0,npoints), wtd_aug_z.covs[[2]], rep(0,npoints), rep(0,npoints)), nrow = npoints, ncol = 7)
  wtd_aug_z.selected_matrix <- matrix(data = c(rep(1,npoints), rep(3, npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), wtd_aug_z.covs[[3]], rep(0,npoints)), nrow = npoints, ncol = 7)
  wtd_aug_z.propSelected_matrix <- matrix(data = c(rep(1,npoints), rep(3, npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), rep(0,npoints), wtd_aug_z.covs[[4]]), nrow = npoints, ncol = 7)
  
  
  #'  Predict across each iteration of coefficients from the top model
  
  
  
  
  
  
  
  
  
  
  