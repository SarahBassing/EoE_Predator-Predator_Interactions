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
  
  #'  ------------------------------------------------
  ####  Load and format various datasets and outputs  ####
  #'  ------------------------------------------------
  
  #'  Load top models
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/elk_july_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/elk_aug_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/wtd_july_topmods.RData")
  load("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Top_models/wtd_aug_topmods.RData")
  
  #'  Add names to lists 
  elk_july_modnames <- c("RN_elk_july_null", "RN_elk_july_selected")
  names(elk_july_topmods) <- elk_july_modnames
  elk_aug_modnames <- c("RN_elk_aug_null", "RN_elk_aug_predicted.propSelected")
  names(elk_aug_topmods) <- elk_aug_modnames
  wtd_july_modnames <- c("RN_wtd_july_selected.propSelected", "RN_wtd_july_predicted.propSelected", "RN_wtd_july_global1", "RN_wtd_july_global2")
  names(wtd_july_topmods) <- wtd_july_modnames
  wtd_aug_modnames <- c("RN_wtd_aug_predicted.propSelected", "RN_wtd_aug_selected.propSelected", "RN_wtd_aug_predicted", "RN_wtd_aug_maxTbio", "RN_wtd_aug_global1", "RN_wtd_aug_selected", "RN_wtd_aug_global2", "RN_wtd_aug_Tbio_max.cv")
  names(wtd_aug_topmods) <- wtd_aug_modnames
  
  #'  Load saved rownames for annual detection histories
  load("./Data/Side_quests/Hilger/stacked_rownames_elk.RData")
  load("./Data/Side_quests/Hilger/stacked_rownames_wtd.RData")
  
  #'  Load saved covariate data (stacked to mirror detection histories)
  load("./Data/Side_quests/Hilger/stacked_stations_elk_july_aug.RData")
  load("./Data/Side_quests/Hilger/stacked_stations_wtd_july_aug.RData")
  
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
  
  
  #'  ---------------------------------
  ####  Retrieve estimated N per site  ####
  #'  ---------------------------------
  #'  Function to grab estimated N per site from each top model
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
  
  
  #'  ----------------------------------------------
  ####  Average parameter estimates in wide format  ####
  #'  ----------------------------------------------
  #'  Define number of significant digits to round to
  rounddig <- 3
  
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
  
  #'  ---------------------------------
  ####  Summary stats for publication  ####
  #'  ---------------------------------
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
  write_csv(DH_summary, "./Outputs/Hilger_RNmodel/Tables/Summary_table_detections.csv")
  
  
  #'  -----------------------------------------
  ####  Predict partial effects of covariates  ####
  #'  -----------------------------------------
  #'  Define number of new observations to predict across
  nobs <- 500
  
  #'  Create design matrix with new z-transform covariate data based on top model
  design_matrix <- function(focal_cov, nobs, lambda_cov, focal_cov_index) {
    #'  Create range of covariate values to predict across
    print(r <- range(focal_cov))
    print(mean.cov <- mean(focal_cov))
    print(sd.cov <- sd(focal_cov))
    cov.pred.orig <- seq(r[1], r[2], length.out = nobs)
    print(head(cov.pred.orig))
    #'  Scale new data to be consistent with data used in model
    cov.pred <- (cov.pred.orig - mean.cov) / sd.cov
    
    #'  Create design matrices for each focal covariate on lambda
    #'  Fill design matrix with categorical or mean value of each covariate
    lambda_cov <- lambda_cov
    lambda_covs <- matrix(lambda_cov, nrow = nobs, ncol = length(lambda_cov), byrow = TRUE)
    #'  Replace focal covariate column with new scaled data
    #'  New data must be in correct column based on order of covariates in original model
    lambda_covs[,focal_cov_index] <- cov.pred
    print(head(lambda_covs))
    
    #'  List new design matrix and the un-transformed new covariate value
    newcov_list <- list(lambda_covs, cov.pred.orig)
    
    return(newcov_list)
  }
  #'  Grab covariates in top models and create new z-transformed observations for 
  #'  each focal covariate (Station list order: July [[1]], August [[2]])
  #'  Intercept = 1; Year = 1 b/c plan to predict lambda for 2022 so need to turn on year effect
  #'  Top Elk July (non-null) model: intercept, year, *total selected* (3 columns)
  elk_july_z.selected_matrix <- design_matrix(station_elk_list[[1]]$selected, nobs = nobs, 
                                        lambda_cov = c(1, 1, 0), focal_cov_index = 3)
  #'  Top Elk Aug (non-null) model: intercept, year, *total predicted*, *proportion selected* (4 columns)
  elk_aug_z.predicted_matrix <- design_matrix(station_elk_list[[2]]$total, nobs = nobs, 
                                        lambda_cov = c(1, 1, 0, 0), focal_cov_index = 3)
  elk_aug_z.propSelected_matrix <- design_matrix(station_elk_list[[2]]$prop_selected, nobs = nobs, 
                                           lambda_cov = c(1, 1, 0, 0), focal_cov_index = 4)
  #'  Top WTD July model: intercept, year, *total selected*, *proportion selected* (4 columns) #mean Tbio, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 columns)
  # wtd_july_z.cvHQ_matrix <- design_matrix(station_wtd_list[[1]]$cv_HQ, nobs = nobs, 
  #                                   lambda_cov = c(1, 1, 0, 0, 0, 0, 0), focal_cov_index = 5)
  wtd_july_z.selected_matrix <- design_matrix(station_wtd_list[[1]]$selected, nobs = nobs, 
                                        lambda_cov = c(1, 1, 0, 0), focal_cov_index = 3)
  wtd_july_z.propSelected_matrix <- design_matrix(station_wtd_list[[1]]$prop_selected, nobs = nobs, 
                                            lambda_cov = c(1, 1, 0, 0), focal_cov_index = 4)
  #'  Top WTD Aug model: intercept, year, *total predicted*, *proportion selected* (4 columns) #*mean Tbio*, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 columns)
  # wtd_aug_z.meanTbio_matrix <- design_matrix(station_wtd_list[[2]]$mean_Tbio_kg_ha, nobs = nobs, 
  #                                      lambda_cov = c(1, 1, 0, 0, 0, 0, 0), focal_cov_index = 3)
  # wtd_aug_z.cvHQ_matrix <- design_matrix(station_wtd_list[[2]]$cv_HQ, nobs = nobs, 
  #                                  lambda_cov = c(1, 1, 0, 0, 0, 0, 0), focal_cov_index = 5)
  wtd_aug_z.predicted_matrix <- design_matrix(station_wtd_list[[2]]$total, nobs = nobs, 
                                       lambda_cov = c(1, 1, 0, 0), focal_cov_index = 3)
  wtd_aug_z.propSelected_matrix <- design_matrix(station_wtd_list[[2]]$prop_selected, nobs = nobs, 
                                           lambda_cov = c(1, 1, 0, 0), focal_cov_index = 4)
  
  #'  Predict effect of focal covariate across each MCMC iterations
  #'  Code adapted from AHM book 1 code by Mike Meredith (https://github.com/mikemeredith/AHM_code/blob/main/AHM1_ch06/AHM1_06.11.R_)
  predict_parital_effect <- function(mod, newdata, param_names, sims_param3, sims_param4, sims_param5, sims_param6, sims_param7, pred_names) {
    #'  Number of draws in the MCMC chain
    ndraws <- mod$mcmc.info$n.samples
    
    #'  Design matrix
    newcovs <- newdata[[1]]
    
    #'  Unstandardized new focal covariate
    unstandared_cov <- newdata[[2]]
    
    #'  Create objects holding all MCMC iterations for intercept and year effect
    #'  These parameters show up in all models except null model
    sims_beta0 <- mod$sims.list$beta0
    sims_year <- mod$sims.list$b.year[,3]      # b.year[,3] = year effect for 2022
    
    #'  Assemble linear predictors across all iterations of the MCMC chains
    #'  Create and name matrix of each parameter estimate per iteration
    sims.list_matrix <- bind_cols(sims_beta0, sims_year, sims_param3, sims_param4, sims_param5, sims_param6, sims_param7)
    # sims.list_matrix <- bind_cols(sims_param3, sims_param4, sims_param5, sims_param6, sims_param7)
    names(sims.list_matrix) <- param_names
    sims.list_matrix <- as.matrix(sims.list_matrix)
    print(dim(sims.list_matrix))
    
    #' Subsample MCMC iterations to reduce computation time
    set.seed(2024)
    sub.sample.draws <- 3000      
    selection <- sort(sample(1:ndraws, sub.sample.draws))
    
    #' Array to hold predictions for every observation 
    lamPred <- array(NA, dim = c(nobs, sub.sample.draws))
    
    #'  Matrix multiplication: design matrix by predictor coefficients across iterations
    #'  For each subsampled MCMC iteration, multiply each coefficient by its 
    #'  respective new covariate value and sum across estimates, then exponentiate 
    #'  to calculate the predicted lambda. 
    #'  This is done for each row (nobs) in newcovs and repeated for each subsampled 
    #'  MCMC iteration, producing a predicted lambda for each value of the new
    #'  focal covariate per subsampled MCMC iteration.
    for(i in 1:sub.sample.draws) {
      MCMCstep <- selection[i]
      lamPred[,i] <- exp(newcovs %*% sims.list_matrix[MCMCstep,1:ncol(newcovs)])
      #'  If you want to do it by hand and prove matrix multiplication is doing what you think
      # lamPred[,i] <- exp(sims_beta0[MCMCstep] + sims_year[MCMCstep] + (sims_param3[MCMCstep] * newcovs[,3]))
    }
    
    #' Get posterior mean, SD, 2.5% and 97.5% quadrats for every row 
    meanlam <- apply(lamPred, 1, mean)   # Mean across rows to get posterior mean
    sdlam <- apply(lamPred, 1, sd)
    quant2.5lam <- apply(lamPred, 1, quantile, probs = 0.025)
    quant97.5lam <- apply(lamPred, 1, quantile, probs = 0.975)
    print(max(meanlam))                  # Check maximum
    
    #'  List prediction and covariates together for plotting
    predicted_lambda <- bind_cols(meanlam, sdlam, quant2.5lam, quant97.5lam, newcovs, unstandared_cov)
    names(predicted_lambda) <- pred_names
    
    return(predicted_lambda)
  }
  #'  Top Elk July model: intercept, year, *total selected* (3 variables)
  elk_july_lamPred_selected <- predict_parital_effect(elk_july_topmods[[2]], newdata = elk_july_z.selected_matrix,
                                                      param_names = c("V1", "V2", "V3"), sims_param3 = elk_july_topmods[[2]]$sims.list$b.selected,
                                                      sims_param4 = NULL, sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                      pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", "Intercept", "Year", "Total_Selected.z", "Total_Selected"))
  #'  Top Elk Aug model: intercept, year, *total predicted*, *proportion selected* (4 variables)
  elk_aug_lamPred_predicted <- predict_parital_effect(elk_aug_topmods[[2]], newdata = elk_aug_z.predicted_matrix,
                                                      param_names = c("V1", "V2", "V3", "V4"), sims_param3 = elk_aug_topmods[[2]]$sims.list$b.predicted,
                                                      sims_param4 = 0, sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                      pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                      "Intercept", "Year", "Total_Predicted.z", "Proportion_Selected.z", "Total_Predicted"))
  elk_aug_lamPred_propSelected <- predict_parital_effect(elk_aug_topmods[[2]], newdata = elk_aug_z.propSelected_matrix,
                                                         param_names = c("V1", "V2", "V3", "V4"), sims_param3 = 0,
                                                         sims_param4 = elk_aug_topmods[[2]]$sims.list$b.prop.selected, 
                                                         sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                         pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                         "Intercept", "Year", "Total_Predicted.z", "Proportion_Selected.z", "Proportion_Selected"))
  
  #'  Top WTD July model: intercept, year, *total selected*, *proportion selected* (4 columns) #mean Tbio, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 variables)
  # wtd_july_lamPred_cvHQ <- predict_parital_effect(wtd_july_topmods[[1]], newdata = wtd_july_z.cvHQ_matrix,
  #                                                 param_names = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), 
  #                                                 sims_param3 = 0, sims_param4 = 0, 
  #                                                 sims_param5 = wtd_july_topmods[[1]]$sims.list$b.cvHQ, 
  #                                                 sims_param6 = 0, sims_param7 = 0,
  #                                                 pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
  #                                                                "Intercept", "Year", "Mean_TBiomass.z", "Max_HQ.z", "CV_HQ.z", 
  #                                                                "Total_Selected.z", "Proportion_Selected.z", "CV_HQ"))
  wtd_july_lamPred_selected <- predict_parital_effect(wtd_july_topmods[[1]], newdata = wtd_july_z.selected_matrix,
                                                      param_names = c("V1", "V2", "V3", "V4"), sims_param3 = wtd_july_topmods[[1]]$sims.list$b.selected, 
                                                      sims_param4 = 0, sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                      pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                     "Intercept", "Year", "Total_Selected.z", "Proportion_Selected.z", "Total_Selected"))
  wtd_july_lamPred_propSelected <- predict_parital_effect(wtd_july_topmods[[1]], newdata = wtd_july_z.propSelected_matrix,
                                                          param_names = c("V1", "V2", "V3", "V4"), sims_param3 = 0, 
                                                          sims_param4 = wtd_july_topmods[[1]]$sims.list$b.prop.selected, 
                                                          sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                          pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                         "Intercept", "Year", "Total_Selected.z", "Proportion_Selected.z", "Proportion_Selected"))
  
  #'  Top WTD Aug model: intercept, year, *total predicted*, *proportion selected* (4 columns) #*mean Tbio*, max HQ, *cv HQ*, *total selected*, *proportion selected* (7 variables)
  # wtd_aug_lamPred_meanTbio <- predict_parital_effect(wtd_aug_topmods[[1]], newdata = wtd_aug_z.meanTbio_matrix,
  #                                                    param_names = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), 
  #                                                    sims_param3 = wtd_aug_topmods[[1]]$sims.list$b.meanTbio, 
  #                                                    sims_param4 = 0, sims_param5 = 0, sims_param6 = 0, sims_param7 = 0,
  #                                                    pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
  #                                                                   "Intercept", "Year", "Mean_TBiomass.z", "Max_HQ.z", "CV_HQ.z", 
  #                                                                   "Total_Selected.z", "Proportion_Selected.z", "Mean_TBiomass"))
  # wtd_aug_lamPred_cvHQ <- predict_parital_effect(wtd_aug_topmods[[1]], newdata = wtd_aug_z.cvHQ_matrix,
  #                                                param_names = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), 
  #                                                sims_param3 = 0, sims_param4 = 0, 
  #                                                sims_param5 = wtd_aug_topmods[[1]]$sims.list$b.cvHQ, 
  #                                                sims_param6 = 0, sims_param7 = 0,
  #                                                pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
  #                                                               "Intercept", "Year", "Mean_TBiomass.z", "Max_HQ.z", "CV_HQ.z", 
  #                                                               "Total_Selected.z", "Proportion_Selected.z", "CV_HQ"))
  wtd_aug_lamPred_predicted <- predict_parital_effect(wtd_aug_topmods[[1]], newdata = wtd_aug_z.predicted_matrix,
                                                      param_names = c("V1", "V2", "V3", "V4"), sims_param3 = wtd_aug_topmods[[1]]$sims.list$b.predicted,
                                                      sims_param4 = 0, sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                      pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                     "Intercept", "Year", "Total_Predicted.z", "Proportion_Selected.z", "Total_Predicted"))
  wtd_aug_lamPred_propSelected <- predict_parital_effect(wtd_aug_topmods[[1]], newdata = wtd_aug_z.propSelected_matrix,
                                                         param_names = c("V1", "V2", "V3", "V4"), sims_param3 = 0, 
                                                         sims_param4 = wtd_aug_topmods[[1]]$sims.list$b.prop.selected, 
                                                         sims_param5 = NULL, sims_param6 = NULL, sims_param7 = NULL,
                                                         pred_names = c("meanLambda", "sdLambda", "llLambda", "ulLambda", 
                                                                        "Intercept", "Year", "Total_Predicted.z", "Proportion_Selected.z", "Proportion_Selected"))
  #'  SAVE!
  write_csv(elk_july_lamPred_selected, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_elk_July_TotalSelected.csv")
  write_csv(elk_aug_lamPred_predicted, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_elk_Aug_TotalPredicted.csv")
  write_csv(elk_aug_lamPred_propSelected, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_elk_Aug_PropSelected.csv")
  # write_csv(wtd_july_lamPred_cvHQ, "./Outputs/Hilger_RNmodel/Tables/Predicted_lambda_wtd_July_cvHighQuality.csv")
  write_csv(wtd_july_lamPred_selected, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_wtd_July_TotalSelected.csv")
  write_csv(wtd_july_lamPred_propSelected, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_wtd_July_PropSelected.csv")
  # write_csv(wtd_aug_lamPred_meanTbio, "./Outputs/Hilger_RNmodel/Tables/Predicted_lambda_wtd_Aug_MeanTotalBiomass.csv")
  # write_csv(wtd_aug_lamPred_cvHQ, "./Outputs/Hilger_RNmodel/Tables/Predicted_lambda_wtd_Aug_cvHighQuality.csv")
  write_csv(wtd_aug_lamPred_predicted, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_wtd_Aug_TotalPredicted.csv")
  write_csv(wtd_aug_lamPred_propSelected, "./Outputs/Hilger_RNmodel/Tables/Predicted lambda tables/Predicted_lambda_wtd_Aug_PropSelected.csv")
  
  
  #'  -------------------------
  ####  Partial Effects Plots  ####
  #'  -------------------------
  ggplot(elk_july_lamPred_selected, aes(x = Total_Selected, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Total selected forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of total selected species on July elk abundance")
  ggplot(elk_aug_lamPred_predicted, aes(x = Total_Predicted, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Total predicted forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of total predicted species on Aug elk abundance")
  ggplot(elk_aug_lamPred_propSelected, aes(x = Proportion_Selected, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Proportion of selected forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of proportion of selected species on Aug elk abundance")
  # ggplot(wtd_july_lamPred_cvHQ, aes(x = CV_HQ, y = meanLambda)) +
  #   geom_line(lwd = 1) + 
  #   geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
  #   theme_classic() +
  #   xlab("CV high quality forage available") + 
  #   ylab("Predicted lambda (index of relative abundance)") +
  #   ggtitle("Partial effect of CV high quality forage on July white-tailed deer abundance")
  ggplot(wtd_july_lamPred_selected, aes(x = Total_Selected, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Total selected forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of total selected species on July white-tailed deer abundance")
  ggplot(wtd_july_lamPred_propSelected, aes(x = Proportion_Selected, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Proportion of selected forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of proportion of selected species on July white-tailed deer abundance")
  # ggplot(wtd_aug_lamPred_meanTbio, aes(x = Mean_TBiomass, y = meanLambda)) +
  #   geom_line(lwd = 1) + 
  #   geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
  #   theme_classic() +
  #   xlab("Mean total forage biomass available (kg/ha)") + 
  #   ylab("Predicted lambda (index of relative abundance)") +
  #   ggtitle("Partial effect of mean total biomass on August white-tailed deer abundance")
  # ggplot(wtd_aug_lamPred_cvHQ, aes(x = CV_HQ, y = meanLambda)) +
  #   geom_line(lwd = 1) + 
  #   geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
  #   theme_classic() +
  #   xlab("CV high quality forage available") + 
  #   ylab("Predicted lambda (index of relative abundance)") +
  #   ggtitle("Partial effect of CV high quality forage on August white-tailed deer abundance")
  ggplot(wtd_aug_lamPred_predicted, aes(x = Total_Predicted, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Total predicted forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of total predicted species on August white-tailed deer abundance")
  ggplot(wtd_aug_lamPred_propSelected, aes(x = Proportion_Selected, y = meanLambda)) +
    geom_line(lwd = 1) + 
    geom_ribbon(aes(ymin = llLambda, ymax = ulLambda), alpha = 0.3) +
    theme_classic() +
    xlab("Proportion of selected forage species present") + 
    ylab("Predicted lambda (index of relative abundance)") +
    ggtitle("Partial effect of proportion of selected species on August white-tailed deer abundance")
  
  
  
  
  