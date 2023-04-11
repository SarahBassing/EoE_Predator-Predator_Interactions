  #'  ---------------------------------
  #'  Result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------
  #'  Script to summarize covariate and detection data, as well as results from 
  #'  top models and organize them into result tables.
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(stringr)
  library(tidyverse)
  
    
  #'  ------------------------
  ####  Data summary tables  ####
  #'  ------------------------

  #'  Load covariate and detection data
  load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_img.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  
  
  #####  Summarize detection data  ####
  #'  ------------------------------
  #'  Function to summarize species-specific detection data
  all_detections <- function(dh1, dh2, spp) {
    dh1 <- as.data.frame(dh1[[1]])
    dh2 <- as.data.frame(dh2[[1]])
    #'  Rename columns in each detection history
    newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
    colnames(dh1) <- newcols; colnames(dh2) <- newcols
    #'  Add year to each data set
    dh1$Year <- "2020"
    dh2$Year <- "2021"
    #'  Bind data and summarize detections
    dh <- rbind(dh1, dh2) %>%
      rownames_to_column(., "NewLocationID") %>%
      relocate(Year, .before = occ1) %>%
             #'  Sum detection events across sampling occasions
      mutate(total_dets = select(., occ1:occ11) %>% rowSums(na.rm = TRUE),
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
  DH_bear <- all_detections(DH_eoe20s_predators[[1]], DH_eoe21s_predators[[1]], spp = "Black bear")
  DH_bob <- all_detections(DH_eoe20s_predators[[2]], DH_eoe21s_predators[[2]], spp = "Bobcat")
  DH_coy <- all_detections(DH_eoe20s_predators[[3]], DH_eoe21s_predators[[3]], spp = "Coyote")
  DH_lion <- all_detections(DH_eoe20s_predators[[4]], DH_eoe21s_predators[[4]], spp = "Mountain lion")
  DH_wolf <- all_detections(DH_eoe20s_predators[[5]], DH_eoe21s_predators[[5]], spp = "Wolf")
  
  #'  Final detection history summary table
  DH_summary <- rbind(DH_bear, DH_bob, DH_coy, DH_lion, DH_wolf) %>%
    dplyr::select(-`Number of operating cameras`)
  
  #'  Save!
  write.csv(DH_summary, "./Outputs/Tables/Summary_table_DH.csv")
  
  
  #####  Summarize covariate data  ####
  #'  ------------------------------
  Covariate <- c("Elevation (m)", "Forest cover (%)", "Elk mean RAI", 
                 "Lagomorph mean RAI", "Moose mean RAI", 
                 "Mule deer mean RAI", "White-tailed deer mean RAI", 
                 "Shannon's diveristy index (H)")
  covs <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
    dplyr::select(c(NewLocationID, Elevation__10m2, perc_forest, elk_perday, lagomorphs_perday, 
                    moose_perday,  muledeer_perday, whitetaileddeer_perday, H)) 
  nobs <- nrow(covs)
  cov_means <- covs %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    gather(key = "Variable", value = "Mean") 
  cov_sd <- covs %>% summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))
  cov_se <- cov_sd/sqrt(nobs); cov_se <- gather(cov_se, key = "Variable", value = "SE")
  cov_min <- covs %>% summarise(across(where(is.numeric), ~ min(.x, na.rm = TRUE))) %>%
    gather(key = "Variable", value = "Min")
  cov_max <- covs %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
    gather(key = "Variable", value = "Max")
  cov_summary <- full_join(cov_means, cov_se, by = "Variable") %>%
    full_join(cov_min, by = "Variable") %>%
    full_join(cov_max, by = "Variable") %>%
    cbind(Covariate) %>%
    relocate(Covariate, .before = "Mean") %>%
    dplyr::select(-Variable) %>%
    mutate(Mean = round(Mean, 2), 
           SE = round(SE, 3), 
           Min = round(Min, 2), 
           Max = round(Max, 2))
  
  #'  Table number of cameras per setup and year
  eoe_covs_20s$Year <- "2020"
  eoe_covs_21s$Year <- "2021"
  covs <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
    mutate(Setup = ifelse(grepl("P", NewLocationID), "Predator", "Ungulate")) %>%
    dplyr::select(c(Year, Setup)) #, GMU
  cam_deployment_summary <- as.data.frame(table(covs)) %>%
    arrange(Year)
  colnames(cam_deployment_summary) <- c("Year", "Camera setup", "Operable cameras (n)")
  
  #'  Save covariate summary tables
  write.csv(cov_summary, "./Outputs/Tables/Summary_table_covariates.csv")
  write.csv(cam_deployment_summary, "./Outputs/Tables/Summary_table_camera_deployment.csv")
  
  
  
  #'  --------------------------------
  ####  Model result summary tables  ####
  #'  --------------------------------
  
  #'  Identify top models
  load("./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.RData")
  print(topmodels)
  
  #'  Load top models              #######  MAKE SURE THESE ARE UP-TO-DATE!!!  #######
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2023-04-08.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-09.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_2023-04-07.RData")
  
  
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species info
  rounddig <- 2
  mod_out <- function(mod, spp1, spp2) {
    out <- as.data.frame(mod$summary) %>%
      rownames_to_column(., "Parameter") %>%
      transmute(
        Parameter = Parameter,
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Mean = format(round(mean, rounddig), nsmall = 2),
        Lower_CRI = format(round(`2.5%`, rounddig), nsmall = 2),
        Upper_CRI = format(round(`97.5%`, rounddig), nsmall = 2),
        Mean = as.numeric(Mean),
        Lower_CRI = as.numeric(Lower_CRI),
        Upper_CRI = as.numeric(Upper_CRI)
      ) %>%
      relocate(Parameter, .before = Mean) %>%
      filter(!str_detect(Parameter, "z")) %>%
      filter(Parameter != "deviance") 
    
    #'  Split out results into occupancy, detection, and mean probability dfs
    occ_out <- out %>%
      filter(str_detect(Parameter, "beta")) %>%
      #'  Remove estimates that were fixed to 0 (psix if no interaction)
      filter(Mean != 0)
    #'  Detection results
    det_out <- out %>%
      filter(str_detect(Parameter, "alpha")) %>%
      #'  Remove estimates that were fixed to 0 (px if no interaction)
      filter(Mean != 0)
    #'  Mean occupancy and detection probabilities
    mean_out <- out %>%
      filter(str_detect(Parameter, "mean"))
    #'  List dfs
    split_out <- list(occ_out, det_out, mean_out)
    
    return(split_out)
  }
  out_wolf.bear_top <- mod_out(wolf.bear.hab, "Wolf", "Black bear")
  out_wolf.bear_2nd <- mod_out(wolf.bear.preydiv, "Wolf", "Black bear")
  out_wolf.coy <- mod_out(wolf.coy.hab, "Wolf", "Coyote")
  out_wolf.lion <- mod_out(wolf.lion.null, "Wolf", "Mountain lion")
  out_lion.bear <- mod_out(lion.bear.null, "Mountain lion", "Black bear")
  out_lion.bob <- mod_out(lion.bob.null, "Mountain lion", "Bobcat")
  out_coy.bob <- mod_out(coy.bob.global, "Coyote", "Bobcat")
  
  
  #####  Occupancy results  ####
  #'  -----------------------
  #'  Switch place-holder parameter names with useful ones for occupancy submodel
  rename_occ_params <- function(out, intx3, intx4, intx5, cov2, cov3, cov4, cov5, cov6, cov7, cov8) {
    renamed_out <- out %>%
      #'  Add species names to appropriate parameters
      mutate(Parameter = str_replace(Parameter, "betaSpp12", "Interaction"),
             Parameter = str_replace(Parameter, "betaSpp1", Species1),
             Parameter = str_replace(Parameter, "betaSpp2", Species2),
             #'  Rename psi/p intercepts and trail setup parameters
             Parameter = str_replace(Parameter, "\\[1]", ": Intercept"),
             Parameter = str_replace(Parameter, "\\[2]", paste(":", cov2)),
             # Parameter = ifelse(grepl("\\[1]", Parameter), "Intercept", Parameter), 
             # Parameter = ifelse(grepl("\\[2]", Parameter), "Trail setup", Parameter),
             #'  rename any covariates on interaction term
             Parameter = str_replace(Parameter, "Interaction\\[3]", paste("Interaction:",intx3)),
             Parameter = str_replace(Parameter, "Interaction\\[4]", paste("Interaction:",intx4)),
             Parameter = str_replace(Parameter, "Interaction\\[5]", paste("Interaction:",intx5)),
             # Parameter = ifelse(Parameter == "Interaction[3]", intx3, Parameter),
             # Parameter = ifelse(Parameter == "Interaction[4]", intx4, Parameter),
             # Parameter = ifelse(Parameter == "Interaction[5]", intx5, Parameter),
             Parameter = str_replace(Parameter, "\\[3]", paste(":", cov3)),
             Parameter = str_replace(Parameter, "\\[4]", paste(":", cov4)),
             Parameter = str_replace(Parameter, "\\[5]", paste(":", cov5)),
             Parameter = str_replace(Parameter, "\\[6]", paste(":", cov6)),
             Parameter = str_replace(Parameter, "\\[7]", paste(":", cov7)),
             Parameter = str_replace(Parameter, "\\[8]", paste(":", cov8)))
             # Parameter = ifelse(grepl("\\[3]", Parameter), cov3, Parameter),
             # Parameter = ifelse(grepl("\\[4]", Parameter), cov4, Parameter),
             # Parameter = ifelse(grepl("\\[5]", Parameter), cov5, Parameter),
             # Parameter = ifelse(grepl("\\[6]", Parameter), cov6, Parameter),
             # Parameter = ifelse(grepl("\\[7]", Parameter), cov7, Parameter),
             # Parameter = ifelse(grepl("\\[8]", Parameter), cov8, Parameter))
    return(renamed_out)
  }
  occ_wolf.bear.top <- rename_occ_params(out_wolf.bear_top[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                         cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                         cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = ifelse(Parameter == "Interaction", paste0(Parameter, ": Intercept"), Parameter),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_wolf.bear.2nd <- rename_occ_params(out_wolf.bear_2nd[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                         cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                         cov6 = "Shannon's H", cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = ifelse(Parameter == "Interaction", paste0(Parameter, ": Intercept"), Parameter),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_wolf.coy <- rename_occ_params(out_wolf.coy[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = ifelse(Parameter == "Interaction", paste0(Parameter, ": Intercept"), Parameter),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  occ_wolf.lion <- rename_occ_params(out_wolf.lion[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  occ_lion.bear <- rename_occ_params(out_lion.bear[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                     cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_lion.bob <- rename_occ_params(out_lion.bob[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  occ_coy.bob <- rename_occ_params(out_coy.bob[[1]], intx3 = "N white-tailed deer", intx4 = "N lagomorph", intx5 = "Shannon's H", 
                                   cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                   cov6 = "N white-tailed deer", cov7 = "N lagomorph", cov8 = "Shannon's H") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  
  #'  Combine all occupancy results (long table)
  top_null_results <- rbind(occ_wolf.lion, occ_lion.bear, occ_lion.bob)
  top_non_null_results <- rbind(occ_wolf.bear.2nd, occ_wolf.coy, occ_coy.bob) # NOTE: using 2nd best supported wolf.bear model
  top_occmod_table_long <- rbind(occ_wolf.bear.2nd, occ_wolf.coy, occ_wolf.lion, occ_lion.bear, occ_lion.bob, occ_coy.bob) # NOTE: using 2nd best supported wolf.bear model
  
  #'  Reformat into a wide table
  top_occmod_table_wide <- top_occmod_table_long %>%
    #'  Remove any extra spaces introduced during data formatting
    mutate(Lower_CRI = str_replace_all(Lower_CRI, " ", ""),
           Upper_CRI = str_replace_all(Upper_CRI, " ", ""),
           #'  Combine into single 95% CRI column
           CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")")) %>%
    dplyr::select(-c(Lower_CRI, Upper_CRI)) %>%
    unite(Mean_CRI, Mean, CRI, sep = " ") %>%
    spread(Parameter, Mean_CRI) %>%
    relocate("Species 1: Intercept", .after = "Species2") %>%
    relocate("Species 1: Trail setup", .after = "Species 1: Intercept") %>%
    relocate("Species 1: Year 2", .after = "Species 1: Trail setup") %>%
    relocate("Species 1: Elevation", .after = "Species 1: Year 2") %>%
    relocate("Species 1: Forest cover", .after = "Species 1: Elevation") %>%
    relocate("Species 1: N white-tailed deer", .after = "Species 1: Forest cover") %>%
    relocate("Species 1: N lagomorph", .after = "Species 1: N white-tailed deer") %>%
    relocate("Species 1: Shannon's H", .after = "Species 1: N lagomorph") %>%
    relocate("Species 2: Intercept", .after = "Species 1: Shannon's H") %>%
    relocate("Species 2: Trail setup", .after = "Species 2: Intercept") %>%
    relocate("Species 2: Year 2", .after = "Species 2: Trail setup") %>%
    relocate("Species 2: Elevation", .after = "Species 2: Year 2") %>%
    relocate("Species 2: Forest cover", .after = "Species 2: Elevation") %>%
    relocate("Species 2: N white-tailed deer", .after = "Species 2: Forest cover") %>%
    relocate("Species 2: N lagomorph", .after = "Species 2: N white-tailed deer") %>%
    relocate("Species 2: Shannon's H", .after = "Species 2: N lagomorph") %>%
    relocate("Interaction: Intercept", .after = "Species 2: Shannon's H") %>%
    relocate("Interaction: Trail setup", .after = "Species 2: Intercept") %>%
    relocate("Interaction: N white-tailed deer", .after = "Species 2: Trail setup") %>%
    relocate("Interaction: N lagomorph", .after = "Species 2: N white-tailed deer") %>%
    relocate("Interaction: Shannon's H", .after = "Species 2: N lagomorph")
  
  #'  Save occupancy results
  write.csv(top_null_results, file = paste0("./Outputs/Tables/top_null_results_", Sys.Date(), ".csv"))
  write.csv(top_non_null_results, file = paste0("./Outputs/Tables/top_non_null_results_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_long, file = paste0("./Outputs/Tables/top_occmod_table_long_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_wide, file = paste0("./Outputs/Tables/top_occmod_table_wide_", Sys.Date(), ".csv"))
  
  
  #####  Detection results  ####
  #'  -----------------------
  #'  Switch place-holder parameter names with useful ones for occupancy submodel
  rename_occ_params <- function(out, cov2, cov3) {
    renamed_out <- out %>%
      #'  Add species names to appropriate parameters
      mutate(Parameter = str_replace(Parameter, "alphaSpp12", "Interaction Spp12"),
             Parameter = str_replace(Parameter, "alphaSpp21", "Interaction Spp21"),
             Parameter = str_replace(Parameter, "alphaSpp1", Species1),
             Parameter = str_replace(Parameter, "alphaSpp2", Species2),
             #'  Rename parameters
             Parameter = str_replace(Parameter, "\\[1]", ": Intercept"),
             Parameter = str_replace(Parameter, "\\[2]", paste(":", cov2)),
             Parameter = str_replace(Parameter, "\\[3]", paste(":", cov3)))
    return(renamed_out)
  }
  det_wolf.bear.top <- rename_occ_params(out_wolf.bear_top[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_wolf.bear.2nd <- rename_occ_params(out_wolf.bear_2nd[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_wolf.coy <- rename_occ_params(out_wolf.coy[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  det_wolf.lion <- rename_occ_params(out_wolf.lion[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  det_lion.bear <- rename_occ_params(out_lion.bear[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_lion.bob <- rename_occ_params(out_lion.bob[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  det_coy.bob <- rename_occ_params(out_coy.bob[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  
  top_detmod_table_long <- rbind(det_wolf.bear.2nd, det_wolf.coy, det_wolf.lion, det_lion.bear, det_lion.bob, det_coy.bob) # NOTE: using 2nd best supported wolf.bear model
  
  #'  Reformat into a wide table
  top_detmod_table_wide <- top_detmod_table_long %>%
    #'  Remove any extra spaces introduced during data formatting
    mutate(Lower_CRI = str_replace_all(Lower_CRI, " ", ""),
           Upper_CRI = str_replace_all(Upper_CRI, " ", ""),
           #'  Combine into single 95% CRI column
           CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")")) %>%
    dplyr::select(-c(Lower_CRI, Upper_CRI)) %>%
    unite(Mean_CRI, Mean, CRI, sep = " ") %>%
    spread(Parameter, Mean_CRI) %>%
    relocate("Species 1: Intercept", .after = "Species2") %>%
    relocate("Species 1: Trail setup", .after = "Species 1: Intercept") %>%
    relocate("Species 1: Sampling effort", .after = "Species 1: Trail setup") %>%
    relocate("Species 2: Intercept", .after = "Species 1: Sampling effort") %>%
    relocate("Species 2: Trail setup", .after = "Species 2: Intercept") %>%
    relocate("Species 2: Sampling effort", .after = "Species 2: Trail setup") 
  
  #'  Save detection results
  write.csv(top_detmod_table_long, file = paste0("./Outputs/Tables/top_detmod_table_long_", Sys.Date(), ".csv"))
  write.csv(top_detmod_table_wide, file = paste0("./Outputs/Tables/top_detmod_table_wide_", Sys.Date(), ".csv"))
