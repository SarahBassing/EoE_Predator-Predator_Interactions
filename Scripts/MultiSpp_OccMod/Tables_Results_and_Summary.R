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
  load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_img_updated_072924.RData")
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
  
  #'  Calculate mean number of detection events per species per year
  (mean_dets <- mean(DH_summary$`Total detection events`))
  (se_dets <- sd(DH_summary$`Total detection events`)/sqrt(nrow(DH_summary)))
  
  #' #'  Save!
  #' write.csv(DH_summary, "./Outputs/Tables/Summary_table_DH.csv")
  
  
  #####  Summarize covariate data  ####
  #'  ------------------------------
  Covariate <- c("Elevation (m)", "Forest cover (%)", "Terrain ruggedness (TRI)", "Elk RAI", "Lagomorph RAI", "Moose RAI", "White-tailed deer RAI") 
                 #"Cattle RAI", "Mule deer RAI", "Shannon's diversity index (H)")
  covs <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
    dplyr::select(c(NewLocationID, Elevation__10m2, perc_forest, TRI, elk_perday, 
                    lagomorphs_perday, moose_perday, whitetaileddeer_perday)) %>%
    filter(NewLocationID != "GMU6_U_23") #livestock_perday, muledeer_perday, H
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
  
  #'  Sampling effort
  #'  Avg number of days each site was operational (generally until camera failure)
  nobs <- nrow(count_eoe20s21s_effort)
  (avg_op_days <- as.data.frame(count_eoe20s21s_effort) %>% 
    `rownames<-`( NULL ) %>%
    replace(is.na(.), 0) %>%
    mutate(sum = rowSums(.)) %>%
    summarise(avg_days = mean(sum),
              sd_days = sd(sum),
              se_days = sd(sum)/sqrt(nobs),
              max_days = max(sum)))
  #'  Avg number of days per sampling occasion camera was operational
  (mean_effort <- mean(as.matrix(count_eoe20s21s_effort), na.rm = TRUE))
  nobs_effort <- dim(count_eoe20s21s_effort)[1] * dim(count_eoe20s21s_effort)[2]
  (se_effort <- sd(as.matrix(count_eoe20s21s_effort), na.rm = TRUE) / sqrt(nobs_effort))
  (min_effort <- min(as.matrix(count_eoe20s21s_effort), na.rm = TRUE))
  (max_effort <- max(as.matrix(count_eoe20s21s_effort), na.rm = TRUE))
  effort <- c(mean_effort, se_effort, min_effort, max_effort) 
  effort <- format(round(effort, 2), nsmall = 2)
  effort <- c("Sampling effort (n days/week)", effort) 

  #'  Append to covariate summary table
  cov_summary <- rbind(cov_summary, effort) %>%
    mutate(Mean = as.numeric(Mean),
           SE = as.numeric(SE),
           Min = as.numeric(Min),
           Max = as.numeric(Max))

  #'  Table number of cameras per setup and year
  eoe_covs_20s$Year <- "2020"
  eoe_covs_21s$Year <- "2021"
  covs <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
    mutate(Setup = ifelse(grepl("P", NewLocationID), "Predator", "Ungulate")) %>%
    dplyr::select(c(Year, Setup)) #, GMU
  cam_deployment_summary <- as.data.frame(table(covs)) %>%
    arrange(Year)
  colnames(cam_deployment_summary) <- c("Year", "Camera setup", "Operable cameras (n)")
  
  #' #'  Save covariate summary tables
  #' write.csv(cov_summary, "./Outputs/Tables/Summary_table_covariates.csv")
  #' write.csv(cam_deployment_summary, "./Outputs/Tables/Summary_table_camera_deployment.csv")
  
  
  
  #'  --------------------------------
  ####  Model result summary tables  ####
  #'  --------------------------------
  
  #'  Identify top models
  load("./Outputs/Tables/DIC_top_models.RData")
  print(topmodels)
  
  #'  Load top models              #######  MAKE SURE THESE ARE UP-TO-DATE!!!  #######
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2024-07-17.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2024-07-21.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2024-08-29.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2024-08-29.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") 
  
  #'  Additional null models to snag mean psi and p from
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_2024-08-27.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(yr)_p(.)_2024-08-27.RData")  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(yr)_p(.)_2024-08-27.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(yr)_p(.)_2024-08-27.RData")  
  
  #'  Co-detecition model based on each top model
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_px(.)_2024-08-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_px(.)_2024-08-07.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_px(.)_2024-08-07.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_px(.)_2024-08-07.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_px(.)_2024-08-07.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_altGoF_2024-08-25.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_altGoF_2024-08-26.RData")
  
  
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
        Mean = round(mean, rounddig),
        Lower_CRI = round(`2.5%`, rounddig),
        Upper_CRI = round(`97.5%`, rounddig),
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
  out_wolf.bear <- mod_out(wolf.bear.hab, "Wolf", "Black bear")
  out_wolf.coy <- mod_out(wolf.coy.hab, "Wolf", "Coyote")
  out_wolf.lion <- mod_out(wolf.lion.null, "Wolf", "Mountain lion")
  out_lion.bear <- mod_out(lion.bear.null, "Mountain lion", "Black bear")
  out_lion.bob <- mod_out(lion.bob.null, "Mountain lion", "Bobcat")
  out_coy.bob <- mod_out(coy.bob.habx, "Coyote", "Bobcat")
  out_bear.coy <- mod_out(bear.coy.habx, "Black bear", "Coyote")
  
  out_wolf.bear_null <- mod_out(wolf.bear.null, "Wolf", "Black bear")
  out_wolf.coy_null <- mod_out(wolf.coy.null, "Wolf", "Coyote")
  out_coy.bob_null <- mod_out(coy.bob.null, "Coyote", "Bobcat")
  out_bear.coy_null <- mod_out(bear.coy.null, "Black bear", "Coyote")
  
  out_wolf.bear.px <- mod_out(wolf.bear.null.px, "Wolf", "Black bear")
  out_wolf.coy.px <- mod_out(wolf.coy.hab.px, "Wolf", "Coyote")
  out_wolf.lion.px <- mod_out(wolf.lion.null.px, "Wolf", "Mountain lion")
  out_lion.bear.px <- mod_out(lion.bear.null.px, "Mountain lion", "Black bear")
  out_lion.bob.px <- mod_out(lion.bob.null.px, "Mountain lion", "Bobcat")
  out_coy.bob.px <- mod_out(coy.bob.habx.px, "Coyote", "Bobcat")
  out_bear.coy.px <- mod_out(bear.coy.habx.px, "Black bear", "Coyote")
  
  
  #####  Occupancy results  ####
  #'  -----------------------
  #'  Switch place-holder parameter names with useful ones for occupancy submodel
  rename_occ_params <- function(out, intx3, intx4, intx5, cov2, cov3, cov4, cov5, cov6) {
    renamed_out <- out %>%
      #'  Add species names to appropriate parameters
      mutate(Lower_CRI = format(Lower_CRI, nsmall = 2),
             Upper_CRI = format(Upper_CRI, nsmall = 2),
             Parameter = str_replace(Parameter, "betaSpp12", "Interaction"),
             Parameter = str_replace(Parameter, "betaSpp1", Species1),
             Parameter = str_replace(Parameter, "betaSpp2", Species2),
             #'  Rename psi/p intercepts and trail setup parameters
             Parameter = str_replace(Parameter, "\\[1]", ": Intercept"),
             Parameter = str_replace(Parameter, "\\[2]", paste(":", cov2)),
             #'  rename any covariates on interaction term
             Parameter = str_replace(Parameter, "Interaction\\[3]", paste("Interaction:",intx3)),
             Parameter = str_replace(Parameter, "Interaction\\[4]", paste("Interaction:",intx4)),
             Parameter = str_replace(Parameter, "Interaction\\[5]", paste("Interaction:",intx5)),
             Parameter = str_replace(Parameter, "\\[3]", paste(":", cov3)),
             Parameter = str_replace(Parameter, "\\[4]", paste(":", cov4)),
             Parameter = str_replace(Parameter, "\\[5]", paste(":", cov5)),
             Parameter = str_replace(Parameter, "\\[6]", paste(":", cov6)))
    return(renamed_out)
  }
  occ_wolf.bear <- rename_occ_params(out_wolf.bear[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                     cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Forest cover", cov5 = "Elevation",  cov6 = "TRI") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_wolf.coy <- rename_occ_params(out_wolf.coy[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Forest cover", cov5 = "Elevation",  cov6 = "TRI") %>%
    mutate(Parameter = ifelse(Parameter == "Interaction", paste0(Parameter, ": Intercept"), Parameter),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  occ_wolf.lion <- rename_occ_params(out_wolf.lion[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, cov6 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  occ_lion.bear <- rename_occ_params(out_lion.bear[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                     cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, cov6 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_lion.bob <- rename_occ_params(out_lion.bob[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, cov6 = NA) %>%
    mutate(Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  occ_coy.bob <- rename_occ_params(out_coy.bob[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                   cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Forest cover", cov5 = "Elevation", cov6 = "TRI") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  occ_bear.coy <- rename_occ_params(out_bear.coy[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                   cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Forest cover", cov5 = "Elevation", cov6 = "TRI") %>%
    mutate(Parameter = str_replace(Parameter, "Black bear", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  
  #'  Combine all occupancy results (long table)
  top_null_results <- rbind(occ_wolf.lion, occ_lion.bear, occ_lion.bob)
  top_non_null_results <- rbind(occ_wolf.bear, occ_wolf.coy, occ_bear.coy, occ_coy.bob) 
  top_occmod_table_long <- rbind(occ_wolf.bear, occ_wolf.lion, occ_wolf.coy, occ_lion.bear, occ_lion.bob, occ_bear.coy, occ_coy.bob) 
  
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
    relocate("Species 1: Forest cover", .after = "Species 1: Year 2") %>%
    relocate("Species 1: Elevation", .after = "Species 1: Forest cover") %>%
    relocate("Species 1: TRI", .after = "Species 1: Elevation") %>%
    relocate("Species 2: Intercept", .after = "Species 1: TRI") %>%
    relocate("Species 2: Trail setup", .after = "Species 2: Intercept") %>%
    relocate("Species 2: Year 2", .after = "Species 2: Trail setup") %>%
    relocate("Species 2: Forest cover", .after = "Species 2: Year 2") %>%
    relocate("Species 2: Elevation", .after = "Species 2: Forest cover") %>%
    relocate("Species 2: TRI", .after = "Species 2: Elevation") %>%
    relocate("Interaction", .after = "Species 2: TRI") %>%
    rename("Interaction: Intercept" = "Interaction")
  
  #'  Alternate approach to tabling results
  spp1_out <- filter(top_occmod_table_long, grepl("Species 1", Parameter)) %>%
    #'  Drop species identifier in parameter name and reformat mean & 95% CRI
    mutate(Parameter = gsub(".*:", "", Parameter),
           #'  Combine lower & upper 95% CRI into single column
           CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")"),
           #'  Make sure there are two decimals for each value
           Mean = format(Mean, nsmall = 2),
           #'  Combine mean and 95% CRI into single column
           Mean_CRI = paste(Mean, CRI),
           Mean_CRI = format(Mean_CRI, nsmall = 2)) %>%
    dplyr::select(-c(Lower_CRI, Upper_CRI, Mean, CRI)) %>%
    #'  Add species identifier to mean & CRI columns
    rename(Mean_CRI_Species1 = Mean_CRI)
  spp2_out <- filter(top_occmod_table_long, grepl("Species 2", Parameter)) %>%
    mutate(Parameter = gsub(".*:", "", Parameter),
           CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")"),
           Mean = format(Mean, nsmall = 2),
           Mean_CRI = paste(Mean, CRI),
           Mean_CRI = format(Mean_CRI, nsmall = 2)) %>%
    dplyr::select(-c(Lower_CRI, Upper_CRI, Mean, CRI)) %>%
    rename(Mean_CRI_Species2 = Mean_CRI)
  interaction_out <- filter(top_occmod_table_long, grepl("Interaction", Parameter)) %>%
    mutate(Parameter = gsub(".*:", "", Parameter),
           CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")"),
           Mean = format(Mean, nsmall = 2),
           Mean_CRI = paste(Mean, CRI),
           Mean_CRI = format(Mean_CRI, nsmall = 2)) %>%
    dplyr::select(-c(Lower_CRI, Upper_CRI, Mean, CRI)) %>%
    rename(Mean_CRI_Interaction = Mean_CRI)
  #'  Join species 1 & 2 estimates together
  top_occmod_table <- full_join(spp1_out, spp2_out, by = c("Species1", "Species2", "Parameter")) %>%
    full_join(interaction_out, by = c("Species1", "Species2", "Parameter"))
  colnames(top_occmod_table) <- c("Species 1", "Species 2", "Parameter", "Species 1 Mean (95% CRI)", 
                                  "Species 2 Mean (95% CRI)", "Interaction Mean (95% CRI)")
  
  #'  Save occupancy results
  write.csv(top_null_results, file = paste0("./Outputs/Tables/top_null_results_", Sys.Date(), ".csv"))
  write.csv(top_non_null_results, file = paste0("./Outputs/Tables/top_non_null_results_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_long, file = paste0("./Outputs/Tables/top_occmod_table_long_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_wide, file = paste0("./Outputs/Tables/top_occmod_table_wide_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table, file = paste0("./Outputs/Tables/top_occmod_table_", Sys.Date(), ".csv"))
  
 
  #####  Detection results  ####
  #'  -----------------------
  #'  Switch place-holder parameter names with useful ones for occupancy submodel
  rename_det_params <- function(out, cov2, cov3) {
    renamed_out <- out %>%
      #'  Add species names to appropriate parameters
      mutate(Lower_CRI = format(Lower_CRI, nsmall = 2),
             Upper_CRI = format(Upper_CRI, nsmall = 2),
             Parameter = str_replace(Parameter, "alphaSpp12", "Interaction Spp12"),
             Parameter = str_replace(Parameter, "alphaSpp21", "Interaction Spp21"),
             Parameter = str_replace(Parameter, "alphaSpp1", Species1),
             Parameter = str_replace(Parameter, "alphaSpp2", Species2),
             #'  Rename parameters
             Parameter = str_replace(Parameter, "\\[1]", ": Intercept"),
             Parameter = str_replace(Parameter, "\\[2]", paste(":", cov2)),
             Parameter = str_replace(Parameter, "\\[3]", paste(":", cov3)))
    return(renamed_out)
  }
  #'  Top model detection results
  det_wolf.bear <- rename_det_params(out_wolf.bear[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(#Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_wolf.coy <- rename_det_params(out_wolf.coy[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  det_wolf.lion <- rename_det_params(out_wolf.lion[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  det_lion.bear <- rename_det_params(out_lion.bear[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_lion.bob <- rename_det_params(out_lion.bob[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  det_coy.bob <- rename_det_params(out_coy.bob[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  det_bear.coy <- rename_det_params(out_bear.coy[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Black bear", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  
  
  top_detmod_table_long <- rbind(det_wolf.bear, det_wolf.coy, det_wolf.lion, det_lion.bear, det_lion.bob, det_bear.coy, det_coy.bob) 
  
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

  
  #'  Top model with co-detection results
  det_wolf.bear.px <- rename_det_params(out_wolf.bear.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_wolf.coy.px <- rename_det_params(out_wolf.coy.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  det_wolf.lion.px <- rename_det_params(out_wolf.lion.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  det_lion.bear.px <- rename_det_params(out_lion.bear.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  det_lion.bob.px <- rename_det_params(out_lion.bob.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  det_coy.bob.px <- rename_det_params(out_coy.bob.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  det_bear.coy.px <- rename_det_params(out_bear.coy.px[[2]], cov2 = "Trail setup", cov3 = "Sampling effort") %>%
    mutate(Parameter = str_replace(Parameter, "Black bear", "Species 1"),
           Parameter = str_replace(Parameter, "Coyote", "Species 2"))
  
  co_detmod_table_long <- rbind(det_wolf.bear.px, det_wolf.coy.px, det_wolf.lion.px, det_lion.bear.px, det_lion.bob.px, det_coy.bob.px) %>%
    mutate(Parameter = ifelse(Parameter == "Interaction Spp12", "Interaction: Intercept", Parameter),
           Parameter = ifelse(Parameter == "Interaction Spp12: Intercept", "Interaction: Intercept", Parameter))
  
  #'  Reformat into a wide table
  co_detmod_table_wide <- co_detmod_table_long %>%
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
    relocate("Species 2: Sampling effort", .after = "Species 2: Trail setup") %>%
    relocate("Interaction: Intercept", .after = "Species 2: Sampling effort")
  
  #'  Save detection results
  write.csv(top_detmod_table_long, file = paste0("./Outputs/Tables/top_detmod_table_long_", Sys.Date(), ".csv"))
  write.csv(top_detmod_table_wide, file = paste0("./Outputs/Tables/top_detmod_table_wide_", Sys.Date(), ".csv"))
  write.csv(co_detmod_table_long, file = paste0("./Outputs/Tables/co_detmod_table_long_", Sys.Date(), ".csv"))
  write.csv(co_detmod_table_wide, file = paste0("./Outputs/Tables/co_detmod_table_wide_", Sys.Date(), ".csv"))

  
  #####  Mean occupancy & detection probability  ####
  #'  --------------------------------------------
  #'  Switch place-holder parameter names with useful ones for occupancy submodel
  rename_mean_psi_p <- function(out) {
    renamed_out <- out %>%
      #'  Add species names to appropriate parameters
      mutate(Mean = format(Mean, nsmall = 2),
             Lower_CRI = format(Lower_CRI, nsmall = 2),
             Upper_CRI = format(Upper_CRI, nsmall = 2),
             Mean = str_replace_all(Mean, " ", ""),
             Lower_CRI = str_replace_all(Lower_CRI, " ", ""),
             Upper_CRI = str_replace_all(Upper_CRI, " ", ""),
             #'  Combine into single 95% CRI column
             CRI = paste0("(", Lower_CRI, ", ", Upper_CRI, ")"),
             Mean = paste(Mean, CRI),
             Species = ifelse(Parameter == "mean.psiSpp1", Species1, Species2),
             Species = ifelse(Parameter == "mean.pSpp1", Species1, Species),
             Parameter = str_replace(Parameter, "mean.psiSpp1", "Mean occupancy"),
             Parameter = str_replace(Parameter, "mean.psiSpp2", "Mean occupancy"),
             Parameter = str_replace(Parameter, "mean.pSpp1", "Mean detection"),
             Parameter = str_replace(Parameter, "mean.pSpp2", "Mean detection"),
             Predator_pair = paste(Species1, "-", Species2)) %>%
      dplyr::select(-c(Species1, Species2, Lower_CRI, Upper_CRI, CRI)) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Predator_pair, .before = Species)
    return(renamed_out)
  }
  # mean_wolf.coy <- rename_mean_psi_p(out_wolf.coy[[3]])
  # mean_coy.bob <- rename_mean_psi_p(out_coy.bob[[3]])
  # mean_bear.coy <- rename_mean_psi_p(out_bear.coy[[3]])
  # mean_wolf.bear <- rename_mean_psi_p(out_wolf.bear[[3]])
  mean_wolf.lion <- rename_mean_psi_p(out_wolf.lion[[3]])
  mean_lion.bear <- rename_mean_psi_p(out_lion.bear[[3]])
  mean_lion.bob <- rename_mean_psi_p(out_lion.bob[[3]])
  
  mean_wolf.bear <- rename_mean_psi_p(out_wolf.bear_null[[3]])
  mean_wolf.coy <- rename_mean_psi_p(out_wolf.coy_null[[3]])
  mean_coy.bob <- rename_mean_psi_p(out_coy.bob_null[[3]])
  mean_bear.coy <- rename_mean_psi_p(out_bear.coy_null[[3]])
  
  #'  Merge all results  
  mean_occ_det <- rbind(mean_wolf.coy, mean_bear.coy, mean_coy.bob, mean_wolf.bear, mean_wolf.lion, mean_lion.bear, mean_lion.bob)
  #'  Split by occupancy vs detection probability
  mean_occ <- filter(mean_occ_det, Parameter == "Mean occupancy") %>%
    rename("Mean Occupancy (95% CRI)" = "Mean") %>%
    dplyr::select(-Parameter)
  mean_det <- filter(mean_occ_det, Parameter == "Mean detection") %>%
    rename("Mean Detection (95% CRI)" = "Mean") %>%
    dplyr::select(-Parameter)
  #'  Re-combine occupancy and detection estimates so there is 1 row per species
  mean_occ_det <- full_join(mean_occ, mean_det, by = c("Predator_pair", "Species")) %>%
    arrange(Species) %>%
    group_by(Species) %>%
    slice(1L) %>%
    ungroup() %>%
    dplyr::select(-Predator_pair)
  ### NOTE: Reports multiple occ and det estimates per species given multiple models with the same species
  ### NOTE: Estimates differ depending on if they come from null model vs model with covariates
  ### NOTE: mean.psi and mean.p for covariate models are occ/det at ungulate cameras (which is LOW compared to predator cameras)
  ### NOTE: mean.psi and mean.p for null models are true mean occ/det (averaged across ungulate & predator cameras)
  ### NOTE: but don't account for any heterogeneity in occ/det. 
  ### NOTE: Currently using estimates from null models since this is unaffected by camera placement
  
  write.csv(mean_occ_det, "./Outputs/Tables/Summary_mean_psi&p.csv")
  