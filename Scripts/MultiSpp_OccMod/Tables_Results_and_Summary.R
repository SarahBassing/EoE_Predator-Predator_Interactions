  #'  ---------------------------------
  #'  Result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  April 2023
  #'  ---------------------------------
  #'  Script to summarize results from top models in table format
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(stringr)
  library(tidyverse)
  
    #'  Identify top models
  load("./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.RData")
  print(topmodels)
  
  #'  Load top models
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2023-04-08.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-09.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(.)_p(.)_2023-04-05.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(.)_p(.)_2023-03-31.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(.)_p(.)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_2023-04-07.RData")

  ####  Summary tables  ####
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
        Estimate = format(round(mean, rounddig), nsmall = 2),
        lower = format(round(`2.5%`, rounddig), nsmall = 2),
        upper = format(round(`97.5%`, rounddig), nsmall = 2)
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      filter(!str_detect(Parameter, "z")) %>%
      filter(Parameter != "deviance") 
    
    #'  Split out results into occupancy, detection, and mean probability dfs
    occ_out <- out %>%
      filter(str_detect(Parameter, "beta")) %>%
      #'  Remove estimates that were fixed to 0 (psix if no interaction)
      filter(Estimate != 0)
    #'  Detection results
    det_out <- out %>%
      filter(str_detect(Parameter, "alpha")) %>%
      #'  Remove estimates that were fixed to 0 (px if no interaction)
      filter(Estimate != 0)
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
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Wolf", "Species 1"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 2"))
  occ_lion.bear <- rename_occ_params(out_lion.bear[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                     cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Black bear", "Species 2"))
  occ_lion.bob <- rename_occ_params(out_lion.bob[[1]], intx3 = NA, intx4 = NA, intx5 = NA, 
                                    cov2 = "Year 2", cov3 = NA, cov4 = NA, cov5 = NA, 
                                    cov6 = NA, cov7 = NA, cov8 = NA) %>%
    mutate(Parameter = paste0(Parameter, ": Intercept"),
           Parameter = str_replace(Parameter, "Mountain lion", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  occ_coy.bob <- rename_occ_params(out_coy.bob[[1]], intx3 = "N white-tailed deer", intx4 = "N lagomorph", intx5 = "Shannon's H", 
                                   cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                   cov6 = "N white-tailed deer", cov7 = "N lagomorph", cov8 = "Shannon's H") %>%
    mutate(Parameter = str_replace(Parameter, "Coyote", "Species 1"),
           Parameter = str_replace(Parameter, "Bobcat", "Species 2"))
  
  
  top_null_results <- rbind(occ_wolf.lion, occ_lion.bear, occ_lion.bob)
  top_non_null_results <- rbind(occ_wolf.bear.2nd, occ_wolf.coy, occ_coy.bob) # NOTE: using 2nd best supported wolf.bear model
  top_occmod_table_long <- rbind(occ_wolf.bear.2nd, occ_wolf.coy, occ_wolf.lion, occ_lion.bear, occ_lion.bob, occ_coy.bob) # NOTE: using 2nd best supported wolf.bear model
  
  
  top_occmod_table_wide <- all_occ_results %>%
    mutate(lower = str_replace_all(lower, " ", ""),
           upper = str_replace_all(upper, " ", ""),
           CRI = paste0("(", lower, ", ", upper, ")")) %>%
    dplyr::select(-c(lower, upper)) %>%
    unite(Mean_CRI, Estimate, CRI, sep = " ") %>%
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
  
  
  
  #'  Save
  write.csv(top_null_results, file = paste0("./Outputs/Tables/top_null_results_", Sys.Date(), ".csv"))
  write.csv(top_non_null_results, file = paste0("./Outputs/Tables/top_non_null_results_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_long, file = paste0("./Outputs/Tables/top_occmod_table_long_", Sys.Date(), ".csv"))
  write.csv(top_occmod_table_wide, file = paste0("./Outputs/Tables/top_occmod_table_wide_", Sys.Date(), ".csv"))
  
  
  
  
  
  
