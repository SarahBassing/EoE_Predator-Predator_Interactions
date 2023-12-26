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
  
  #'  Function to format model results table
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
        `95% CRI` = paste0(Lower_CRI, " - ", Upper_CRI),
        overlap0 = ifelse(overlap0 == 0, FALSE, TRUE)) %>%
      dplyr::select(-c(Lower_CRI, Upper_CRI)) %>%
      #'  Give parameters more meaningful names
      mutate(
        Parameter = ifelse(Parameter == "beta0", "Intercept [GMU10A]", Parameter),
        Parameter = ifelse(Parameter == "beta1", "Percent forest", Parameter),
        Parameter = ifelse(Parameter == "beta2", "Elevation", Parameter),
        Parameter = ifelse(Parameter == "beta3", "Elevation^2", Parameter),
        Parameter = ifelse(Parameter == "beta4[2]", "GMU6", Parameter),              #beta4[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "beta4[3]", "GMU1", Parameter),
        Parameter = ifelse(Parameter == "alpha0", "Intercept [Random]", Parameter),
        Parameter = ifelse(Parameter == "alpha2[2]", "Setup [Road]", Parameter),     #alpha2[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "mu.lambda", "mean lambda", Parameter),
        Parameter = ifelse(Parameter == "mu.r", "mean r", Parameter),
        Parameter = ifelse(Parameter == "mean.psi", "mean Pr(occupancy)", Parameter),
        overlap0 = ifelse(str_detect(Parameter, "Intercept "), NA, overlap0))  %>%
      #'  Save only coefficients and mean derived parameters
      filter(Parameter == "Intercept [GMU10A]" | Parameter == "Percent forest" | Parameter == "Elevation" |
               Parameter == "Elevation^2" | Parameter == "GMU6" | Parameter == "GMU1" | Parameter == "Intercept [Random]" |
               Parameter == "Setup [Road]" | Parameter == "mean lambda" | Parameter == "mean r" | Parameter == "mean Pr(occupancy)")
    
    print(out)
    return(out)
  }
  bear_tbls <- mapply(mod_out_summary_table, mod = bear_list, spp = "Black bear", yr = yr_list, SIMPLIFY = FALSE)
  bob_tbls <- mapply(mod_out_summary_table, mod = bob_list, spp = "Bobcat", yr = yr_list, SIMPLIFY = FALSE)
  coy_tbls <- mapply(mod_out_summary_table, mod = coy_list, spp = "Coyote", yr = yr_list, SIMPLIFY = FALSE)
  lion_tbls <- mapply(mod_out_summary_table, mod = lion_list, spp = "Mountain lion", yr = yr_list, SIMPLIFY = FALSE)
  wolf_tbls <- mapply(mod_out_summary_table, mod = wolf_list, spp = "Wolf", yr = yr_list, SIMPLIFY = FALSE)
  elk_tbls <- mapply(mod_out_summary_table, mod = elk_list, spp = "Elk", yr = yr_list, SIMPLIFY = FALSE)
  moose_tbls <- mapply(mod_out_summary_table, mod = moose_list, spp = "Moose", yr = yr_list, SIMPLIFY = FALSE)
  wtd_tbls <- mapply(mod_out_summary_table, mod = wtd_list, spp = "White-tailed deer", yr = yr_list, SIMPLIFY = FALSE)
  lago_tbls <- mapply(mod_out_summary_table, mod = lago_list, spp = "Lagomorphs", yr = yr_list, SIMPLIFY = FALSE)
  
  #'  Unlist species-specific result tables
  bear_tbl <- do.call("rbind", bear_tbls)
  bob_tbl <- do.call("rbind", bob_tbls)
  coy_tbl <- do.call("rbind", coy_tbls)
  lion_tbl <- do.call("rbind", lion_tbls)
  wolf_tbl <- do.call("rbind", wolf_tbls)
  elk_tbl <- do.call("rbind", elk_tbls)
  moose_tbl <- do.call("rbind", moose_tbls)
  wtd_tbl <- do.call("rbind", wtd_tbls)
  lago_tbl <- do.call("rbind", lago_tbls)
  
  RN_result_tbl <- rbind(bear_tbl, bob_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl, lago_tbl)
  RN_result_tbl_list <- list(bear_tbl, bob_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl, lago_tbl)
  
  #'  Save
  write_csv(RN_result_tbl, file = "./Outputs/Relative_Abundance/RN_model/RN_result_tbl.csv")
  save(RN_result_tbl_list, file = "./Outputs/Relative_Abundance/RN_model/RN_result_tbl_list.R")
  