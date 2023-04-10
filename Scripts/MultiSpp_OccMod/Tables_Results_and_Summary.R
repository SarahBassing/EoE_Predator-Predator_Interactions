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
        Estimate = round(mean, rounddig),
        lower = round(`2.5%`, rounddig),
        upper = round(`97.5%`, rounddig)
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
  occ_coy.bob <- rename_occ_params(out_coy.bob[[1]], intx3 = "N white-tailed deer", intx4 = "N lagomorph", intx5 = "Shannon's H", 
                                   cov2 = "Trail setup", cov3 = "Year 2", cov4 = "Elevation", cov5 = "Forest cover", 
                                   cov6 = "N white-tailed deer", cov7 = "N lagomorph", cov8 = "Shannon's H")
  
  
    renamed_out <- out %>%
    #'  rename interaction intercept first
    mutate(Parameter = ifelse(grepl("betaSpp12\\[1]", Parameter), "Co-occ intercept", Parameter),
           #'  rename psi/p intercepts and trail setup parameters
           Parameter = ifelse(grepl("\\[1]", Parameter), "Intercept", Parameter), 
           Parameter = ifelse(grepl("\\[2]", Parameter), "Trail setup", Parameter),
           #'  rename any covariates on interaction term
           Parameter = ifelse(Parameter == "betaSpp12[3]", "N white-tailed deer", Parameter),
           Parameter = ifelse(Parameter == "betaSpp12[4]", "N lagomorph", Parameter),
           Parameter = ifelse(Parameter == "betaSpp12[5]", "Shannon's H", Parameter),
           Parameter = ifelse(grepl("\\[3]", Parameter), "Year 2", Parameter),
           Parameter = ifelse(grepl("\\[4]", Parameter), "Elevation", Parameter),
           Parameter = ifelse(grepl("\\[5]", Parameter), "Forest cover", Parameter),
           Parameter = ifelse(grepl("\\[6]", Parameter), "N white-tailed deer", Parameter),
           Parameter = ifelse(grepl("\\[7]", Parameter), "N lagomorph", Parameter),
           Parameter = ifelse(grepl("\\[8]", Parameter), "Shannon's H", Parameter))
  return(renamed_out)
  
