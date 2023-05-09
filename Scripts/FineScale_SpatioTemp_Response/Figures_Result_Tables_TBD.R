  #'  ---------------------------------
  #'  TBD figures and result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ---------------------------------
  #'  Script to create figures and result tables based on results from analyses 
  #'  that estimated the effect of recent competitor presence and prey availability 
  #'  on the time between detections of predators.
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(ggplot2)
  library(stringr)
  library(tidyverse)
  library(khroma)
  library(patchwork)
  
  #'  Load top models
  load("./Outputs/Time_btwn_Detections/tbd.comp.bear_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_X_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.lion_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.wolf_preyRAI.RData")
  
  #'  Snag and rename coefficents from each model
  coefs <- function(mod_out, spp, prey1, prey2, prey3, comp1, comp2, comp3, comp4) {
    Species <- spp
    Estimate <- round(unlist(mod_out$mean), 2)
    lci <- round(unlist(mod_out$q2.5), 2)
    uci <- round(unlist(mod_out$q97.5), 2)
    CI <- paste(" ", lci, "-", uci) # need that extra space in front b/c excel thinks this is an equation otherwise
    out <- as.data.frame(cbind(Species, Estimate, CI, lci, uci))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "95% CI", "lci", "uci")
    renamed_out <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta.prey1", paste("Prey RAI:", prey1), Parameter),
             Parameter = ifelse(Parameter == "beta.prey2", paste("Prey RAI:", prey2), Parameter),
             Parameter = ifelse(Parameter == "beta.prey3", paste("Prey RAI:", prey3), Parameter),
             Parameter = ifelse(Parameter == "beta.competitor1", paste("Competitor:", comp1), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor2", paste("Competitor:", comp2), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor3", paste("Competitor:", comp3), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor4", paste("Competitor:", comp4), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd1", paste("Mean TBD:", comp1), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd2", paste("Mean TBD:", comp2), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd3", paste("Mean TBD:", comp3), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd4", paste("Mean TBD:", comp4), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd1", paste0("Competitor:Prey (", comp1, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd2", paste0("Competitor:Prey (", comp2, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd3", paste0("Competitor:Prey (", comp3, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd4", paste0("Competitor:Prey (", comp4, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago1", paste0("Competitor:Prey (", comp1, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago2", paste0("Competitor:Prey (", comp2, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago3", paste0("Competitor:Prey (", comp3, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago4", paste0("Competitor:Prey (", comp4, " x ", prey2, ")"), Parameter),
             Parameter = gsub("spp.tbd.", "tbd ", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter)) %>%
      filter(Estimate != 0) %>%
      mutate(Estimate = as.numeric(Estimate),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed_out)
  }
  tbd.bear.out <- coefs(tbd.bear.preyabund, spp = "Black bear", prey1 = "Elk", prey2 = "White-tailed deer")
  tbd.bob.out <- coefs(tbd.bob.compID.preyabund, spp = "Bobcat", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.coy.out <- coefs(tbd.coy.compIDxpreyabund, spp = "Coyote", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Black bear", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.lion.out <- coefs(tbd.lion.preyabund, spp = "Mountain lion", prey1 = "Elk", prey2 = "White-tailed deer")
  tbd.wolf.out <- coefs(tbd.wolf.preyabund, spp = "Wolf", prey1 = "Elk", prey2 = "Moose", prey3 = "White-tailed deer")
  
  #'  Pull out just coefficient estimates
  bear.coefs <- tbd.bear.out[1:4,1:4]
  bob.coefs <- tbd.bob.out[1:11,1:4]
  coy.coefs <- tbd.coy.out[1:17,1:4]
  lion.coefs <- tbd.lion.out[1:4,1:4]
  wolf.coefs <- tbd.wolf.out[1:5,1:4]
  
  #'  Pull out mean TBD estimates
  bear.mean.tbd
  bob.mean.tbd
  coy.mean.tbd
  lion.mean.tbd
  wolf.mean.tbd
  
  #'  Pull out predicted TBD values
  bear.tbd.predictions <- tbd.bear.out[]
  bob.tbd.predictions <- tbd.bob.out[]
  coy.tbd.predictions <- tbd.coy.out[]
  lion.tbd.predictions <- tbd.lion.out[]
  wolf.tbd.predictions <- tbd.wolf.out[]
  
  tbd.coefs <- rbind(bear.coefs, bob.coefs, coy.coefs, lion.coefs, wolf.coefs)
  # write.csv(tbd.coefs, "./Outputs/Tables/TBD_coefficient_estimates_allSpp.csv")
  