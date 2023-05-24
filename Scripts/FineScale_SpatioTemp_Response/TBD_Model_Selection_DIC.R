  #'  ---------------------------------
  #'  TBD model selection & result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ---------------------------------
  #'  Script to identify most supported model using DIC and format result tables
  #'  for analyses that estimated the time between detections of competitors.
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(AICcmodavg)
  library(tidyverse)
  
  #'  -------------------------------
  ####  COMPETITOR - FOCAL PREDATOR  ####
  #'  -------------------------------
  #'  Model selection for models estimating time between detections of two different predators
  
  #####  Black bear models  ####
  load("./Outputs/Time_btwn_Detections/tbd.bear_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_sppID_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_sppID_X_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_sppID_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bear_sppID_X_preyRAI.RData") # Converged poorly; over-parameterized
  load("./Outputs/Time_btwn_Detections/tbd.bear_global.RData")          # Converged poorly; over-parameterized
  #'  List for model selection
  bear_tbd_list <- list(tbd.bear.null, tbd.bear.sppID, tbd.bear.div, tbd.bear.preyabund, tbd.bear.sppID.div, tbd.bear.sppIDxdiv, tbd.bear.sppID.preyabund) #tbd.bear.sppIDxpreyabund, tbd.bear.global
  bear_tbd_name <- c("tbd.bear.null", "tbd.bear.sppID", "tbd.bear.div", "tbd.bear.preyabund", "tbd.bear.sppID.div", "tbd.bear.sppIDxdiv", "tbd.bear.sppID.preyabund") #"tbd.bear.sppIDxpreyabund", "tbd.bear.global"
  
  #####  Bobcat models  ####
  load("./Outputs/Time_btwn_Detections/tbd.bob_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_sppID_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_sppID_X_preydiv.RData")  # Converged poorly; over-parameterized
  load("./Outputs/Time_btwn_Detections/tbd.bob_sppID_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_sppID_X_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.bob_global.RData")           # Converged poorly; over-parameterized
  #'  List for model selection
  bob_tbd_list <- list(tbd.bob.null, tbd.bob.sppID, tbd.bob.div, tbd.bob.preyabund, tbd.bob.sppID.div, tbd.bob.sppID.preyabund, tbd.bob.sppIDxpreyabund) # tbd.bob.sppIDxdiv, tbd.bob.global 
  bob_tbd_name <- c("tbd.bob.null", "tbd.bob.sppID", "tbd.bob.div", "tbd.bob.preyabund", "tbd.bob.sppID.div", "tbd.bob.sppID.preyabund", "tbd.bob.sppIDxpreyabund") # "tbd.bob.sppIDxdiv", "tbd.bob.global" 
    
  #####  Coyote models  ####
  load("./Outputs/Time_btwn_Detections/tbd.coy_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_sppID_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_sppID_X_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_sppID_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_sppID_X_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.coy_global.RData")
  #'  List for model selection
  coy_tbd_list <- list(tbd.coy.null, tbd.coy.sppID, tbd.coy.div, tbd.coy.preyabund, tbd.coy.sppID.div, tbd.coy.sppIDxdiv, tbd.coy.sppID.preyabund, tbd.coy.sppIDxpreyabund, tbd.coy.global) 
  coy_tbd_name <- c("tbd.coy.null", "tbd.coy.sppID", "tbd.coy.div", "tbd.coy.preyabund", "tbd.coy.sppID.div", "tbd.coy.sppIDxdiv", "tbd.coy.sppID.preyabund", "tbd.coy.sppIDxpreyabund", "tbd.coy.global") 
  
  #####  Mountain lion models  ####
  load("./Outputs/Time_btwn_Detections/tbd.lion_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.lion_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.lion_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.lion_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.lion_sppID_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.lion_sppID_X_preydiv.RData") # Converged poorly; over-parameterized
  load("./Outputs/Time_btwn_Detections/tbd.lion_sppID_preyRAI.RData")   
  load("./Outputs/Time_btwn_Detections/tbd.lion_sppID_X_preyRAI.RData") # Converged poorly; over-parameterized
  load("./Outputs/Time_btwn_Detections/tbd.lion_global.RData")          # Converged poorly; over-parameterized
  #'  List for model selection
  lion_tbd_list <- list(tbd.lion.null, tbd.lion.sppID, tbd.lion.div, tbd.lion.preyabund, tbd.lion.sppID.div, tbd.lion.sppID.preyabund) #tbd.lion.sppIDxdiv, tbd.lion.sppIDxpreyabund, tbd.lion.global) 
  lion_tbd_name <- c("tbd.lion.null", "tbd.lion.sppID", "tbd.lion.div", "tbd.lion.preyabund", "tbd.lion.sppID.div", "tbd.lion.sppID.preyabund") #"tbd.lion.sppIDxdiv", "tbd.lion.sppIDxpreyabund", "tbd.lion.global") 
    
  #####  Wolf models  ####
  load("./Outputs/Time_btwn_Detections/tbd.wolf_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_sppID_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_sppID_X_preydiv.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_sppID_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.wolf_sppID_X_preyRAI.RData") # Converged poorly; over-parameterized
  load("./Outputs/Time_btwn_Detections/tbd.wolf_global.RData")          # Converged poorly; over-parameterized
  #'  List for model selection
  wolf_tbd_list <- list(tbd.wolf.null, tbd.wolf.sppID, tbd.wolf.div, tbd.wolf.preyabund, tbd.wolf.sppID.div, tbd.wolf.sppIDxdiv, tbd.wolf.sppID.preyabund) #, tbd.wolf.sppIDxpreyabund, tbd.wolf.global) 
  wolf_tbd_name <- c("tbd.wolf.null", "tbd.wolf.sppID", "tbd.wolf.div", "tbd.wolf.preyabund", "tbd.wolf.sppID.div", "tbd.wolf.sppIDxdiv", "tbd.wolf.sppID.preyabund") #, "tbd.wolf.sppIDxpreyabund", "tbd.wolf.global") 
    
  #####  Model selection using DIC  ####
  (topmod_beartbd <- dictab(cand.set = bear_tbd_list, modnames = bear_tbd_name, sort = TRUE)) 
  (topmod_bobtbd <- dictab(cand.set = bob_tbd_list, modnames = bob_tbd_name, sort = TRUE)) 
  (topmod_coytbd <- dictab(cand.set = coy_tbd_list, modnames = coy_tbd_name, sort = TRUE)) 
  (topmod_liontbd <- dictab(cand.set = lion_tbd_list, modnames = lion_tbd_name, sort = TRUE)) 
  (topmod_wolftbd <- dictab(cand.set = wolf_tbd_list, modnames = wolf_tbd_name, sort = TRUE)) 
  
  #'  Best supported model per species-pair
  (topmodels <- rbind(topmod_beartbd[1,], topmod_bobtbd[1,], topmod_coytbd[1,], topmod_liontbd[1,], topmod_wolftbd[1,]))
  
  
  
  
  
  
  
  
  
  
  
  #'  Full table of models ranked by DIC for all species-pairs
  model_list_DIC <- rbind(topmod_beartbd, topmod_bobtbd, topmod_coytbd, topmod_liontbd, topmod_wolftbd) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    #'  Split model name based on placement of multiple periods
    #'  https://stackoverflow.com/questions/26265400/use-regex-in-r-to-retrieve-string-before-second-occurence-of-a-period
    mutate(Species = sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species = str_replace(Species, "tbd.", ""),
           Species = str_replace(Species, "bear", "Black bear"),
           Species = str_replace(Species, "lion", "Mountain lion"),
           Species = str_replace(Species, "coy", "Coyote"),
           Species = str_replace(Species, "bob", "Bobcat"),
           Species = str_replace(Species, "wolf", "Wolf"),
           Model_name = sub("(^[^.]+[.][^.]+)(.+$)", "\\2", Modnames),
           Model_name = str_replace(Model_name, ".", ""),
           Model = gsub("null", "Model 1", Model_name),           
           Model = gsub("compIDxdiv", "Model 6", Model),
           Model = gsub("compIDxpreyabund", "Model 8", Model),
           Model = gsub("compID.div", "Model 5", Model),
           Model = gsub("compID.preyabund", "Model 7", Model),
           Model = gsub("compID", "Model 2", Model),
           Model = gsub("div", "Model 3", Model),
           Model = gsub("preyabund", "Model 4", Model),
           Model = gsub("global", "Model 9", Model),
           Model_name = gsub("null", "Null", Model_name),
           Model_name = gsub("compIDxdiv", "Competitor ID * prey diversity", Model_name),
           Model_name = gsub("compIDxpreyabund", "Competitor ID * prey abundance", Model_name),
           Model_name = gsub("compID.div", "Competitor ID + prey diversity", Model_name),
           Model_name = gsub("compID.preyabund", "Competitor ID + prey abundance", Model_name),
           Model_name = gsub("compID", "Competitor ID", Model_name),
           Model_name = gsub("div", "Prey diversity", Model_name),
           Model_name = gsub("preyabund", "Prey abundance", Model_name),
           Model_name = gsub("global", "Global", Model_name)) %>%
    relocate(Species, .before = DIC) %>%
    relocate(Model, .after = Species) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_list_DIC) <- c("Species", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/Tables/DIC_TBD_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/Tables/DIC_TBD_model_selection_results.csv")
  save(topmodels, file = "./Outputs/Tables/DIC_TBD_top_models.RData")
  save(model_list_DIC, file = "./Outputs/Tables/DIC_TBD_model_selection_results.RData")
  
  
  #'  Next stop, create figures to visualize these results
  
  
  
