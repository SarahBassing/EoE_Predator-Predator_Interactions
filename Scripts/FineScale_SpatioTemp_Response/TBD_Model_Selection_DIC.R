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
  (topmodels <- rbind(topmod_beartbd[1,], topmod_bobtbd[2,], topmod_coytbd[1,], topmod_liontbd[1,], topmod_wolftbd[1,])) #topmod_bobtbd[1,], 
  #'  Note: currently using 2nd most supported model for bobcat since w/in 0.65 deltaDIC of top model 
  
  
  #'  ------------------------
  #### PREY - FOCAL PREDATOR  ####
  #'  ------------------------
  #'  Model selection for models estimating time between detections of prey followed by predators
  
  #####  Black bear models  ####
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bear_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bear_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bear_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bear_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bear_global.RData")          
  #'  List for model selection
  bear_nt_tbd_list <- list(tbd.nt.bear.null, tbd.nt.bear.sppID, tbd.nt.bear.div, tbd.nt.bear.preyabund, tbd.nt.bear.global) 
  bear_nt_tbd_name <- c("tbd.nt.bear.null", "tbd.nt.bear.sppID", "tbd.nt.bear.div", "tbd.nt.bear.preyabund", "tbd.nt.bear.global") 
  
  #####  Bobcat models  ####
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bob_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bob_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bob_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bob_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.bob_global.RData")           
  #'  List for model selection
  bob_nt_tbd_list <- list(tbd.nt.bob.null, tbd.nt.bob.sppID, tbd.nt.bob.div, tbd.nt.bob.preyabund, tbd.nt.bob.global)  
  bob_nt_tbd_name <- c("tbd.nt.bob.null", "tbd.nt.bob.sppID", "tbd.nt.bob.div", "tbd.nt.bob.preyabund", "tbd.nt.bob.global")  
  
  #####  Coyote models  ####
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.coy_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.coy_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.coy_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.coy_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.coy_global.RData")
  #'  List for model selection
  coy_nt_tbd_list <- list(tbd.nt.coy.null, tbd.nt.coy.sppID, tbd.nt.coy.div, tbd.nt.coy.preyabund, tbd.nt.coy.global) 
  coy_nt_tbd_name <- c("tbd.nt.coy.null", "tbd.nt.coy.sppID", "tbd.nt.coy.div", "tbd.nt.coy.preyabund", "tbd.nt.coy.global") 
  
  #####  Mountain lion models  ####
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.lion_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.lion_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.lion_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.lion_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.lion_global.RData")          
  #'  List for model selection
  lion_nt_tbd_list <- list(tbd.nt.lion.null, tbd.nt.lion.sppID, tbd.nt.lion.div, tbd.nt.lion.preyabund, tbd.nt.lion.global) 
  lion_nt_tbd_name <- c("tbd.nt.lion.null", "tbd.nt.lion.sppID", "tbd.nt.lion.div", "tbd.nt.lion.preyabund", "tbd.nt.lion.global") 
  
  #####  Wolf models  ####
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_intercept_only.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_sppID.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_preydiversity.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_global.RData")          
  #'  List for model selection
  wolf_nt_tbd_list <- list(tbd.nt.wolf.null, tbd.nt.wolf.sppID, tbd.nt.wolf.div, tbd.nt.wolf.preyabund, tbd.nt.wolf.global) 
  wolf_nt_tbd_name <- c("tbd.nt.wolf.null", "tbd.nt.wolf.sppID", "tbd.nt.wolf.div", "tbd.nt.wolf.preyabund", "tbd.nt.wolf.global") 
  
  #####  Model selection using DIC  ####
  (topmod_nt_beartbd <- dictab(cand.set = bear_nt_tbd_list, modnames = bear_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_bobtbd <- dictab(cand.set = bob_nt_tbd_list, modnames = bob_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_coytbd <- dictab(cand.set = coy_nt_tbd_list, modnames = coy_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_liontbd <- dictab(cand.set = lion_nt_tbd_list, modnames = lion_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_wolftbd <- dictab(cand.set = wolf_nt_tbd_list, modnames = wolf_nt_tbd_name, sort = TRUE)) 
  
  #'  Best supported model per species-pair
  (topmodels_nt <- rbind(topmod_nt_beartbd[1,], topmod_nt_bobtbd[1,], topmod_nt_coytbd[1,], topmod_nt_liontbd[1,], topmod_nt_wolftbd[1,]))
  
  
  #'  -----------------
  ####  RESULT TABLES  ####
  #'  -----------------
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
           Model = gsub("sppIDxdiv", "Model 6", Model),
           Model = gsub("sppIDxpreyabund", "Model 8", Model),
           Model = gsub("sppID.div", "Model 5", Model),
           Model = gsub("sppID.preyabund", "Model 7", Model),
           Model = gsub("sppID", "Model 2", Model),
           Model = gsub("div", "Model 3", Model),
           Model = gsub("preyabund", "Model 4", Model),
           Model = gsub("global", "Model 9", Model),
           Model_name = gsub("null", "Null", Model_name),
           Model_name = gsub("sppIDxdiv", "Competitor ID * prey diversity", Model_name),
           Model_name = gsub("sppIDxpreyabund", "Competitor ID * prey abundance", Model_name),
           Model_name = gsub("sppID.div", "Competitor ID + prey diversity", Model_name),
           Model_name = gsub("sppID.preyabund", "Competitor ID + prey abundance", Model_name),
           Model_name = gsub("sppID", "Competitor ID", Model_name),
           Model_name = gsub("div", "Prey diversity", Model_name),
           Model_name = gsub("preyabund", "Prey abundance", Model_name),
           Model_name = gsub("global", "Global", Model_name)) %>%
    relocate(Species, .before = DIC) %>%
    relocate(Model, .after = Species) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_list_DIC) <- c("Species", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Full table of models ranked by DIC for all species-pairs
  model_nt_list_DIC <- rbind(topmod_nt_beartbd, topmod_nt_bobtbd, topmod_nt_coytbd, topmod_nt_liontbd, topmod_nt_wolftbd) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    #'  Split model name based on placement of multiple periods
    #'  https://stackoverflow.com/questions/26265400/use-regex-in-r-to-retrieve-string-before-second-occurence-of-a-period
    mutate(Species = sub( "(^[^.]+[.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species = str_replace(Species, "tbd.nt.", ""),
           Species = str_replace(Species, "bear", "Black bear"),
           Species = str_replace(Species, "lion", "Mountain lion"),
           Species = str_replace(Species, "coy", "Coyote"),
           Species = str_replace(Species, "bob", "Bobcat"),
           Species = str_replace(Species, "wolf", "Wolf"),
           Species = gsub("\\..*", "", Species),
           Model_name = gsub("^.*\\.", "", Modnames),
           Model = gsub("null", "Model 1", Model_name),           
           Model = gsub("sppID", "Model 2", Model),
           Model = gsub("div", "Model 3", Model),
           Model = gsub("preyabund", "Model 4", Model),
           Model = gsub("global", "Model 5", Model),
           Model_name = gsub("null", "Null", Model_name),
           Model_name = gsub("sppID.div", "Species ID + prey diversity", Model_name),
           Model_name = gsub("sppID.preyabund", "Species ID + prey abundance", Model_name),
           Model_name = gsub("sppID", "Species ID", Model_name),
           Model_name = gsub("div", "Prey diversity", Model_name),
           Model_name = gsub("preyabund", "Prey abundance", Model_name),
           Model_name = gsub("global", "Global", Model_name)) %>%
    relocate(Species, .before = DIC) %>%
    relocate(Model, .after = Species) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_nt_list_DIC) <- c("Species", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/Tables/DIC_TBD_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/Tables/DIC_TBD_model_selection_results.csv")
  save(topmodels, file = "./Outputs/Tables/DIC_TBD_top_models.RData")
  save(model_list_DIC, file = "./Outputs/Tables/DIC_TBD_model_selection_results.RData")
  
  write.csv(topmodels_nt, file = "./Outputs/Tables/DIC_TBD_top_nt_models.csv")
  write.csv(model_nt_list_DIC, file = "./Outputs/Tables/DIC_TBD_nt_model_selection_results.csv")
  save(topmodels_nt, file = "./Outputs/Tables/DIC_TBD_top_nt_models.RData")
  save(model_nt_list_DIC, file = "./Outputs/Tables/DIC_TBD_nt_model_selection_results.RData")
  
  
  #'  Next stop, create figures to visualize these results
  
  
  
