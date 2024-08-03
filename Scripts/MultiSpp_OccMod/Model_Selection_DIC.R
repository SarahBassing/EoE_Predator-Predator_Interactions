  #'  ---------------------------------
  #'  Model selection & result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ---------------------------------
  #'  Script to identify most supported model using DIC and format result tables.
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(AICcmodavg)
  library(tidyverse)
  
  #'  Load model outputs and list
  #'  Wolf-Bear                       
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_2024-07-16.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2024-07-17.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-19.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_newDIV_2024-07-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-19.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-20.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-20.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(global)_psix(global)_p(setup_effort)_2024-07-20.RData")  
  wolfbear_list <- list(wolf.bear.null, wolf.bear.hab, wolf.bear.preyabund, wolf.bear.preydiv, wolf.bear.habx, wolf.bear.preyabundx, wolf.bear.preydivx, wolf.bear.global) 
  wolfbear_name <- c("wolf.bear.null", "wolf.bear.hab", "wolf.bear.preyabund", "wolf.bear.preydiv", "wolf.bear.habx", "wolf.bear.preyabundx", "wolf.bear.preydivx", "wolf.bear.global") 
  
  #'  Wolf-Coyote                     
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(yr)_p(.)_2024-07-20.RData")  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2024-07-21.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-21.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-21.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-22.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-22.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-23.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(global)_psix(global)_p(setup_effort)_2024-07-23.RData")  
  wolfcoy_list <- list(wolf.coy.null, wolf.coy.hab, wolf.coy.preyabund, wolf.coy.preydiv, wolf.coy.habx, wolf.coy.preyabundx, wolf.coy.preydivx, wolf.coy.global) 
  wolfcoy_name <- c("wolf.coy.null", "wolf.coy.hab", "wolf.coy.preyabund", "wolf.coy.preydiv", "wolf.coy.habx", "wolf.coy.preyabundx", "wolf.coy.preydivx", "wolf.coy.global")
  
  #'  Wolf-Lion                           
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-07-23.RData")  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_p(setup_effort)_2024-07-24.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-24.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-25.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-25.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-26.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-26.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(global)_psix(global)_p(setup_effort)_2024-07-27.RData") 
  wolflion_list <- list(wolf.lion.null, wolf.lion.hab, wolf.lion.preyabund, wolf.lion.preydiv, wolf.lion.habx, wolf.lion.preyabundx, wolf.lion.preydivx, wolf.lion.global) 
  wolflion_name <- c("wolf.lion.null", "wolf.lion.hab", "wolf.lion.preyabund", "wolf.lion.preydiv", "wolf.lion.habx", "wolf.lion.preyabundx", "wolf.lion.preydivx", "wolf.lion.global") 
  
  #'  Lion-Bear                           
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2024-07-28.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_p(setup_effort)_2024-07-28.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-28.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-29.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-29.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(global)_psix(global)_p(setup_effort)_2024-07-31.RData") 
  lionbear_list <- list(lion.bear.null, lion.bear.hab, lion.bear.preyabund, lion.bear.preydiv, lion.bear.habx, lion.bear.preyabundx, lion.bear.preydivx, lion.bear.global)
  lionbear_name <- c("lion.bear.null", "lion.bear.hab", "lion.bear.preyabund", "lion.bear.preydiv", "lion.bear.habx", "lion.bear.preyabundx", "lion.bear.preydivx", "lion.bear.global")
  
  #'  Lion-Bobcat                         
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2024-07-26.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_p(setup_effort)_2024-07-27.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-27.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-28.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-28.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-28.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-29.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(global)_psix(global)_p(setup_effort)_2024-07-29.RData") 
  lionbob_list <- list(lion.bob.null, lion.bob.hab, lion.bob.preyabund, lion.bob.preydiv, lion.bob.habx, lion.bob.preyabundx, lion.bob.preydivx, lion.bob.global)
  lionbob_name <- c("lion.bob.null", "lion.bob.hab", "lion.bob.preyabund", "lion.bob.preydiv", "lion.bob.habx", "lion.bob.preyabundx", "lion.bob.preydivx", "lion.bob.global")
  
  #'  Coyote-Bobcat                # April 10-11 runs ni = 100000, nt = 10; April 4 runs ni = 50000, nt = 10 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(yr)_p(.)_2024-07-22.RData")  #2023-04-10
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_p(setup_effort)_2024-07-23.RData") #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-23.RData") #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-24.RData")  #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-25.RData") #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-26.RData") #2023-04-11
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_2024-07-26.RData") #2023-04-11
  coybob_list <- list(coy.bob.null, coy.bob.hab, coy.bob.preyabund, coy.bob.preydiv, coy.bob.habx, coy.bob.preyabundx, coy.bob.preydivx, coy.bob.global) 
  coybob_name <- c("coy.bob.null", "coy.bob.hab", "coy.bob.preyabund", "coy.bob.preydiv", "coy.bob.habx", "coy.bob.preyabundx", "coy.bob.preydivx", "coy.bob.global")
  
  #'  Wolf-Bear                       
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(yr)_p(.)_2024-07-22.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_p(setup_effort)_2024-07-23.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_preyabund_yr)_p(setup_effort)_2024-07-23.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_preydiversity_yr)_p(setup_effort)_2024-07-24.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2024-07-25.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2024-07-25.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(global)_psix(global)_p(setup_effort)_2024-07-26.RData")  
  bearcoy_list <- list(bear.coy.null, bear.coy.hab, bear.coy.preyabund, bear.coy.preydiv, bear.coy.habx, bear.coy.preyabundx, bear.coy.preydivx, bear.coy.global) 
  bearcoy_name <- c("bear.coy.null", "bear.coy.hab", "bear.coy.preyabund", "bear.coy.preydiv", "bear.coy.habx", "bear.coy.preyabundx", "bear.coy.preydivx", "bear.coy.global") 
  
  
  ####  Create model selection table using DIC, deltaDIC, and model weights  ####
  (topmod_wolfbear <- dictab(cand.set = wolfbear_list, modnames = wolfbear_name, sort = TRUE)) 
  (topmod_wolfcoy <- dictab(cand.set = wolfcoy_list, modnames = wolfcoy_name, sort = TRUE)) 
  (topmod_wolflion <- dictab(cand.set = wolflion_list, modnames = wolflion_name, sort = TRUE)) 
  (topmod_lionbear <- dictab(cand.set = lionbear_list, modnames = lionbear_name, sort = TRUE)) 
  (topmod_lionbob <- dictab(cand.set = lionbob_list, modnames = lionbob_name, sort = TRUE)) 
  (topmod_coybob <- dictab(cand.set = coybob_list, modnames = coybob_name, sort = TRUE)) 
  (topmod_bearcoy <- dictab(cand.set = bearcoy_list, modnames = bearcoy_name, sort = TRUE)) 
  
  
  #'  Best supported model per species-pair
  (topmodels <- rbind(topmod_wolfbear[1,], topmod_wolfcoy[1,], topmod_wolflion[1,], topmod_lionbear[1,], topmod_lionbob[1,], topmod_coybob[1,], topmod_bearcoy[1,]))
  
  #'  Full table of models ranked by DIC for all species-pairs
  model_list_DIC <- rbind(topmod_wolfbear, topmod_wolfcoy, topmod_coybob, topmod_wolflion, topmod_lionbear, topmod_lionbob) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    #'  Split model name based on placement of multiple periods
    #'  https://stackoverflow.com/questions/26265400/use-regex-in-r-to-retrieve-string-before-second-occurence-of-a-period
    mutate(Species_pair = sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species_pair = str_replace(Species_pair, "\\.", " - "), 
           Species_pair = str_replace(Species_pair, "bear", "Black bear"),
           Species_pair = str_replace(Species_pair, "lion", "Mountain lion"),
           Species_pair = str_replace(Species_pair, "coy", "Coyote"),
           Species_pair = str_replace(Species_pair, "bob", "Bobcat"),
           Species_pair = str_replace(Species_pair, "wolf", "Wolf"),
           Model = gsub(".*null", "Model 1", Modnames), 
           Model = gsub(".*habx", "Model 5", Model), 
           Model = gsub(".*preyabundx", "Model 6", Model),
           Model = gsub(".*preydivx", "Model 7", Model),
           Model = gsub(".*hab", "Model 2", Model), 
           Model = gsub(".*preyabund", "Model 3", Model), 
           Model = gsub(".*preydiv", "Model 4", Model), 
           Model = gsub(".*global", "Model 8", Model),
           Model_name = gsub(".*null", "Null", Modnames),
           Model_name = gsub(".*habx", "Habitat with interaction", Model_name),
           Model_name = gsub(".*preyabundx", "Prey abundance with interaction", Model_name),
           Model_name = gsub(".*preydivx", "Prey diversity with interaction", Model_name),
           Model_name = gsub(".*hab", "Habitat, no interaction", Model_name),
           Model_name = gsub(".*preyabund", "Prey abundance, no interaction", Model_name),
           Model_name = gsub(".*preydiv", "Prey diversity, no interaction", Model_name),
           Model_name = gsub(".*global", "Global", Model_name)) %>%
    relocate(Species_pair, .before = DIC) %>%
    relocate(Model, .after = Species_pair) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_list_DIC) <- c("Predator pair", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/Tables/DIC_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/Tables/DIC_model_selection_results.csv")
  save(topmodels, file = "./Outputs/Tables/DIC_top_models.RData")
  save(model_list_DIC, file = "./Outputs/Tables/DIC_model_selection_results.RData")
  
  
    
    
    
    
  