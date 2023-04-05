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
  #'  Wolf-Bear                       # NOTE: 30,000 iterations, thinning rate of 5
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(.)_p(.)_2023-03-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2023-04-03.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preyabund_yr)_p(setup_effort)_2023-03-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-03-30.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-03-31.RData") # consider rerunning
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2023-03-31.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2023-03-31.RData")
  wolfbear_list <- list(wolf.bear.null, wolf.bear.hab, wolf.bear.preyabund, wolf.bear.preydiv, wolf.bear.habx, wolf.bear.preyabundx, wolf.bear.preydivx) 
  wolfbear_name <- c("wolf.bear.null", "wolf.bear.hab", "wolf.bear.preyabund", "wolf.bear.preydiv", "wolf.bear.habx", "wolf.bear.preyabundx", "wolf.bear.preydivx") 
  
  #'  Wolf-Coyote                     # April runs with 50000 iterations, thinning rate of 10
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(.)_p(.)_2023-04-05.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preyabund_yr)_p(setup_effort)_2023-04-05.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-05.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2023-04-04.RData") 
  wolfcoy_list <- list(wolf.coy.null, wolf.coy.hab, wolf.coy.preyabund, wolf.coy.preydiv, wolf.coy.habx, wolf.coy.preyabundx, wolf.coy.preydivx) 
  wolfcoy_name <- c("wolf.coy.null", "wolf.coy.hab", "wolf.coy.preyabund", "wolf.coy.preydiv", "wolf.coy.habx", "wolf.coy.preyabundx", "wolf.coy.preydivx")
  
  #'  Wolf-Lion     #  April runs with 50000 iterations, thinning rate of 10
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(.)_p(.)_2023-04-05.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preyabund_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preyabund_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2023-04-05.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2023-04-05a.RData") # currently re-running though
  wolflion_list <- list(wolf.lion.null, wolf.lion.hab, wolf.lion.preyabund, wolf.lion.preydiv, wolf.lion.habx, wolf.lion.preyabundx, wolf.lion.preydivx) 
  wolflion_name <- c("wolf.lion.null", "wolf.lion.hab", "wolf.lion.preyabund", "wolf.lion.preydiv", "wolf.lion.habx", "wolf.lion.preyabundx", "wolf.lion.preydivx") 
  
  #'  Lion-Bear
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(.)_p(.)_2023-03-31.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preyabund_yr)_p(setup_effort)_2023-03-31.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-01.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-04-01.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort).RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).RData")
  lionbear_list <- list(lion.bear.null, lion.bear.hab, lion.bear.preyabund, lion.bear.preydiv, lion.bear.habx)#, lion.bear.preyabundx, lion.bear.preydivx) 
  lionbear_name <- c("lion.bear.null", "lion.bear.hab", "lion.bear.preyabund", "lion.bear.preydiv", "lion.bear.habx")#, "lion.bear.preyabundx", "lion.bear.preydivx") 
  
  #'  Lion-Bobcat     # still running interaction  models   --- consider rerunning and providing starting values for mean psi & p
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(.)_p(.)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preyabund_yr)_p(setup_effort)_2023-04-05.RData") # did not converge even with 50,000 iterations, thinning at 10
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-05.RData") # converged but not perfect
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(.)_p(setup_effort).RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort).RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).RData")
  lionbob_list <- list(lion.bob.null, lion.bob.hab)#, lion.bob.preyabund, lion.bob.preydiv, lion.bob.habx, lion.bob.preyabundx, lion.bob.preydivx) 
  lionbob_name <- c("lion.bob.null", "lion.bob.hab")#, "lion.bob.preyabund", "lion.bob.preydiv", "lion.bob.habx", "lion.bob.preyabundx", "lion.bob.preydivx") 
  
  #'  Coyote-Bobcat                # NOTE: 50000 iterations and thinning rate of 10 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(.)_p(.)_2023-04-01.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preyabund_yr)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-04-04.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_2023-04-04.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_2023-04-04.RData")
  coybob_list <- list(coy.bob.null, coy.bob.hab, coy.bob.preyabund, coy.bob.preydiv, coy.bob.habx, coy.bob.preyabundx, coy.bob.preydivx) 
  coybob_name <- c("coy.bob.null", "coy.bob.hab", "coy.bob.preyabund", "coy.bob.preydiv", "coy.bob.habx", "coy.bob.preyabundx", "coy.bob.preydivx") 
  
  ####  Create model selection table using DIC, deltaDIC, and model weights  ####
  (topmod_wolfbear <- dictab(cand.set = wolfbear_list, modnames = wolfbear_name, sort = TRUE)) # wolf.bear.habx best supported model even though NS interaction term
  (topmod_wolfcoy <- dictab(cand.set = wolfcoy_list, modnames = wolfcoy_name, sort = TRUE)) # wolf.coy.hab best supported model even though interaction is signif, followed by wolf.coy.habx - habx likely getting dinged for poorer model convergence
  (topmod_wolflion <- dictab(cand.set = wolflion_list, modnames = wolflion_name, sort = TRUE)) # wolf.lion.null best supported model, followed by wolf.lion.habx ... probably b/c no habitat variables matter for lions apparently 
  (topmod_lionbear <- dictab(cand.set = lionbear_list, modnames = lionbear_name, sort = TRUE)) # lion.bear.null best supported model but still running prey interaction mods
  (topmod_lionbob <- dictab(cand.set = lionbob_list, modnames = lionbob_name, sort = TRUE)) # lion.bob.null currently best supported
  (topmod_coybob <- dictab(cand.set = coybob_list, modnames = coybob_name, sort = TRUE)) # coy.bob.habx best supported model, followed by coy.bob.preyabund & coy.bob.preyabundx 
  
  
  #'  Best supported model per species-pair
  topmodels <- rbind(topmod_coybob[1,])
  topmodels <- rbind(topmod_wolfbear[1,], topmod_wolfcoy[1,], topmod_wolflion[1,], topmod_lionbear[1,], topmod_lionbob[1,], topmod_coybob[1,])
  
  #'  Full table of models ranked by DIC for all species-pairs
  model_list_DIC <- rbind(topmod_coybob)
  model_list_DIC <- rbind(topmod_wolfbear, topmod_wolfcoy, topmod_wolflion, topmod_lionbear, topmod_lionbob, topmod_coybob)
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/MultiSpp_OccMod_Outputs/DIC_model_selection_results.csv")
  save(topmodels, file = "./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.RData")
  save(model_list_DIC, file = "./Outputs/MultiSpp_OccMod_Outputs/DIC_model_selection_results.RData")
  
  
    
    
    
    
  