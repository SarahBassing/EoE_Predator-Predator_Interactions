  #'  --------------------------------------------
  #'  Multi-Species Occupancy Models in unmarked 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  --------------------------------------------
  #'  Script to source data formatting scripts & run multi-species, single-season
  #'  occupancy models for EoE predator detections. Focusing on black bear, bobcat,
  #'  coyote, mountain lion, and wolf interactions. Cameras ran summer 2020, winter
  #'  2020-2021, and winter 2021. Co-occurrence models test whether predator
  #'  co-occurrence is non-independent and whether their occurrence, co-occurrence,
  #'  and detection are influenced by other variables of interest.
  #'  
  #'  Summer primary period is considered July 1 - Spet. 15, equating to 11 1-wk 
  #'  sampling periods. 
  #'   
  #'  
  #'  Encounter histories are generated with the Detection_histories_for_occmod.R
  #'  DH and covariate data formated for unmarked with Format_data_3spp_occmod_unmarked.R
  #'  --------------------------------------------
  
  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(unmarked)
  library(MuMIn)
  library(condformat)
  library(tidyverse)
  
  #'  Source script that formats detection and covariate data for unmarked
  source("./Scripts/MultiSpp_OccMod/Format_data_3spp_occmod_unmarked.R")
  
  
  ####  Multi-Species Occupancy models  ####
  #'  ==================================
  #'  Multi-species occupancy model --> occuMulti (pg 83) in unmarked manual
  #'  Occupancy formulas should match number/order of columns in fDesign matrix
  #'  i.e., one formula for each species/interaction of interest
  #'  Detection formulas should match number/order of species in list of DH
  #'  Use ~1 to estimate intercept without covariates
  #'  Use ~0 to fix a natural parameter to 0
  #'  E.g., occFormulas <- c("~1", "~1", "~1", "~0", "~1", "~1", "~0") estimates 
  #'  intercept for 1st order natural parameters (3 spp), fixes first 2nd order  
  #'  natural parameter to 0 but estimates intercept for other 2nd order parameters, 
  #'  and fixes 3rd order natural parameter to 0 (only applies if 3+ spp)
  #'  Covariates: Can use different covariates on different natural parameters, 
  #'  E.g., covs on 1st order parameters to explain single-spp occurrence 
  #'  regardless of other spp, covs on 2nd order parameters to explain co-occ
  #'  Use AIC for model selection.
  #'  =============================
  
  ####  DETECTION SUBMODEL  ####
  #'  Detection formulas
  detFormulas_trail <- c("~CameraFacing", "~CameraFacing", "~CameraFacing")
  detFormulas_setup <- c("~Setup", "~Setup", "~Setup") 
  detFormulas_height <- c("~Height", "~Height", "~Height")
  detFormulas_wolfact <- c("~wolf_activity", "~wolf_activity", "~1") 
  
  ####  OCCUPANCY SUBMODEL  ####
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null1 <- c("~1", "~1", "~1", "~0", "~0", "~0", "~0") # no interactions
  occFormulas_null2 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~0") # 2-spp interactions
  
  #'  Question 2: Does co-occurrence or conditional occupancy vary with covariates?
  occFormulas_hab1 <- c("~Elev", "~Elev", "~Elev",
                        "~Elev", "~Elev", "~Elev", "~0")
  occFormulas_hab2 <- c("~Elev + PercForest", "~Elev + PercForest", "~Elev + PercForest",
                       "~Elev + PercForest", "~Elev + PercForest", "~Elev + PercForest", "~0")
  occFormulas_group1 <- c("~1", "~1", "~1",
                         "~MinGroupSize", "~MinGroupSize", "~MinGroupSize", "~0")
  occFormulas_group2 <- c("~1", "~1", "~1", "~1", "~1", "~MinGroupSize", "~0")
  occFormulas_prey <- c("~Nungulate", "~Nungulate", "~Nungulate",
                        "~Nungulate", "~Nungulate", "~Nungulate", "~0")
  occFormulas_diversity <- c("~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd",
                             "~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd", "~0")
  occFormulas_anthro <- c("~Dist2Suburbs + Nhuman + Nlivestock", "~Dist2Suburbs + Nhuman + Nlivestock", "~Dist2Suburbs + Nhuman + Nlivestock",
                          "~Dist2Suburbs + Nhuman + Nlivestock", "~Dist2Suburbs + Nhuman + Nlivestock", "~Dist2Suburbs + Nhuman + Nlivestock", "~0")
  

  ####  Detection Sub-Model Selection V1  ####
  #'  Function run and compare detection sub-models for each species pairing
  choose_det_submod_wolf <- function(det_submod, umf) {
    print(det_trail <- occuMulti(det_submod[[1]], occFormulas_null2, umf, silent = TRUE))
    print(det_setup <- occuMulti(det_submod[[2]], occFormulas_null2, umf, silent = TRUE))
    print(det_height <- occuMulti(det_submod[[3]], occFormulas_null2, umf, silent = TRUE))
    #' List of fitted models
    det_fld <- fitList(det_trail, det_setup, det_height) 
    #' Model selection
    print(modSel(det_fld))
  }
  det_submod <- list(detFormulas_trail, detFormulas_setup, detFormulas_height)
  
  ####  Wolf-Bear-Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_lion_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbrl_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_lion_20s_umf, silent = TRUE))
  (wbrl_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_lion_20s_umf, silent = TRUE))
  (wbrl_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bear_lion_20s_umf, silent = TRUE))
  # (wbrl_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bear_lion_20s_umf, silent = TRUE))
  # (wbrl_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, wolf_bear_lion_20s_umf, silent = TRUE))
  # (wbrl_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, wolf_bear_lion_20s_umf, silent = TRUE))
  (wbrl_20s_prey <- occuMulti(detFormulas_setup, occFormulas_prey, wolf_bear_lion_20s_umf, silent = TRUE))
  (wbrl_20s_diversity <- occuMulti(detFormulas_setup, occFormulas_diversity, wolf_bear_lion_20s_umf, silent = TRUE))
  (wbrl_20s_anthro <- occuMulti(detFormulas_setup, occFormulas_anthro, wolf_bear_lion_20s_umf, silent = TRUE))
  
  wbrl_20s_occ_fld <- fitList(wbrl_20s_null1, wbrl_20s_null2, wbrl_20s_hab1, wbrl_20s_prey, wbrl_20s_diversity, wbrl_20s_anthro)
  #' Model selection
  modSel(wbrl_20s_occ_fld)
  summary(wbrl_20s_hab1)
  
  ####  Wolf-Bear-Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_lion_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbrl_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_lion_21s_umf, silent = TRUE))
  (wbrl_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_lion_21s_umf, silent = TRUE))
  (wbrl_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bear_lion_21s_umf, silent = TRUE))
  # (wbrl_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bear_lion_21s_umf, silent = TRUE))
  # (wbrl_21s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, wolf_bear_lion_21s_umf, silent = TRUE))
  # (wbrl_21s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, wolf_bear_lion_21s_umf, silent = TRUE))
  (wbrl_21s_prey <- occuMulti(detFormulas_setup, occFormulas_prey, wolf_bear_lion_21s_umf, silent = TRUE))
  # (wbrl_21s_diversity <- occuMulti(detFormulas_setup, occFormulas_diversity, wolf_bear_lion_21s_umf, silent = TRUE))
  # (wbrl_21s_anthro <- occuMulti(detFormulas_setup, occFormulas_anthro, wolf_bear_lion_21s_umf, silent = TRUE))
  
  wbrl_21s_occ_fld <- fitList(wbrl_21s_null1, wbrl_21s_null2, wbrl_21s_hab1, wbrl_21s_prey)
  #' Model selection
  modSel(wbrl_21s_occ_fld)
  summary(wbrl_21s_hab1)
  
  ####  Wolf-Bobcat-Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_lion_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbl_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bob_lion_20s_umf, silent = TRUE))
  (wbl_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bob_lion_20s_umf, silent = TRUE))
  (wbl_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bob_lion_20s_umf, silent = TRUE))
  (wbl_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bob_lion_20s_umf, silent = TRUE))
  # (wbl_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, wolf_bob_lion_20s_umf, silent = TRUE))
  # (wbl_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, wolf_bob_lion_20s_umf, silent = TRUE))
  (wbl_20s_prey <- occuMulti(detFormulas_setup, occFormulas_prey, wolf_bob_lion_20s_umf, silent = TRUE))
  # (wbl_20s_diversity <- occuMulti(detFormulas_setup, occFormulas_diversity, wolf_bob_lion_20s_umf, silent = TRUE))
  (wbl_20s_anthro <- occuMulti(detFormulas_setup, occFormulas_anthro, wolf_bob_lion_20s_umf, silent = TRUE)) # suspiciously large coefficient estimates
  
  wbl_20s_occ_fld <- fitList(wbl_20s_null1, wbl_20s_null2, wbl_20s_hab1, wbl_20s_hab2, wbl_20s_prey) #, wbl_20s_anthro
  #' Model selection
  modSel(wbl_20s_occ_fld)
  summary(wbl_20s_hab1)
  
  ####  Wolf-Bobcat-Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_lion_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wbl_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_bob_lion_21s_umf, silent = TRUE))
  (wbl_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_bob_lion_21s_umf, silent = TRUE))
  # (wbl_21s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_bob_lion_21s_umf, silent = TRUE))
  # (wbl_21s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_bob_lion_21s_umf, silent = TRUE))
  # (wbl_21s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1, wolf_bob_lion_21s_umf, silent = TRUE))
  # (wbl_21s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2, wolf_bob_lion_21s_umf, silent = TRUE))
  (wbl_21s_prey <- occuMulti(detFormulas_trail, occFormulas_prey, wolf_bob_lion_21s_umf, silent = TRUE))
  (wbl_21s_diversity <- occuMulti(detFormulas_trail, occFormulas_diversity, wolf_bob_lion_21s_umf, silent = TRUE)) 
  # (wbl_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_bob_lion_21s_umf, silent = TRUE)) 
  
  wbl_21s_occ_fld <- fitList(wbl_21s_null1, wbl_21s_null2, wbl_21s_hab1, wbl_21s_prey, wbl_21s_diversity) #, wbl_21s_anthro
  #' Model selection
  modSel(wbl_21s_occ_fld)
  summary(wbl_21s_hab1)
  
  ####  Wolf-Coyote-Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_lion_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wcl_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_lion_20s_umf, silent = TRUE))
  (wcl_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_lion_20s_umf, silent = TRUE))
  (wcl_20s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_coy_lion_20s_umf, silent = TRUE))
  # (wcl_20s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_coy_lion_20s_umf, silent = TRUE))
  # (wcl_20s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1, wolf_coy_lion_20s_umf, silent = TRUE))
  # (wcl_20s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2, wolf_coy_lion_20s_umf, silent = TRUE))
  (wcl_20s_prey <- occuMulti(detFormulas_trail, occFormulas_prey, wolf_coy_lion_20s_umf, silent = TRUE))
  (wcl_20s_diversity <- occuMulti(detFormulas_trail, occFormulas_diversity, wolf_coy_lion_20s_umf, silent = TRUE))
  (wcl_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_coy_lion_20s_umf, silent = TRUE))
  
  wcl_20s_occ_fld <- fitList(wcl_20s_null1, wcl_20s_null2, wcl_20s_hab1, wcl_20s_prey, wcl_20s_diversity, wcl_20s_anthro)
  #' Model selection
  modSel(wcl_20s_occ_fld)
  summary(wcl_20s_prey)
  
  ####  Wolf-Coyote-Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_lion_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wcl_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_lion_21s_umf, silent = TRUE))
  (wcl_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_lion_21s_umf, silent = TRUE))
  (wcl_21s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_coy_lion_21s_umf, silent = TRUE))
  (wcl_21s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_coy_lion_21s_umf, silent = TRUE))
  # (wcl_21s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1, wolf_coy_lion_21s_umf, silent = TRUE))
  # (wcl_21s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2, wolf_coy_lion_21s_umf, silent = TRUE))
  (wcl_21s_prey <- occuMulti(detFormulas_trail, occFormulas_prey, wolf_coy_lion_21s_umf, silent = TRUE))
  (wcl_21s_diversity <- occuMulti(detFormulas_trail, occFormulas_diversity, wolf_coy_lion_21s_umf, silent = TRUE)) # suspiciously large coefficient estimates
  (wcl_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_coy_lion_21s_umf, silent = TRUE)) # suspiciously large coefficient estimates
  
  wcl_21s_occ_fld <- fitList(wcl_21s_null1, wcl_21s_null2, wcl_21s_hab1, wcl_21s_hab2, wcl_21s_prey) #, wcl_21s_diversity, wcl_21s_anthro
  #' Model selection
  modSel(wcl_21s_occ_fld)
  summary(wcl_21s_prey) #wcl_21s_diversity
  
  
  
  ####  Detection Sub-Model Selection V2  ####
  #'  Function run and compare detection sub-models for each species pairing
  choose_det_submod <- function(det_submod, umf) {
    print(det_trail <- occuMulti(det_submod[[1]], occFormulas_null2, umf, silent = TRUE))
    print(det_setup <- occuMulti(det_submod[[2]], occFormulas_null2, umf, silent = TRUE))
    print(det_height <- occuMulti(det_submod[[3]], occFormulas_null2, umf, silent = TRUE))
    print(det_wolfact <- occuMulti(det_submod[[4]], occFormulas_null2, umf, silent = TRUE))
    #' List of fitted models
    det_fld <- fitList(det_trail, det_setup, det_height, det_wolfact) 
    #' Model selection
    print(modSel(det_fld))
  }
  det_submod <- list(detFormulas_trail, detFormulas_setup, detFormulas_height, detFormulas_wolfact)
  
  
  ####  Lion - Bear Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bear_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  
  
  
  
  
  
  
  
  
  