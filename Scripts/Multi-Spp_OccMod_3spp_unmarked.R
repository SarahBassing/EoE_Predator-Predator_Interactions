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
  #'  Winter primary period is considered ....
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
  source("./Scripts/Format_data_3spp_occmod_unmarked.R")
  
  
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
  #'  
  #'  Testing hypothesis that co-occurrence is non-independent and that cattle/
  #'  hunter activity impacts occurrence and/or co-occurrence patterns between
  #'  predators and prey.
  #'  
  #'  Include a consistent set of additional covariates to account for habitat
  #'  variation and other factors we know influence occurrence and detection.
  #'  Use AIC for model selection.
  #'  =============================
  
  ####  DETECTION SUBMODEL  ####
  #'  Detection formulas
  detFormulas_trail <- c("~CameraFacing", "~CameraFacing", "~CameraFacing")
  detFormulas_setup <- c("~Setup", "~Setup", "~Setup") 
  
  ####  OCCUPANCY SUBMODEL  ####
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null1 <- c("~1", "~1", "~1", "~0", "~0", "~0", "~0") # no interactions
  occFormulas_null2 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~0") # 2-spp interactions
  occFormulas_null3 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~1") # 3-spp interaction
  
  #'  Question 2: Does co-occurrence or conditional occupancy vary with covariates?
  occFormulas_gmu1 <- c("~GMU",  "~GMU", "~GMU",
                        "~1", "~1", "~1", "~1")
  occFormulas_gmu2 <- c("~GMU",  "~GMU", "~GMU",
                        "~GMU", "~GMU", "~GMU", "~1")
  occFormulas_gmu3 <- c("~GMU",  "~GMU", "~GMU",
                        "~GMU", "~GMU", "~GMU", "~GMU")
  occFormulas_hab1 <- c("~perc_forest", "~perc_forest", "~perc_forest",
                        "~1", "~1", "~1", "~1")
  occFormulas_hab2 <- c("~perc_forest", "~perc_forest", "~perc_forest",
                        "~perc_forest", "~perc_forest", "~perc_forest", "~1")
  occFormulas_hab3 <- c("~perc_forest", "~perc_forest", "~perc_forest",
                        "~perc_forest", "~perc_forest", "~perc_forest", "~perc_forest")

  
  ####  Apex Predator Smr21  ####
  (apex_trail <- occuMulti(detFormulas_trail, occFormulas_null1, apex_smr21_UMF, silent = TRUE))
  (apex_setup <- occuMulti(detFormulas_setup, occFormulas_null1, apex_smr21_UMF, silent = TRUE))
  #' List of fitted models
  apex_det_fld <- fitList(apex_trail, apex_setup)
  #' Model selection
  modSel(apex_det_fld)
  
  (apex_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, apex_smr21_UMF, silent = TRUE))
  (apex_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, apex_smr21_UMF, silent = TRUE))
  # (apex_null3 <- occuMulti(detFormulas_setup, occFormulas_null3, apex_smr21_UMF, silent = TRUE))
  (apex_gmu1 <- occuMulti(detFormulas_setup, occFormulas_gmu1, apex_smr21_UMF, silent = TRUE))
  # (apex_gmu2 <- occuMulti(detFormulas_setup, occFormulas_gmu2, apex_smr21_UMF, silent = TRUE))
  # (apex_gmu3 <- occuMulti(detFormulas_setup, occFormulas_gmu3, apex_smr21_UMF, silent = TRUE))
  # (apex_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, apex_smr21_UMF, silent = TRUE))
  (apex_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, apex_smr21_UMF, silent = TRUE))
  (apex_hab3 <- occuMulti(detFormulas_setup, occFormulas_hab3, apex_smr21_UMF, silent = TRUE))
  
  apex_occ_fld <- fitList(apex_null1, apex_null2, apex_gmu1, apex_hab2, apex_hab3)
  #' Model selection
  modSel(apex_occ_fld)
  summary(apex_hab2)
  
  
  ####  COY-LION-WOLF UMF  ####
  (clw_trail <- occuMulti(detFormulas_trail, occFormulas_null1, coy_lion_wolf_UMF, silent = TRUE))
  (clw_setup <- occuMulti(detFormulas_setup, occFormulas_null1, coy_lion_wolf_UMF, silent = TRUE))
  #' List of fitted models
  clw_det_fld <- fitList(clw_trail, clw_setup)
  #' Model selection
  modSel(clw_det_fld)
  
  (clw_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coy_lion_wolf_UMF, silent = TRUE))
  (clw_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_null3 <- occuMulti(detFormulas_trail, occFormulas_null3, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_gmu1 <- occuMulti(detFormulas_trail, occFormulas_gmu1, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_gmu2 <- occuMulti(detFormulas_trail, occFormulas_gmu2, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_gmu3 <- occuMulti(detFormulas_trail, occFormulas_gmu3, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, coy_lion_wolf_UMF, silent = TRUE))
  # (clw_hab3 <- occuMulti(detFormulas_trail, occFormulas_hab3, coy_lion_wolf_UMF, silent = TRUE))
  clw_occ_fld <- fitList(clw_null1, clw_null2)
  #' Model selection
  modSel(clw_occ_fld)
  summary(clw_null2)
  

  ####  BOB-LION-WOLF UMF  ####
  (blw_trail <- occuMulti(detFormulas_trail, occFormulas_null1, bob_lion_wolf_UMF, silent = TRUE))
  (blw_setup <- occuMulti(detFormulas_setup, occFormulas_null1, bob_lion_wolf_UMF, silent = TRUE))
  #' List of fitted models
  blw_det_fld <- fitList(blw_trail, blw_setup)
  #' Model selection
  modSel(blw_det_fld)
  
  (blw_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_lion_wolf_UMF, silent = TRUE))
  (blw_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_lion_wolf_UMF, silent = TRUE))
  (blw_null3 <- occuMulti(detFormulas_trail, occFormulas_null3, bob_lion_wolf_UMF, silent = TRUE))
  (blw_gmu1 <- occuMulti(detFormulas_trail, occFormulas_gmu1, bob_lion_wolf_UMF, silent = TRUE))
  # (blw_gmu2 <- occuMulti(detFormulas_trail, occFormulas_gmu2, bob_lion_wolf_UMF, silent = TRUE))
  # (blw_gmu3 <- occuMulti(detFormulas_trail, occFormulas_gmu3, bob_lion_wolf_UMF, silent = TRUE))
  (blw_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_lion_wolf_UMF, silent = TRUE))
  # (blw_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, bob_lion_wolf_UMF, silent = TRUE))
  # (blw_hab3 <- occuMulti(detFormulas_trail, occFormulas_hab3, bob_lion_wolf_UMF, silent = TRUE))
  blw_occ_fld <- fitList(blw_null1, blw_null2, blw_null3, blw_gmu1, blw_hab1)
  #' Model selection
  modSel(blw_occ_fld)
  summary(blw_hab1)
  
  
  
  
  
  
  
  
  
  
  