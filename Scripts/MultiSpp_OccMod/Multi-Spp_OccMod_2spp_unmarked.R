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
  source("./Scripts/MultiSpp_OccMod/Format_data_2spp_occmod_unmarked.R")
  
  
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
  detFormulas_trail <- c("~CameraFacing", "~CameraFacing")
  detFormulas_setup <- c("~Setup", "~Setup") 
  detFormulas_wolfact <- c("~wolf_activity", "~wolf_activity") 
  
  ####  OCCUPANCY SUBMODEL  ####
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null1 <- c("~1", "~1", "~0") # no interactions
  occFormulas_null2 <- c("~1", "~1", "~1") # 2-spp interactions
  
  #'  Question 2: Does co-occurrence or conditional occupancy vary with covariates?
  occFormulas_gmu1 <- c("~GMU",  "~GMU", "~1")
  occFormulas_gmu2 <- c("~GMU",  "~GMU", "~GMU")

  occFormulas_hab1 <- c("~perc_forest", "~perc_forest", "~1")
  occFormulas_hab2 <- c("~perc_forest", "~perc_forest", "~perc_forest")
  
  occFormulas_group1 <- c("~min_group_size", "~min_group_size", "~1")
  occFormulas_group2 <- c("~min_group_size", "~min_group_size", "~min_group_size")
  
  occFormulas_habgroup <- c("~perc_forest", "~perc_forest", "~min_group_size")
  
  occFormula_bblmort1 <- c("~Bear_mort_km2", "~Lion_mort_n", "~1")
  occFormula_bblmort2 <- c("~Bear_mort_km2", "~Lion_mort_n", "~Lion_mort_n")
  occFormula_blmort1 <- c("~Bob_mort_km2", "~Lion_mort_n", "~1")
  occFormula_blmort2 <- c("~Bob_mort_km2", "~Lion_mort_n", "~Lion_mort_n")
  occFormula_clmort1 <- c("~1", "~Lion_mort_n", "~1")
  occFormula_clmort2 <- c("~1", "~Lion_mort_n", "~Lion_mort_n")
  occFormula_wolfmort <- c("~1", "~1", "~Wolf_mort_km2")
  
  
  ####  Bear-Lion Smr21  ####
  (bbl_trail <- occuMulti(detFormulas_trail, occFormulas_null1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_setup <- occuMulti(detFormulas_setup, occFormulas_null1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_wolfact <- occuMulti(detFormulas_wolfact, occFormulas_null1, bear_lion_smr20_UMF, silent = TRUE))
  #' List of fitted models
  bbl_det_fld <- fitList(bbl_trail, bbl_setup, bbl_wolfact)
  #' Model selection
  modSel(bbl_det_fld)
  
  (bbl_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_gmu1 <- occuMulti(detFormulas_setup, occFormulas_gmu1, bear_lion_smr20_UMF, silent = TRUE))
  # (bbl_gmu2 <- occuMulti(detFormulas_setup, occFormulas_gmu2, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_habgroup <- occuMulti(detFormulas_setup, occFormulas_habgroup, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_mort1 <- occuMulti(detFormulas_setup, occFormula_bblmort1, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_mort2 <- occuMulti(detFormulas_setup, occFormula_bblmort2, bear_lion_smr20_UMF, silent = TRUE))
  (bbl_mort3 <- occuMulti(detFormulas_setup, occFormula_wolfmort, bear_lion_smr20_UMF, silent = TRUE))
  bbl_occ_fld <- fitList(bbl_null1, bbl_null2, bbl_gmu1, bbl_hab1, bbl_hab2, bbl_group1, bbl_group2, bbl_habgroup, bbl_mort1, bbl_mort2, bbl_mort3)
  #' Model selection
  modSel(bbl_occ_fld)
  summary(bbl_habgroup)
  
  
  ####  COY-LION  ####
  (cl_trail <- occuMulti(detFormulas_trail, occFormulas_null1, coy_lion_smr20_UMF, silent = TRUE))
  (cl_setup <- occuMulti(detFormulas_setup, occFormulas_null1, coy_lion_smr20_UMF, silent = TRUE))
  (cl_wolfact <- occuMulti(detFormulas_wolfact, occFormulas_null1, coy_lion_smr20_UMF, silent = TRUE))
  #' List of fitted models
  cl_det_fld <- fitList(cl_trail, cl_setup, cl_wolfact)
  #' Model selection
  modSel(cl_det_fld)
  
  (cl_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coy_lion_smr20_UMF, silent = TRUE))
  (cl_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, coy_lion_smr20_UMF, silent = TRUE))
  (cl_gmu1 <- occuMulti(detFormulas_trail, occFormulas_gmu1, coy_lion_smr20_UMF, silent = TRUE))
  (cl_gmu2 <- occuMulti(detFormulas_trail, occFormulas_gmu2, coy_lion_smr20_UMF, silent = TRUE))
  (cl_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_lion_smr20_UMF, silent = TRUE))
  (cl_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, coy_lion_smr20_UMF, silent = TRUE))
  # (cl_group1 <- occuMulti(detFormulas_trail, occFormulas_group1, coy_lion_smr20_UMF, silent = TRUE))
  # (cl_group2 <- occuMulti(detFormulas_trail, occFormulas_group2, coy_lion_smr20_UMF, silent = TRUE))
  (cl_habgroup <- occuMulti(detFormulas_setup, occFormulas_habgroup, coy_lion_smr20_UMF, silent = TRUE))
  (cl_mort1 <- occuMulti(detFormulas_setup, occFormula_clmort1, coy_lion_smr20_UMF, silent = TRUE))
  # (cl_mort2 <- occuMulti(detFormulas_setup, occFormula_clmort2, coy_lion_smr20_UMF, silent = TRUE))
  (cl_mort3 <- occuMulti(detFormulas_setup, occFormula_wolfmort, coy_lion_smr20_UMF, silent = TRUE))
  cl_occ_fld <- fitList(cl_null1, cl_null2, cl_gmu1, cl_gmu2, cl_hab1, cl_hab2, cl_habgroup, cl_mort1, cl_mort3)
  #' Model selection
  modSel(cl_occ_fld)
  summary(cl_habgroup)
  
  
  ####  BOB-LION  ####
  (bl_trail <- occuMulti(detFormulas_trail, occFormulas_null1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_setup <- occuMulti(detFormulas_setup, occFormulas_null1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_wolfact <- occuMulti(detFormulas_wolfact, occFormulas_null1, bob_lion_smr20_UMF, silent = TRUE))
  #' List of fitted models
  bl_det_fld <- fitList(bl_trail, bl_setup, bl_wolfact)
  #' Model selection
  modSel(bl_det_fld)
  
  (bl_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, bob_lion_smr20_UMF, silent = TRUE))
  (bl_gmu1 <- occuMulti(detFormulas_setup, occFormulas_gmu1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_gmu2 <- occuMulti(detFormulas_setup, occFormulas_gmu2, bob_lion_smr20_UMF, silent = TRUE))
  (bl_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, bob_lion_smr20_UMF, silent = TRUE))
  (bl_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, bob_lion_smr20_UMF, silent = TRUE))
  (bl_habgroup <- occuMulti(detFormulas_setup, occFormulas_habgroup, bob_lion_smr20_UMF, silent = TRUE))
  (bl_mort1 <- occuMulti(detFormulas_setup, occFormula_blmort1, bob_lion_smr20_UMF, silent = TRUE))
  (bl_mort2 <- occuMulti(detFormulas_setup, occFormula_blmort2, bob_lion_smr20_UMF, silent = TRUE))
  (bl_mort3 <- occuMulti(detFormulas_setup, occFormula_wolfmort, bob_lion_smr20_UMF, silent = TRUE))
  bl_occ_fld <- fitList(bl_null1, bl_null2, bl_gmu1, bl_gmu2, bl_hab1, bl_hab2, bl_group1, bl_group2, bl_habgroup, bl_mort1, bl_mort2, bl_mort3)
  #' Model selection
  modSel(bl_occ_fld)
  summary(bl_group2)
  
  
  
  
  
  
  
  
  
  
