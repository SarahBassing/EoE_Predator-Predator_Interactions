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
  #'  Summer primary period is considered July 1 - Sept. 15, equating to 11 1-wk 
  #'  sampling periods. 
  #'  Winter primary period is considered Dec 1 - Feb. 1, equating to 8 1-wk
  #'  sampling periods.
  #'   
  #'  
  #'  Encounter histories are generated with the Detection_histories_for_occmod.R
  #'  DH and covariate data formatted for unmarked with Format_data_2spp_occmod_unmarked.R
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
  #'  Use AIC for model selection.
  #'  =============================
  
  ####  DETECTION SUBMODEL  ####
  #'  Detection formulas
  detFormulas_trail <- c("~CameraFacing", "~CameraFacing")
  detFormulas_setup <- c("~Setup", "~Setup") 
  detFormulas_height <- c("~Height", "~Height")
  detFormulas_wolfact <- c("~wolf_activity", "~wolf_activity") 
  
  ####  OCCUPANCY SUBMODEL  ####
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null1 <- c("~1", "~1", "~0") # no interactions
  occFormulas_null2 <- c("~1", "~1", "~1") # 2-spp interactions
  
  #'  Question 2: Does co-occurrence or conditional occupancy vary with covariates?
  occFormulas_hab1 <- c("~Elev + PercForest", "~Elev + PercForest", "~1")
  occFormulas_hab2 <- c("~Elev + PercForest", "~Elev + PercForest", "~Elev + PercForest")
  
  occFormulas_prey1 <- c("~Nungulate", "~Nungulate", "~1")
  occFormulas_prey2 <- c("~Nungulate", "~Nungulate", "~Nungulate")
  occFormulas_prey3 <- c("~Elev + PercForest + Nungulate", "~Elev + PercForest + Nungulate", "~Elev + PercForest + Nungulate")
  
  occFormulas_diversity1 <- c("~Nelk + Nmoose + Nmd + Nwtd", "~Nelk + Nmoose + Nmd + Nwtd", "~1")
  occFormulas_diversity2 <- c("~Nelk + Nmoose + Nmd + Nwtd", "~Nelk + Nmoose + Nmd + Nwtd", "~Nelk + Nmoose + Nmd + Nwtd")
  occFormulas_diversity3 <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd")
  
  #'  Log distance to nearest road b/c expect behavioral response to dampen with growing distance from road
  #'  Not logging distance to suburbia b/v expect a relatively linear decline in mortality risk further from suburban/urband areas
  # occFormulas_anthro1 <- c("~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~1")
  # occFormulas_anthro2 <- c("~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  # occFormulas_anthro3 <- c("~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  # 
  # occFormulas_anthro1b <- c("~Dist2Burbs + logNearestRd + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~1")
  # occFormulas_anthro2b <- c("~Dist2Burbs + logNearestRd + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  # occFormulas_anthro3b <- c("~Elev + PercForest + Dist2Burbs + logNearestRd + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  #'  Anthropogenic factors hypothesized to be associated with mortality risk
  occFormulas_anthro1 <- c("~Dist2Burbs + Nlivestock", "~Dist2Burbs + Nlivestock", "~1")
  occFormulas_anthro2 <- c("~Dist2Burbs + Nlivestock", "~Dist2Burbs + Nlivestock", "~Dist2Burbs + Nlivestock")
  occFormulas_anthro3 <- c("~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock")
  
  #'  Anthropogenic factors hypothesized to be associated with disturbance
  occFormulas_anthro1b <- c("~logNearestRd + Nhuman", "~logNearestRd + Nhuman", "~1")
  occFormulas_anthro2b <- c("~logNearestRd + Nhuman", "~logNearestRd + Nhuman", "~logNearestRd + Nhuman")
  occFormulas_anthro3b <- c("~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman")
  
  occFormula_global1a <- c("~Elev + PercForest + Nungulate + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nungulate + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~1")
  occFormula_global2a <- c("~Elev + PercForest + Nungulate + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nungulate + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nungulate + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  occFormula_global1b <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~1")
  occFormula_global2b <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  
  
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
  det_submod <- list(detFormulas_trail, detFormulas_setup, detFormulas_height, detFormulas_wolfact)
  
  
  ####  Wolf - Bear Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey1 <- occuMulti(detFormulas_setup, occFormulas_prey1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey2 <- occuMulti(detFormulas_setup, occFormulas_prey2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey3 <- occuMulti(detFormulas_setup, occFormulas_prey3, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity1 <- occuMulti(detFormulas_setup, occFormulas_diversity1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity2 <- occuMulti(detFormulas_setup, occFormulas_diversity2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity3 <- occuMulti(detFormulas_setup, occFormulas_diversity3, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global1a <- occuMulti(detFormulas_setup, occFormula_global1a, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global2a <- occuMulti(detFormulas_setup, occFormula_global2a, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global1b <- occuMulti(detFormulas_setup, occFormula_global1b, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global2b <- occuMulti(detFormulas_setup, occFormula_global2b, wolf_bear_20s_umf, silent = TRUE))
  wbr_20s_occ_fld <- fitList(wbr_20s_null1, wbr_20s_null2, wbr_20s_hab1, wbr_20s_hab2, 
                             wbr_20s_prey1, wbr_20s_prey2, wbr_20s_prey3, wbr_20s_diversity1, 
                             wbr_20s_diversity2, wbr_20s_diversity3, wbr_20s_anthro1, 
                             wbr_20s_anthro2, wbr_20s_anthro3, wbr_20s_global1a, 
                             wbr_20s_global2a, wbr_20s_global1b, wbr_20s_global2b) 
  #' Model selection
  modSel(wbr_20s_occ_fld)
  summary(wbr_20s_global1a)
  summary(wbr_20s_hab1)
  
  ####  Wolf - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_prey1 <- occuMulti(detFormulas_setup, occFormulas_prey1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_prey2 <- occuMulti(detFormulas_setup, occFormulas_prey2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_prey3 <- occuMulti(detFormulas_setup, occFormulas_prey3, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_diversity1 <- occuMulti(detFormulas_setup, occFormulas_diversity1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_diversity2 <- occuMulti(detFormulas_setup, occFormulas_diversity2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_diversity3 <- occuMulti(detFormulas_setup, occFormulas_diversity3, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_global1a <- occuMulti(detFormulas_setup, occFormula_global1a, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_global2a <- occuMulti(detFormulas_setup, occFormula_global2a, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_global1b <- occuMulti(detFormulas_setup, occFormula_global1b, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_global2b <- occuMulti(detFormulas_setup, occFormula_global2b, wolf_bear_21s_umf, silent = TRUE))  
  wbr_21s_occ_fld <- fitList(wbr_21s_null1, wbr_21s_null2, wbr_21s_hab1, wbr_21s_hab2, 
                             wbr_21s_prey1, wbr_21s_prey2, wbr_21s_prey3, wbr_21s_diversity1,
                             wbr_21s_diversity2, wbr_21s_diversity3, wbr_21s_anthro1,
                             wbr_21s_anthro3, wbr_21s_global1a, #wbr_21s_anthro2, wbr_21s_global2a, , wbr_21s_global2b
                             wbr_21s_global1b) 
  #' Model selection
  modSel(wbr_21s_occ_fld)
  summary(wbr_21s_global1a)
  summary(wbr_21s_hab1)
  
  
  ####  Wolf - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bob_20s_umf, silent = TRUE))
  # (wb_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, wolf_bob_20s_umf, silent = TRUE))
  # (wb_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, wolf_bob_20s_umf, silent = TRUE))
  # (wb_20s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, wolf_bob_20s_umf, silent = TRUE))
  # (wb_20s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, wolf_bob_20s_umf, silent = TRUE))
  wb_20s_occ_fld <- fitList(wb_20s_null1, wb_20s_null2, wb_20s_hab1, wb_20s_hab2) #, wb_20s_group1, wb_20s_group2, wb_20s_habgroup1, wb_20s_habgroup2
  #' Model selection
  modSel(wb_20s_occ_fld)
  summary(wb_20s_hab1)
  
  ####  Wolf - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_bob_21s_umf, silent = TRUE))
  # (wb_21s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1_nowolf, wolf_bob_21s_umf, silent = TRUE))
  # (wb_21s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2_nowolf, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_habgroup1 <- occuMulti(detFormulas_trail, occFormulas_habgroup1_nowolf, wolf_bob_21s_umf, silent = TRUE))
  # (wb_21s_habgroup2 <- occuMulti(detFormulas_trail, occFormulas_habgroup2_nowolf, wolf_bob_21s_umf, silent = TRUE))
  wb_21s_occ_fld <- fitList(wb_21s_null1, wb_21s_null2, wb_21s_hab1, wb_21s_hab2, wb_21s_habgroup1) #, wb_21s_group1, wb_21s_group2, wb_21s_habgroup1, wb_21s_habgroup2
  #' Model selection
  modSel(wb_21s_occ_fld)
  summary(wb_21s_hab1)
  
  
  ####  Wolf - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro1b <- occuMulti(detFormulas_setup, occFormulas_anthro1b, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro2b <- occuMulti(detFormulas_setup, occFormulas_anthro2b, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_anthro3b <- occuMulti(detFormulas_setup, occFormulas_anthro3b, wolf_coy_20s_umf, silent = TRUE))
  wc_20s_occ_fld <- fitList(wc_20s_null1, wc_20s_null2, wc_20s_hab1, wc_20s_hab2) #, wc_20s_group1, wc_20s_group2, wc_20s_habgroup1, wc_20s_habgroup2
  #' Model selection
  modSel(wc_20s_occ_fld)
  summary(wc_20s_hab1) # wc_20s_hab2 is actually top model but much larger standard errors - way less certain about estimates
  
  ####  Wolf - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_hab2 <- occuMulti(detFormulas_trail, occFormulas_hab2, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro1b <- occuMulti(detFormulas_setup, occFormulas_anthro1b, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro2b <- occuMulti(detFormulas_setup, occFormulas_anthro2b, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthro3b <- occuMulti(detFormulas_setup, occFormulas_anthro3b, wolf_coy_21s_umf, silent = TRUE))
  
  wc_21s_occ_fld <- fitList(wc_21s_null1, wc_21s_null2, wc_21s_hab1, wc_21s_hab2, wc_21s_group1) #, wc_21s_group1, wc_21s_group2, wc_21s_habgroup1, wc_21s_habgroup2
  #' Model selection
  modSel(wc_21s_occ_fld)
  summary(wc_21s_group1)
  summary(wc_21s_hab1)
  
  
  ####  Wolf - Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_prey1 <- occuMulti(detFormulas_setup, occFormulas_prey1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_prey2 <- occuMulti(detFormulas_setup, occFormulas_prey2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_prey3 <- occuMulti(detFormulas_setup, occFormulas_prey3, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_diversity1 <- occuMulti(detFormulas_setup, occFormulas_diversity1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_diversity2 <- occuMulti(detFormulas_setup, occFormulas_diversity2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_diversity3 <- occuMulti(detFormulas_setup, occFormulas_diversity3, wolf_lion_20s_umf, silent = TRUE)) # suspect coefficient magnitude
  (wl_20s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_lion_20s_umf, silent = TRUE)) # suspect coefficient magnitude
  (wl_20s_anthro1b <- occuMulti(detFormulas_setup, occFormulas_anthro1b, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_anthro2b <- occuMulti(detFormulas_setup, occFormulas_anthro2b, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_anthro3b <- occuMulti(detFormulas_setup, occFormulas_anthro3b, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_global1a <- occuMulti(detFormulas_setup, occFormula_global1a, wolf_lion_20s_umf, silent = TRUE))
  #(wl_20s_global2a <- occuMulti(detFormulas_setup, occFormula_global2a, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_global1b <- occuMulti(detFormulas_setup, occFormula_global1b, wolf_lion_20s_umf, silent = TRUE))
  #(wl_20s_global2b <- occuMulti(detFormulas_setup, occFormula_global2b, wolf_lion_20s_umf, silent = TRUE))
  wl_20s_occ_fld <- fitList(wl_20s_null1, wl_20s_null2, wl_20s_hab1, wl_20s_hab2, 
                            wl_20s_prey1, wl_20s_prey2, wl_20s_prey3, wl_20s_diversity1,
                            wl_20s_diversity2, wl_20s_diversity3, wl_20s_anthro1,
                            wl_20s_anthro2, wl_20s_anthro3, wl_20s_global1a,
                            wl_20s_global1b) # wl_20s_global2b, wl_20s_global2b
  #' Model selection
  modSel(wl_20s_occ_fld)
  summary(wl_20s_global1b)
  
  ####  Wolf - Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_prey1 <- occuMulti(detFormulas_setup, occFormulas_prey1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_prey2 <- occuMulti(detFormulas_setup, occFormulas_prey2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_prey3 <- occuMulti(detFormulas_setup, occFormulas_prey3, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_diversity1 <- occuMulti(detFormulas_setup, occFormulas_diversity1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_diversity2 <- occuMulti(detFormulas_setup, occFormulas_diversity2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_diversity3 <- occuMulti(detFormulas_setup, occFormulas_diversity3, wolf_lion_21s_umf, silent = TRUE)) # suspect coefficient magnitude
  (wl_21s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_lion_21s_umf, silent = TRUE)) # suspect coefficient magnitude
  (wl_21s_anthro1b <- occuMulti(detFormulas_setup, occFormulas_anthro1b, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_anthro2b <- occuMulti(detFormulas_setup, occFormulas_anthro2b, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_anthro3b <- occuMulti(detFormulas_setup, occFormulas_anthro3b, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_global1a <- occuMulti(detFormulas_setup, occFormula_global1a, wolf_lion_21s_umf, silent = TRUE))
  (wl_20s_global2a <- occuMulti(detFormulas_setup, occFormula_global2a, wolf_lion_21s_umf, silent = TRUE))
  (wl_20s_global1b <- occuMulti(detFormulas_setup, occFormula_global1b, wolf_lion_21s_umf, silent = TRUE))
  (wl_20s_global2b <- occuMulti(detFormulas_setup, occFormula_global2b, wolf_lion_21s_umf, silent = TRUE))  wl_21s_occ_fld <- fitList(wl_21s_null1, wl_21s_null2, wl_21s_hab1, wl_21s_hab2) #, wl_21s_group1, wl_21s_group2, wl_21s_habgroup1, wl_21s_habgroup2
  #' Model selection
  modSel(wl_21s_occ_fld)
  summary(wl_21s_hab1)
  
  
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
  (lbr_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2, lion_bear_20s_umf, silent = TRUE))
  lbr_20s_occ_fld <- fitList(lbr_20s_null1, lbr_20s_null2, lbr_20s_hab1, lbr_20s_hab2, lbr_20s_group1, lbr_20s_group2, lbr_20s_habgroup1, lbr_20s_habgroup2) 
  #' Model selection
  modSel(lbr_20s_occ_fld)
  summary(lbr_20s_habgroup1)
  
  ####  Lion - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = lion_bear_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (lbr_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, lion_bear_21s_umf, silent = TRUE))
  # (lbr_21s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, lion_bear_21s_umf, silent = TRUE))
  # (lbr_21s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, lion_bear_21s_umf, silent = TRUE))
  # (lbr_21s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, lion_bear_21s_umf, silent = TRUE))
  # (lbr_21s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, lion_bear_21s_umf, silent = TRUE))
  lbr_21s_occ_fld <- fitList(lbr_21s_null1, lbr_21s_null2, lbr_21s_hab1, lbr_21s_hab2) #, lbr_21s_group1, lbr_21s_group2, lbr_21s_habgroup1, lbr_21s_habgroup2
  #' Model selection
  modSel(lbr_21s_occ_fld)
  summary(lbr_21s_hab1)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
