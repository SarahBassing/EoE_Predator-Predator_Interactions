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
  #'  Define univariate sub-models
  #'  Detection sub-models
  detFormulas_trail <- c("~CameraFacing", "~CameraFacing")
  detFormulas_setup <- c("~Setup", "~Setup") 
  detFormulas_height <- c("~Height", "~Height")
  detFormulas_wolfact <- c("~wolf_activity", "~wolf_activity") 
  #'  Ooccupancy sub-models
  occFormulas_elev <- c("~Elev", "~Elev", "~Elev")
  occFormulas_elev2 <- c("~Elev + I(Elev^2)", "~Elev + I(Elev^2)", "~Elev + I(Elev^2)")
  occFormulas_pforest <- c("~PercForest", "~PercForest", "~PercForest")
  occFormulas_group <- c("~1", "~1", "~MinGroupSize")
  occFormulas_group_nowolf <- c("~1", "~MinGroupSize", "~MinGroupSize")
  occFormulas_ungulate <- c("~Nungulate", "~Nungulate", "~Nungulate")
  occFormulas_elk <- c("~Nelk", "~Nelk", "~Nelk")
  occFormulas_moose <- c("~Nmoose", "~Nmoose", "~Nmoose")
  occFormulas_md <- c("~Nmd", "~Nmd", "~Nmd")
  occFormulas_wtd <- c("~Nwtd", "~Nwtd", "~Nwtd")
  occFormulas_suburbs <- c("~Dist2Suburbs", "~Dist2Suburbs", "~Dist2Suburbs")
  occFormulas_rural <- c("~Dist2Rural", "~Dist2Rural", "~Dist2Rural")
  occFormulas_dist2rd <- c("~NearestRd", "~NearestRd", "~NearestRd")
  occFormulas_human <- c("~Nhuman", "~Nhuman", "~Nhuman")
  occFormulas_livestock <- c("~Nlivestock", "~Nlivestock", "~Nlivestock")
  
  #'  List sub-models
  det_submodels <- list(detFormulas_trail, detFormulas_setup, detFormulas_height, detFormulas_wolfact)
  occ_submodels <- list(occFormulas_elev, occFormulas_elev2, occFormulas_pforest, occFormulas_group,
                     occFormulas_group_nowolf, occFormulas_ungulate, occFormulas_elk, occFormulas_moose,
                     occFormulas_md, occFormulas_wtd, occFormulas_suburbs, occFormulas_rural,
                     occFormulas_dist2rd, occFormulas_human, occFormulas_livestock)
  
  #'  Function to run univariate models to get familiar with important variables
  univariate_mods <- function(umf, det_submod, occ_submod) {
    #'  Define null sub-models
    detFormulas_null <- c("~1", "~1")
    occFormulas_null <- c("~1", "~1", "~1")
    
    #'  Run univariate models while holding other paramters constant with null sub-model
    print(det_trail <- occuMulti(det_submod[[1]], occFormulas_null, umf, silent = TRUE))
    print(det_setup <- occuMulti(det_submod[[2]], occFormulas_null, umf, silent = TRUE))
    print(det_height <- occuMulti(det_submod[[3]], occFormulas_null, umf, silent = TRUE))
    #print(det_activity <- occuMulti(det_submod[[4]], occFormulas_null, umf, silent = TRUE))
    
    print(occ_elev <- occuMulti(detFormulas_null, occ_submod[[1]], umf, silent = TRUE))
    print(occ_elev2 <- occuMulti(detFormulas_null, occ_submod[[2]], umf, silent = TRUE))
    print(occ_pforest <- occuMulti(detFormulas_null, occ_submod[[3]], umf, silent = TRUE))
    print(occ_group <- occuMulti(detFormulas_null, occ_submod[[4]], umf, silent = TRUE))
    # print(occ_group_nowolf <- occuMulti(detFormulas_null, occ_submod[[5]], umf, silent = TRUE))
    print(occ_ungulate <- occuMulti(detFormulas_null, occ_submod[[6]], umf, silent = TRUE))
    print(occ_elk <- occuMulti(detFormulas_null, occ_submod[[7]], umf, silent = TRUE))
    print(occ_moose <- occuMulti(detFormulas_null, occ_submod[[8]], umf, silent = TRUE))
    print(occ_md <- occuMulti(detFormulas_null, occ_submod[[9]], umf, silent = TRUE))
    print(occ_wtd <- occuMulti(detFormulas_null, occ_submod[[10]], umf, silent = TRUE))
    print(occ_suburbs <- occuMulti(detFormulas_null, occ_submod[[11]], umf, silent = TRUE))
    print(occ_rural <- occuMulti(detFormulas_null, occ_submod[[12]], umf, silent = TRUE))
    print(occ_dist2rd <- occuMulti(detFormulas_null, occ_submod[[13]], umf, silent = TRUE))
    print(occ_human <- occuMulti(detFormulas_null, occ_submod[[14]], umf, silent = TRUE))
    print(occ_livestock <- occuMulti(detFormulas_null, occ_submod[[15]], umf, silent = TRUE))
  }
  ####  Wolf models  ####
  wolf_bear_20s_univar <- univariate_mods(wolf_bear_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, ungulate, moose, wtd, livestock, dist2suburbs
  wolf_bear_21s_univar <- univariate_mods(wolf_bear_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, ungulate, moose, wtd, livestock, dist2rural
  
  wolf_bob_20s_univar <- univariate_mods(wolf_bob_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, Elevation2, ungulate, elk, moose, wtd, dist2suburbs, human, livestock
  wolf_bob_20w_univar <- univariate_mods(wolf_bob_20w_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  none...
  wolf_bob_21s_univar <- univariate_mods(wolf_bob_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, ungulate, elk, moose, md, human, livestock
  
  wolf_coy_20s_univar <- univariate_mods(wolf_coy_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, ungulate, elk, moose, wtd, dist2suburbs, dist2rural, nearestrd, human, livestock
  wolf_coy_20w_univar <- univariate_mods(wolf_coy_20w_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  nearestrd
  wolf_coy_21s_univar <- univariate_mods(wolf_coy_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, ungulate, elk, moose, mule deer, wtd, Dist2suburbs, nearestrd, human, livestock
  
  wolf_lion_20s_univar <- univariate_mods(wolf_lion_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, ungulate, elk, moose, md, wtd, dist2rural, nearestrd (but likely convergence issues), human, livestock
  wolf_lion_20w_univar <- univariate_mods(wolf_lion_20w_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  none...
  wolf_lion_21s_univar <- univariate_mods(wolf_lion_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, Elevation2, forest, ungulate, elk, moose, human (but suspected convergence issues), livestock
  
  ####  Lion models  ####
  lion_bear_20s_univar <- univariate_mods(lion_bear_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, group size, Dist2Suburbs, Dist2Rural, livestock
  lion_bear_21s_univar <- univariate_mods(lion_bear_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, Elevation2 (sorta), forest, group size, wtd, Dist2rural, human
  
  lion_bob_20s_univar <- univariate_mods(lion_bob_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Group size, ungulate, elk, md, human (but looks a little funky)
  lion_bob_20w_univar <- univariate_mods(lion_bob_20w_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  none...
  lion_bob_21s_univar <- univariate_mods(lion_bob_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  group size, elk, Dist2Suburbs, human
  
  lion_coy_20s_univar <- univariate_mods(lion_coy_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  group size, elk, Dist2Suburbs, humans
  lion_coy_20w_univar <- univariate_mods(lion_coy_20w_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  nearestrd
  lion_coy_21s_univar <- univariate_mods(lion_coy_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, group size, ungulate, elk, Dist2Suburbs, nearestrd, livestock
  
  ####  Bear models  ####
  bear_bob_20s_univar <- univariate_mods(bear_bob_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, forest, group size, wtd, human
  bear_bob_21s_univar <- univariate_mods(bear_bob_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation, Elevation2 (ish), forest, group size, ungulate, elk, wtd, Dist2Rural 
  
  bear_coy_20s_univar <- univariate_mods(bear_coy_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  Elevation2, forest, group size, elk, moose, md, wtd, Dist2Suburbs, Dst2Rural, nearestrd, human, livestock
  bear_coy_21s_univar <- univariate_mods(bear_coy_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #' Elevation, Elevation2, forest, group size, ungulate, elk, wtd, Dist2Suburbs, human (although pretty large intercepts)
  
  ####  Meso models  ####
  bob_coy_20s_univar <- univariate_mods(bob_coy_20s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #'  forest, group size, elk, moose, md, wtd, Dist2Suburbs, dist2rural, nearestrd, human
  bob_coy_21s_univar <- univariate_mods(bob_coy_21s_umf, det_submod = det_submodels, occ_submod = occ_submodels)
  #' Elevation, Elevation2, forest, group size, ungulate, elk, wtd, dist2suburbs, dist2rural, nearestrd, human (although large intercepts)
  
  
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
  
  occFormulas_group1 <- c("~MinGroupSize", "~MinGroupSize", "~1")
  occFormulas_group2 <- c("~MinGroupSize", "~MinGroupSize", "~MinGroupSize")
  occFormulas_group <- c("~1", "~1", "~MinGroupSize")
  
  occFormulas_group1_nowolf <- c("~1", "~MinGroupSize", "~1")
  occFormulas_group2_nowolf <- c("~1", "~MinGroupSize", "~MinGroupSize")
  
  occFormulas_habgroup1 <- c("~Elev + PercForest + MinGroupSize", "~Elev + PercForest + MinGroupSize", "~1")
  occFormulas_habgroup2 <- c("~Elev + PercForest + MinGroupSize", "~Elev + PercForest + MinGroupSize", "~Elev + PercForest + MinGroupSize")
  
  occFormulas_habgroup1_nowolf <- c("~Elev + PercForest", "~Elev + PercForest + MinGroupSize", "~1")
  occFormulas_habgroup2_nowolf <- c("~Elev + PercForest", "~Elev + PercForest + MinGroupSize", "~Elev + PercForest + MinGroupSize")
  
  occFormulas_prey1 <- c("~Nungulate", "~Nungulate", "~1")
  occFormulas_prey2 <- c("~Nungulate", "~Nungulate", "~Nungulate")
  occFormulas_prey3 <- c("~Elev + PercForest + Nungulate", "~Elev + PercForest + Nungulate", "~Elev + PercForest + Nungulate")
  
  occFormulas_diversity1 <- c("~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd", "~1")
  occFormulas_diversity2 <- c("~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd", "~Nelk + Nmd + Nwtd")
  occFormulas_diversity3 <- c("~Elev + PercForest + Nelk + Nmd + Nwtd", "~Elev + PercForest + Nelk + Nmd + Nwtd", "~Elev + PercForest + Nelk + Nmd + Nwtd")
  
  occFormulas_anthro1 <- c("~Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~1")
  occFormulas_anthro2 <- c("~Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~Dist2Suburbs + NearestRd + Nhuman + Nlivestock")
  occFormulas_anthro3 <- c("~Elev + PercForest + Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Suburbs + NearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Suburbs + NearestRd + Nhuman + Nlivestock")
  
  # occFormula_topbottom1
  # occFormula_topbottom2
  
  
  
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
  # (wbr_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, wolf_bear_20s_umf, silent = TRUE))
  # (wbr_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, wolf_bear_20s_umf, silent = TRUE))
  # (wbr_20s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, wolf_bear_20s_umf, silent = TRUE))
  # (wbr_20s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey1 <- occuMulti(detFormulas_setup, occFormulas_prey1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey2 <- occuMulti(detFormulas_setup, occFormulas_prey2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_prey3 <- occuMulti(detFormulas_setup, occFormulas_prey3, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity1 <- occuMulti(detFormulas_setup, occFormulas_diversity1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity2 <- occuMulti(detFormulas_setup, occFormulas_diversity2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_diversity3 <- occuMulti(detFormulas_setup, occFormulas_diversity3, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro1 <- occuMulti(detFormulas_setup, occFormulas_anthro1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro2 <- occuMulti(detFormulas_setup, occFormulas_anthro2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro3 <- occuMulti(detFormulas_setup, occFormulas_anthro3, wolf_bear_20s_umf, silent = TRUE))
  wbr_20s_occ_fld <- fitList(wbr_20s_null1, wbr_20s_null2, wbr_20s_hab1, wbr_20s_hab2) #, wbr_20s_group1, wbr_20s_group2, wbr_20s_habgroup1, wbr_20s_habgroup2
  #' Model selection
  modSel(wbr_20s_occ_fld)
  summary(wbr_20s_hab1)
  
  ####  Wolf - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, wolf_bear_21s_umf, silent = TRUE))
  # (wbr_21s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, wolf_bear_21s_umf, silent = TRUE))
  wbr_21s_occ_fld <- fitList(wbr_21s_null1, wbr_21s_null2, wbr_21s_hab1, wbr_21s_hab2) #, wbr_21s_group1, wbr_21s_group2, wbr_21s_habgroup1, wbr_21s_habgroup2
  #' Model selection
  modSel(wbr_21s_occ_fld)
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
  # (wc_20s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1_nowolf, wolf_coy_20s_umf, silent = TRUE))
  # (wc_20s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2_nowolf, wolf_coy_20s_umf, silent = TRUE))
  # (wc_20s_habgroup1 <- occuMulti(detFormulas_trail, occFormulas_habgroup1_nowolf, wolf_coy_20s_umf, silent = TRUE))
  # (wc_20s_habgroup2 <- occuMulti(detFormulas_trail, occFormulas_habgroup2_nowolf, wolf_coy_20s_umf, silent = TRUE))
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
  (wc_21s_group1 <- occuMulti(detFormulas_trail, occFormulas_group1_nowolf, wolf_coy_21s_umf, silent = TRUE))
  # (wc_21s_group2 <- occuMulti(detFormulas_trail, occFormulas_group2_nowolf, wolf_coy_21s_umf, silent = TRUE))
  # (wc_21s_habgroup1 <- occuMulti(detFormulas_trail, occFormulas_habgroup1_nowolf, wolf_coy_21s_umf, silent = TRUE))
  # (wc_21s_habgroup2 <- occuMulti(detFormulas_trail, occFormulas_habgroup2_nowolf, wolf_coy_21s_umf, silent = TRUE))
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
  # (wl_20s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, wolf_lion_20s_umf, silent = TRUE))
  # (wl_20s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, wolf_lion_20s_umf, silent = TRUE))
  # (wl_20s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, wolf_lion_20s_umf, silent = TRUE))
  # (wl_20s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, wolf_lion_20s_umf, silent = TRUE))
  wl_20s_occ_fld <- fitList(wl_20s_null1, wl_20s_null2, wl_20s_hab1, wl_20s_hab2) #, wl_20s_group1, wl_20s_group2, wl_20s_habgroup1, wl_20s_habgroup2
  #' Model selection
  modSel(wl_20s_occ_fld)
  summary(wl_20s_hab1)
  
  ####  Wolf - Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab1 <- occuMulti(detFormulas_setup, occFormulas_hab1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab2 <- occuMulti(detFormulas_setup, occFormulas_hab2, wolf_lion_21s_umf, silent = TRUE))
  # (wl_21s_group1 <- occuMulti(detFormulas_setup, occFormulas_group1_nowolf, wolf_lion_21s_umf, silent = TRUE))
  # (wl_21s_group2 <- occuMulti(detFormulas_setup, occFormulas_group2_nowolf, wolf_lion_21s_umf, silent = TRUE))
  # (wl_21s_habgroup1 <- occuMulti(detFormulas_setup, occFormulas_habgroup1_nowolf, wolf_lion_21s_umf, silent = TRUE))
  # (wl_21s_habgroup2 <- occuMulti(detFormulas_setup, occFormulas_habgroup2_nowolf, wolf_lion_21s_umf, silent = TRUE))
  wl_21s_occ_fld <- fitList(wl_21s_null1, wl_21s_null2, wl_21s_hab1, wl_21s_hab2) #, wl_21s_group1, wl_21s_group2, wl_21s_habgroup1, wl_21s_habgroup2
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
