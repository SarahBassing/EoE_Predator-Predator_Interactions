  #'  ---------------------------
  #'  Univariate Multi-Spp Occupancy Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2023
  #'  --------------------------------------------
  #'  Script to source data formatting scripts & run univariate-ish multi-species, 
  #'  single-season occupancy models for EoE predator detections. Goal is to 
  #'  review general relationships (+/-) of each variable with respect to detection
  #'  or occurrence/co-occurrence parameters to understand general patterns. Not
  #'  using this analysis to guide model building, only using this to compare
  #'  final models to make sure relationships didn't change unexpectedly, indicating
  #'  potential model issues (e.g., collinearity).
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
  occFormulas_suburbs <- c("~Dist2Burbs", "~Dist2Burbs", "~Dist2Burbs")
  occFormulas_rural <- c("~Dist2Rrl", "~Dist2Rrl", "~Dist2Rrl")
  occFormulas_dist2rd <- c("~logNearestRd", "~logNearestRd", "~logNearestRd")
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
  
    
    
    
    