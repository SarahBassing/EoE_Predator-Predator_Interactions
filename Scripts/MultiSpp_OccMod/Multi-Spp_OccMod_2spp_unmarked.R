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
  #'  Null models
  occFormulas_null1 <- c("~1", "~1", "~0") # no interactions
  occFormulas_null2 <- c("~1", "~1", "~1") # 2-spp interactions
  
  #'  Habitat filtering with basic habitat variables
  occFormulas_hab <- c("~Elev + PercForest", "~Elev + PercForest", "~Elev + PercForest")
  
  #'  Relative abundance of prey in function groups (big deer = elk & moose; 
  #'  small deer = mule & white-tailed deer; lagomorph = hares & rabbits)
  occFormulas_preygroups <- c("~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph")
  
  #'  Relative abundance of prey species (captures ful diversity of ungulate prey)
  occFormulas_preydiversity <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph")
  
  #'  Anthropogenic factors hypothesized to be associated with mortality risk
  occFormulas_anthromort <- c("~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock")
  
  #'  Anthropogenic factors hypothesized to be associated with disturbance
  #'  Log distance to nearest road b/c expect behavioral response to dampen with growing distance from road
  #'  Not logging distance to suburbia b/v expect a relatively linear decline in mortality risk further from suburban/urban areas
  occFormulas_anthrodist <- c("~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman")
  
  #'  Global models (1 = prey functional groups; 2 = prey species diversity)
  occFormula_global1 <- c("~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  occFormula_global2 <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  
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
  
  
  ####  Wolf - Bear Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_bear_20s_umf, silent = TRUE))
  wbr_20s_occ_fld <- fitList(wbr_20s_null1, wbr_20s_null2, wbr_20s_hab, wbr_20s_preygroups, wbr_20s_preydiversity, 
                             wbr_20s_anthromort, wbr_20s_anthrodist, wbr_20s_global1, wbr_20s_global2) 
  #' Model selection
  modSel(wbr_20s_occ_fld)
  summary(wbr_20s_anthromort)

  ####  Wolf - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_bear_21s_umf, silent = TRUE)) #fails with lagomorphs
  (wbr_21s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_bear_21s_umf, silent = TRUE)) 
  (wbr_21s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_bear_21s_umf, silent = TRUE)) #fail with lagomorphs
  (wbr_21s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_bear_21s_umf, silent = TRUE)) #fail with lagomorphs
  wbr_21s_occ_fld <- fitList(wbr_21s_null1, wbr_21s_null2, wbr_21s_hab, wbr_21s_preydiversity, #wbr_21s_preygroups, 
                             wbr_21s_anthromort, wbr_21s_anthrodist) #, wbr_21s_global1, wbr_21s_global2
  #' Model selection
  modSel(wbr_21s_occ_fld)
  summary(wbr_21s_anthromort) 
  
  
  ####  Wolf - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_bob_20s_umf, silent = TRUE)) #questionable coeffs
  (wb_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_bob_20s_umf, silent = TRUE)) # fail
  wb_20s_occ_fld <- fitList(wb_20s_null1, wb_20s_null2, wb_20s_hab, wb_20s_preygroups, wb_20s_preydiversity,
                            wb_20s_anthromort, wb_20s_anthrodist, wb_20s_global1) #, wb_20s_global2
  #' Model selection
  modSel(wb_20s_occ_fld)
  summary(wb_20s_preydiversity)  # global1 but questionable coeffs
  summary(wb_20s_anthromort) # ranked very close second to prey diversity
  
  ####  Wolf - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_bob_21s_umf, silent = TRUE)) #questionable coeffs
  (wb_21s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_bob_21s_umf, silent = TRUE)) #fail
  (wb_21s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_bob_21s_umf, silent = TRUE)) #questionable coeffs
  wb_21s_occ_fld <- fitList(wb_21s_null1, wb_21s_null2, wb_21s_hab, wb_21s_preygroups, wb_21s_preydiversity, 
                            wb_21s_anthromort, wb_21s_anthrodist, wb_21s_global2) #wb_21s_global1 
  #' Model selection
  modSel(wb_21s_occ_fld)
  summary(wb_21s_preydiversity) # top ranked model but somewhat questionable coeffs
  summary(wb_21s_anthromort)
  
  
  ####  Wolf - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (wc_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (wc_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_coy_20s_umf, silent = TRUE)) #fail
  (wc_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (wc_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (wc_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_coy_20s_umf, silent = TRUE)) #fail
  (wc_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_coy_20s_umf, silent = TRUE)) #fail
  wc_20s_occ_fld <- fitList(wc_20s_null1, wc_20s_null2, wc_20s_hab, wc_20s_preygroups, #wc_20s_preydiversity,
                            wc_20s_anthromort, wc_20s_anthrodist) #wc_20s_global1, wc_20s_global2
  #' Model selection
  modSel(wc_20s_occ_fld)
  summary(wc_20s_preygroups) 
  
  ####  Wolf - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_coy_21s_umf, silent = TRUE))  
  (wc_21s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_coy_21s_umf, silent = TRUE))  #questionable coeffs without lagomorphs
  (wc_21s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (wc_21s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (wc_21s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_coy_21s_umf, silent = TRUE)) #questionable coeffs
  wc_21s_occ_fld <- fitList(wc_21s_null1, wc_21s_null2, wc_21s_hab, wc_21s_preygroups, wc_21s_preydiversity,
                            wc_21s_anthromort, wc_21s_anthrodist, wc_21s_global1, wc_21s_global2) 
  #' Model selection
  modSel(wc_21s_occ_fld)
  summary(wc_21s_global2)
  summary(wc_21s_preygroups) # global models ranked better but all had questionable coefficient estimates
  
  
  ####  Wolf - Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_lion_20s_umf, silent = TRUE)) #fail
  (wl_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_lion_20s_umf, silent = TRUE)) 
  (wl_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_lion_20s_umf, silent = TRUE))
  wl_20s_occ_fld <- fitList(wl_20s_null1, wl_20s_null2, wl_20s_hab, wl_20s_preygroups, #wl_20s_preydiversity, 
                            wl_20s_anthromort, wl_20s_anthrodist, wl_20s_global1, wl_20s_global2) 
  #' Model selection
  modSel(wl_20s_occ_fld)
  summary(wl_20s_global2)
  
  ####  Wolf - Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, wolf_lion_21s_umf, silent = TRUE)) 
  (wl_21s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, wolf_lion_21s_umf, silent = TRUE)) 
  (wl_21s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, wolf_lion_21s_umf, silent = TRUE)) #questionable coeffs
  (wl_21s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, wolf_lion_21s_umf, silent = TRUE)) #fail
  (wl_21s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, wolf_lion_21s_umf, silent = TRUE)) #fail
  wl_21s_occ_fld <- fitList(wl_21s_null1, wl_21s_null2, wl_21s_hab, wl_21s_preygroups, wl_21s_preydiversity, 
                            wl_21s_anthromort, wl_21s_anthrodist)#, wl_21s_global1, wl_21s_global2) 
  #' Model selection
  modSel(wl_21s_occ_fld)
  summary(wl_21s_anthromort)  # anthrodist model ranked better but had questionable coefficient estimates
  
  
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
  (lbr_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, lion_bear_20s_umf, silent = TRUE)) 
  (lbr_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, lion_bear_20s_umf, silent = TRUE)) 
  (lbr_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, lion_bear_20s_umf, silent = TRUE)) #questionable coeffs
  lbr_20s_occ_fld <- fitList(lbr_20s_null1, lbr_20s_null2, lbr_20s_hab, lbr_20s_preygroups, lbr_20s_preydiversity, 
                             lbr_20s_anthromort, lbr_20s_anthrodist, lbr_20s_global1, lbr_20s_global2) 
  #' Model selection
  modSel(lbr_20s_occ_fld)
  summary(lbr_20s_anthromort)
  
  ####  Lion - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bear_21s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (lbr_21s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, lion_bear_21s_umf, silent = TRUE)) #questionabl coeffs
  (lbr_21s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, lion_bear_21s_umf, silent = TRUE)) #questionable coeffs
  (lbr_21s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, lion_bear_21s_umf, silent = TRUE)) 
  (lbr_21s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, lion_bear_21s_umf, silent = TRUE)) #questionable coeffs
  (lbr_21s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, lion_bear_21s_umf, silent = TRUE)) #fail
  lbr_21s_occ_fld <- fitList(lbr_21s_null1, lbr_21s_null2, lbr_21s_hab, lbr_21s_preygroups, lbr_21s_preydiversity, 
                             lbr_21s_anthromort, lbr_21s_anthrodist, lbr_21s_global1) #, lbr_21s_global2   
  #' Model selection
  modSel(lbr_21s_occ_fld)
  summary(lbr_21s_hab)
  
  
  ####  Lion - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bob_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (lb_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, lion_bob_20s_umf, silent = TRUE)) #questionable coeffs
  (lb_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, lion_bob_20s_umf, silent = TRUE)) 
  (lb_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, lion_bob_20s_umf, silent = TRUE)) #questionable coeffs
  (lb_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, lion_bob_20s_umf, silent = TRUE)) #fail
  (lb_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, lion_bob_20s_umf, silent = TRUE)) #fail
  lb_20s_occ_fld <- fitList(lb_20s_null1, lb_20s_null2, lb_20s_hab, lb_20s_preygroups, lb_20s_preydiversity, 
                             lb_20s_anthromort, lb_20s_anthrodist) #, lb_20s_global1, lb_20s_global2 
  #' Model selection
  modSel(lb_20s_occ_fld)
  summary(lb_20s_anthrodist) # but really this and prey diversity are questionable - next best is null1
  
  ####  Lion - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bob_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (lb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, lion_bob_21s_umf, silent = TRUE)) #fail
  (lb_21s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, lion_bob_21s_umf, silent = TRUE)) #fail
  lb_21s_occ_fld <- fitList(lb_21s_null1, lb_21s_null2, lb_21s_hab, lb_21s_preygroups, lb_21s_preydiversity, 
                             lb_21s_anthromort, lb_21s_global1)   #lb_21s_anthrodist, lb_21s_global2 
  #' Model selection
  modSel(lb_21s_occ_fld)
  summary(lb_21s_hab) # but nothing is actually significant other than intercepts?!?!?!
  
  
  ####  Lion - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_coy_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (lc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_coy_20s_umf, silent = TRUE)) #fail
  (lc_20s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, lion_coy_20s_umf, silent = TRUE)) 
  (lc_20s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, lion_coy_20s_umf, silent = TRUE)) 
  (lc_20s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, lion_coy_20s_umf, silent = TRUE)) 
  (lc_20s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, lion_coy_20s_umf, silent = TRUE)) # questionable coeffs
  lc_20s_occ_fld <- fitList(lc_20s_null1, lc_20s_null2, lc_20s_hab, lc_20s_preygroups, #lc_20s_preydiversity, 
                            lc_20s_anthromort, lc_20s_anthrodist, lc_20s_global1, lc_20s_global2) 
  #' Model selection
  modSel(lc_20s_occ_fld)
  summary(lc_20s_global2) 
  
  ####  Lion - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_coy_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (lc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_coy_21s_umf, silent = TRUE)) 
  (lc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_coy_21s_umf, silent = TRUE)) 
  (lc_21s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, lion_coy_21s_umf, silent = TRUE)) 
  (lc_21s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, lion_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (lc_21s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, lion_coy_21s_umf, silent = TRUE)) #questionable coeffs 
  (lc_21s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, lion_coy_21s_umf, silent = TRUE)) #questionable coeffs
  lc_21s_occ_fld <- fitList(lc_21s_null1, lc_21s_null2, lc_21s_hab, lc_21s_preygroups, lc_21s_preydiversity, 
                            lc_21s_anthromort, lc_21s_anthrodist, lc_21s_global1, lc_21s_global2)    
  #' Model selection
  modSel(lc_21s_occ_fld)
  summary(lc_21s_preygroups) # global models actually ranked better but questionable coeffs
  
  
  ####  Bear - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_bob_20s_umf)  # setup
  
  #'  Run occupancy models using best supported detection sub-model
  (brb_20s_null1 <- occuMulti(detFormulas_setup, occFormulas_null1, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_null2 <- occuMulti(detFormulas_setup, occFormulas_null2, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_hab <- occuMulti(detFormulas_setup, occFormulas_hab, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_preygroups <- occuMulti(detFormulas_setup, occFormulas_preygroups, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_preydiversity <- occuMulti(detFormulas_setup, occFormulas_preydiversity, bear_bob_20s_umf, silent = TRUE)) 
  (brb_20s_anthromort <- occuMulti(detFormulas_setup, occFormulas_anthromort, bear_bob_20s_umf, silent = TRUE)) 
  (brb_20s_anthrodist <- occuMulti(detFormulas_setup, occFormulas_anthrodist, bear_bob_20s_umf, silent = TRUE)) #questionable coeffs
  (brb_20s_global1 <- occuMulti(detFormulas_setup, occFormula_global1, bear_bob_20s_umf, silent = TRUE)) #questionable coeffs
  (brb_20s_global2 <- occuMulti(detFormulas_setup, occFormula_global2, bear_bob_20s_umf, silent = TRUE)) #questionable coeffs
  brb_20s_occ_fld <- fitList(brb_20s_null1, brb_20s_null2, brb_20s_hab, brb_20s_preygroups, brb_20s_preydiversity, 
                             brb_20s_anthromort, brb_20s_anthrodist, brb_20s_global1, brb_20s_global2) 
  #' Model selection
  modSel(brb_20s_occ_fld)
  summary(brb_20s_hab) # anthrodist ranked higher but questionable coeffs
  
  ####  Bear - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_bob_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (brb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, bear_bob_21s_umf, silent = TRUE)) #fail
  brb_21s_occ_fld <- fitList(brb_21s_null1, brb_21s_null2, brb_21s_hab, brb_21s_preygroups, brb_21s_preydiversity, 
                             brb_21s_anthromort, brb_21s_anthrodist, brb_21s_global1) #brb_21s_global2   
  #' Model selection
  modSel(brb_21s_occ_fld)
  summary(brb_21s_preygroups) 
  
  
  ####  Bear - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_coy_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (brc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bear_coy_20s_umf, silent = TRUE)) #fail
  (brc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bear_coy_20s_umf, silent = TRUE)) #fail
  (brc_20s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, bear_coy_20s_umf, silent = TRUE)) 
  (brc_20s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, bear_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (brc_20s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, bear_coy_20s_umf, silent = TRUE)) #fail
  (brc_20s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, bear_coy_20s_umf, silent = TRUE)) #fail
  brc_20s_occ_fld <- fitList(brc_20s_null1, brc_20s_null2, brc_20s_hab, #brc_20s_preygroups, brc_20s_preydiversity, 
                             brc_20s_anthromort, brc_20s_anthrodist) #, brc_20s_global1, brc_20s_global2) 
  #' Model selection
  modSel(brc_20s_occ_fld)
  summary(brc_20s_anthromort) # anthrodist actually ranked better but questionable coeffs
  
  ####  Bear - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_coy_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (brc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bear_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (brc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bear_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (brc_21s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, bear_coy_21s_umf, silent = TRUE)) 
  (brc_21s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, bear_coy_21s_umf, silent = TRUE)) #questionable coeffs
  (brc_21s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, bear_coy_21s_umf, silent = TRUE)) # fail
  (brc_21s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, bear_coy_21s_umf, silent = TRUE)) #questionable coeffs
  brc_21s_occ_fld <- fitList(brc_21s_null1, brc_21s_null2, brc_21s_hab, brc_21s_preygroups, brc_21s_preydiversity, 
                             brc_21s_anthromort, brc_21s_anthrodist, brc_21s_global2) #brc_21s_global1, 
  #' Model selection
  modSel(brc_21s_occ_fld)
  summary(brc_21s_global2) # but definitely questionable coeffs
  
  
  ####  Bobcat - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bob_coy_20s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (bc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bob_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (bc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bob_coy_20s_umf, silent = TRUE)) #questionable coeffs
  (bc_20s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, bob_coy_20s_umf, silent = TRUE)) #fail
  (bc_20s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, bob_coy_20s_umf, silent = TRUE)) # fail
  (bc_20s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, bob_coy_20s_umf, silent = TRUE)) # fail
  bc_20s_occ_fld <- fitList(bc_20s_null1, bc_20s_null2, bc_20s_hab, bc_20s_preygroups, bc_20s_preydiversity, 
                            bc_20s_anthrodist) #, bc_20s_anthromort, bc_20s_global1, bc_20s_global2) 
  #' Model selection
  modSel(bc_20s_occ_fld)
  summary(bc_20s_preydiversity) 
  
  ####  Bear - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bob_coy_21s_umf)  # trail
  
  #'  Run occupancy models using best supported detection sub-model
  (bc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_anthromort <- occuMulti(detFormulas_trail, occFormulas_anthromort, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_anthrodist <- occuMulti(detFormulas_trail, occFormulas_anthrodist, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_global1 <- occuMulti(detFormulas_trail, occFormula_global1, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_global2 <- occuMulti(detFormulas_trail, occFormula_global2, bob_coy_21s_umf, silent = TRUE)) 
  bc_21s_occ_fld <- fitList(bc_21s_null1, bc_21s_null2, bc_21s_hab, bc_21s_preygroups, bc_21s_preydiversity, 
                            bc_21s_anthromort, bc_21s_anthrodist, bc_21s_global1, bc_21s_global2) 
  #' Model selection
  modSel(bc_21s_occ_fld)
  summary(bc_21s_global2) 
  
  
  
  
  
  
  
  
  
  
