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
  library(stringr)
  
  #'  Source script that formats detection and covariate data for unmarked
  # source("./Scripts/MultiSpp_OccMod/Format_data_2spp_occmod_unmarked.R")
  source("./Scripts/MultiSpp_OccMod/Format_data_2spp_occmod_unmarked_predatorcams_only.R")
  
  
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
  detFormulas_effort <- c("~effort", "~effort")
  
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
  
  #'  Anthropogenic factors
  occFormulas_anthro <- c("~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  #'  Anthropogenic factors hypothesized to be associated with mortality risk
  occFormulas_anthromort <- c("~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock", "~Elev + PercForest + Dist2Burbs + Nlivestock")
  
  #'  Anthropogenic factors hypothesized to be associated with disturbance
  #'  Log distance to nearest road b/c expect behavioral response to dampen with growing distance from road
  #'  Not logging distance to suburbia b/v expect a relatively linear decline in mortality risk further from suburban/urban areas
  occFormulas_anthrodist <- c("~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman", "~Elev + PercForest + logNearestRd + Nhuman")
  
  #'  Global models (1 = prey functional groups; 2 = prey species diversity)
  occFormulas_global1 <- c("~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nbig_deer + Nsmall_deer + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  occFormulas_global2 <- c("~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock", "~Elev + PercForest + Nelk + Nmoose + Nmd + Nwtd + Nlagomorph + Dist2Burbs + logNearestRd + Nhuman + Nlivestock")
  
  
  ####  Detection Sub-Model Selection V1  ####
  #'  Function run and compare detection sub-models for each species pairing
  choose_det_submod_wolf <- function(det_submod, umf) {
    print(det_trail <- occuMulti(det_submod[[1]], occFormulas_null2, umf, silent = TRUE))
    # print(det_setup <- occuMulti(det_submod[[2]], occFormulas_null2, umf, silent = TRUE))
    print(det_height <- occuMulti(det_submod[[3]], occFormulas_null2, umf, silent = TRUE))
    print(det_effort <- occuMulti(det_submod[[4]], occFormulas_null2, umf, silent = TRUE))
    #' List of fitted models
    det_fld <- fitList(det_trail, det_height, det_effort) #det_setup, 
    
    #' Model selection
    print(modSel(det_fld))
  }
  det_submod <- list(detFormulas_trail, detFormulas_setup, detFormulas_height, detFormulas_effort)
  
  
  ####  Wolf - Bear Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_20s_umf)  # effort (P only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_20s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, wolf_bear_20s_umf, silent = TRUE))
  (wbr_20s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, wolf_bear_20s_umf, silent = TRUE)) #, method = "L-BFGS-B", control = list(maxit = 10000)
  (wbr_20s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, wolf_bear_20s_umf, silent = TRUE))
  wbr_20s_occ_fld <- fitList(wbr_20s_null1, wbr_20s_null2, wbr_20s_hab, wbr_20s_preygroups, wbr_20s_preydiversity, 
                             wbr_20s_anthro, wbr_20s_global1, wbr_20s_global2)  
  #' Model selection
  modSel(wbr_20s_occ_fld)
  summary(wbr_20s_hab)
  wbr_20s_top <- wbr_20s_hab

  ####  Wolf - Bear Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bear_21s_umf)  # effort(P only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wbr_21s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, wolf_bear_21s_umf, silent = TRUE))
  (wbr_21s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, wolf_bear_21s_umf, silent = TRUE, method = "L-BFGS-B")) 
  (wbr_21s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, wolf_bear_21s_umf, silent = TRUE, method = "L-BFGS-B")) 
  (wbr_21s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, wolf_bear_21s_umf, silent = TRUE)) 
  (wbr_21s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, wolf_bear_21s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) 
  (wbr_21s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, wolf_bear_21s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) # fails with P-only data
  wbr_21s_occ_fld <- fitList(wbr_21s_null1, wbr_21s_null2, wbr_21s_hab, wbr_21s_preygroups, wbr_21s_preydiversity,  
                             wbr_21s_anthro, wbr_21s_global1)  #, wbr_21s_global2 
  #' Model selection
  modSel(wbr_21s_occ_fld)
  summary(wbr_21s_anthro) 
  wbr_21s_top <- wbr_21s_anthro
  
  ####  Wolf - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_20s_umf)  # height (P only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_20s_null1 <- occuMulti(detFormulas_height, occFormulas_null1, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_null2 <- occuMulti(detFormulas_height, occFormulas_null2, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_hab <- occuMulti(detFormulas_height, occFormulas_hab, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_preygroups <- occuMulti(detFormulas_height, occFormulas_preygroups, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_preydiversity <- occuMulti(detFormulas_height, occFormulas_preydiversity, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_anthro <- occuMulti(detFormulas_height, occFormulas_anthro, wolf_bob_20s_umf, silent = TRUE))
  (wb_20s_global1 <- occuMulti(detFormulas_height, occFormulas_global1, wolf_bob_20s_umf, silent = TRUE)) # P/U: , method = "L-BFGS-B", control = list(maxit = 10000)
  (wb_20s_global2 <- occuMulti(detFormulas_height, occFormulas_global2, wolf_bob_20s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) # fails with P-only data
  wb_20s_occ_fld <- fitList(wb_20s_null1, wb_20s_null2, wb_20s_hab, wb_20s_preygroups, wb_20s_preydiversity,
                            wb_20s_anthro, wb_20s_global1) #, wb_20s_global2
  #' Model selection
  modSel(wb_20s_occ_fld)
  summary(wb_20s_global1)  
  wb_20s_top <- wb_20s_global1  #wb_20s_global2 (P/U)
  
  ####  Wolf - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_bob_21s_umf)  # trail (P only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, wolf_bob_21s_umf, silent = TRUE)) #, starts=c(-1, 1, 0.5, 0.5, 0.5, 2.5, 0.5, 0.5, 0.5, 0.5, -0.1, 2.5, 0.5, -1, 0.5, -0.5, -1, -3.1, -2.5, -2.5,-2, -2.5)
  (wb_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, wolf_bob_21s_umf, silent = TRUE)) 
  (wb_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_bob_21s_umf, silent = TRUE))
  (wb_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, wolf_bob_21s_umf, silent = TRUE, control = list(maxit = 10000))) 
  (wb_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, wolf_bob_21s_umf, silent = TRUE)) #questionable intercepts/coeffs
  wb_21s_occ_fld <- fitList(wb_21s_null1, wb_21s_null2, wb_21s_hab, wb_21s_preygroups, wb_21s_preydiversity, 
                            wb_21s_anthro, wb_21s_global1, wb_21s_global2)  
  #' Model selection
  modSel(wb_21s_occ_fld)
  summary(wb_21s_preydiversity) 
  wb_21s_top <- wb_21s_preydiversity  # wb_21s_global1 (P/U)
  
  
  ####  Wolf - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_20s_umf)  # trail (p-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_20s_umf, silent = TRUE))
  (wc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_coy_20s_umf, silent = TRUE)) #extreme int/coeffs but providing start values doesn't change coeffs #, starts = c(-3, 2, 2.5, -1, -0.5, -1, 3, -2, -2.5, -2, 1.5, 1.5, -2, 1.5, 1.5)
  (wc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, wolf_coy_20s_umf, silent = TRUE)) #fails even with starting values, maxit, and different optimizer
  (wc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, wolf_coy_20s_umf, silent = TRUE)) #fails even with starting values, maxit, and different optimizer
  (wc_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_coy_20s_umf, silent = TRUE)) #, starts = c(-2.5, 2, 1, 1.5, 1, -2, 1, 1, 0.1, -1, -1, -1, 2, 1, 2, -2, 1, -1.5, -1.5, 2, -1, -2, 1.5, 1, -1, 1.2, 1)
  (wc_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, wolf_coy_20s_umf, silent = TRUE)) #fail even with starting values, maxit, and different optimizer
  (wc_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, wolf_coy_20s_umf, silent = TRUE)) #fail even with starting values, maxit, and different optimizer
  wc_20s_occ_fld <- fitList(wc_20s_null1, wc_20s_null2, wc_20s_hab, wc_20s_anthro) 
                            #wc_20s_preygroups, wc_20s_preydiversity, wc_20s_global1, wc_20s_global2 
  #' Model selection
  modSel(wc_20s_occ_fld)
  summary(wc_20s_anthro) 
  wc_20s_top <- wc_20s_anthro
  
  ####  Wolf - Coyote Summer 2021  ####     
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_coy_21s_umf)  # trail (P/U, P only)
  
  #'  Run occupancy models using best supported detection sub-model
  (wc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_coy_21s_umf, silent = TRUE))
  (wc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, wolf_coy_21s_umf, silent = TRUE))  
  (wc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, wolf_coy_21s_umf, silent = TRUE)) 
  (wc_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_coy_21s_umf, silent = TRUE))  #questionable coeffs with P/U data  #, starts = c(-1.3, 0.8, -.7, -.3, -.5, 1, 1.5, 1.2, -.7, -.5, 1.1, -.8, 2.1, 2.0, 1, .9, .9, .53, .9, -1.3, -1, -2.5, 1.9, 0.7, -.6, .8, 1.1)
  (wc_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, wolf_coy_21s_umf, silent = TRUE)) #questionable coeffs with P only & P/U
  (wc_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, wolf_coy_21s_umf, silent = TRUE)) #fail with P only data #questionable coeffs with P/U
  wc_21s_occ_fld <- fitList(wc_21s_null1, wc_21s_null2, wc_21s_hab, wc_21s_preygroups, wc_21s_preydiversity,
                            wc_21s_anthro, wc_21s_global1)  #, wc_21s_global2
  #' Model selection
  modSel(wc_21s_occ_fld)
  summary(wc_21s_global1)
  wc_21s_top <- wc_21s_global1 #wc_21s_global2 with P/U
  
  
  ####  Wolf - Lion Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_20s_umf)  # trail but nothing's signif; effort fails (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, wolf_lion_20s_umf, silent = TRUE))
  (wl_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, wolf_lion_20s_umf, silent = TRUE, method = "L-BFGS-B"))
  (wl_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, wolf_lion_20s_umf, silent = TRUE)) #fails even with different optimizer/maxit/new start values
  (wl_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, wolf_lion_20s_umf, silent = TRUE)) 
  (wl_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, wolf_lion_20s_umf, silent = TRUE, method = "L-BFGS-B"))
  (wl_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, wolf_lion_20s_umf, silent = TRUE))
  wl_20s_occ_fld <- fitList(wl_20s_null1, wl_20s_null2, wl_20s_hab, wl_20s_preygroups, #wl_20s_preydiversity, 
                            wl_20s_anthro, wl_20s_global1, wl_20s_global2) 
  #' Model selection
  modSel(wl_20s_occ_fld)
  summary(wl_20s_hab)
  wl_20s_top <- wl_20s_hab #wl_20s_global2 with P/U data
  
  ####  Wolf - Lion Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod_wolf(det_submod, umf = wolf_lion_21s_umf)  # effort (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (wl_21s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, wolf_lion_21s_umf, silent = TRUE))
  (wl_21s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, wolf_lion_21s_umf, silent = TRUE)) 
  (wl_21s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, wolf_lion_21s_umf, silent = TRUE)) 
  (wl_21s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, wolf_lion_21s_umf, silent = TRUE)) #questionable coeffs with P/U data
  (wl_21s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, wolf_lion_21s_umf, silent = TRUE)) #questionable coeffs with P/U data
  (wl_21s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, wolf_lion_21s_umf, silent = TRUE)) #fails with P-only data; questionable coeffs with P/U data
  wl_21s_occ_fld <- fitList(wl_21s_null1, wl_21s_null2, wl_21s_hab, wl_21s_preygroups, wl_21s_preydiversity, 
                            wl_21s_anthro, wl_21s_global1) #, wl_21s_global2
  #' Model selection
  modSel(wl_21s_occ_fld)
  summary(wl_21s_anthro)  
  wl_21s_top <- wl_21s_anthro
  
  
  ####  Detection Sub-Model Selection V2  ####
  #'  Function run and compare detection sub-models for each species pairing
  choose_det_submod <- function(det_submod, umf) {
    print(det_trail <- occuMulti(det_submod[[1]], occFormulas_null2, umf, silent = TRUE))
    # print(det_setup <- occuMulti(det_submod[[2]], occFormulas_null2, umf, silent = TRUE))
    print(det_height <- occuMulti(det_submod[[3]], occFormulas_null2, umf, silent = TRUE))
    print(det_wolfact <- occuMulti(det_submod[[4]], occFormulas_null2, umf, silent = TRUE))
    print(det_effort <- occuMulti(det_submod[[5]], occFormulas_null2, umf, silent = TRUE))
    #' List of fitted models
    det_fld <- fitList(det_trail, det_height, det_wolfact, det_effort) #det_setup, 
    #' Model selection
    print(modSel(det_fld))
  }
  det_submod <- list(detFormulas_trail, detFormulas_setup, detFormulas_height, detFormulas_wolfact, detFormulas_effort)
  
  
  ####  Lion - Bear Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bear_20s_umf)  # trail but nothing significant, effort fails (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lbr_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_bear_20s_umf, silent = TRUE))
  (lbr_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_bear_20s_umf, silent = TRUE)) 
  (lbr_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, lion_bear_20s_umf, silent = TRUE)) 
  (lbr_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, lion_bear_20s_umf, silent = TRUE)) 
  (lbr_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, lion_bear_20s_umf, silent = TRUE)) #fails with P-only data; questionable coeffs with P/U data
  lbr_20s_occ_fld <- fitList(lbr_20s_null1, lbr_20s_null2, lbr_20s_hab, lbr_20s_preygroups, lbr_20s_preydiversity, 
                             lbr_20s_anthro, lbr_20s_global1) #, lbr_20s_global2
  #' Model selection
  modSel(lbr_20s_occ_fld)
  summary(lbr_20s_anthro)  # hab is next best model at deltaAIC 3.22
  lbr_20s_top <- lbr_20s_anthro #lbr_20s_global2 with P/U data but questionable coeffs (hab is next best supported)
  
  ####  Lion - Bear Summer 2021  ####           
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bear_21s_umf)  # effort (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lbr_21s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, lion_bear_21s_umf, silent = TRUE))
  (lbr_21s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, lion_bear_21s_umf, silent = TRUE, method = "L-BFGS-B")) 
  (lbr_21s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, lion_bear_21s_umf, silent = TRUE, control = list(maxit = 10000))) 
  (lbr_21s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, lion_bear_21s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) #fails with P-only data
  (lbr_21s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, lion_bear_21s_umf, silent = TRUE, control = list(maxit = 10000))) # fails with P-only data
  lbr_21s_occ_fld <- fitList(lbr_21s_null1, lbr_21s_null2, lbr_21s_hab, lbr_21s_preygroups, lbr_21s_preydiversity, 
                             lbr_21s_anthro)    #lbr_21s_global1, lbr_21s_global2
  #' Model selection
  modSel(lbr_21s_occ_fld)
  summary(lbr_21s_hab)
  lbr_21s_top <- lbr_21s_hab
  
  
  ####  Lion - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bob_20s_umf)  # height, effort fails (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lb_20s_null1 <- occuMulti(detFormulas_height, occFormulas_null1, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_null2 <- occuMulti(detFormulas_height, occFormulas_null2, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_hab <- occuMulti(detFormulas_height, occFormulas_hab, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_preygroups <- occuMulti(detFormulas_height, occFormulas_preygroups, lion_bob_20s_umf, silent = TRUE))
  (lb_20s_preydiversity <- occuMulti(detFormulas_height, occFormulas_preydiversity, lion_bob_20s_umf, silent = TRUE)) 
  (lb_20s_anthro <- occuMulti(detFormulas_height, occFormulas_anthro, lion_bob_20s_umf, silent = TRUE)) #questionable coeffs with P-only; , starts = c(1, .5, .5, 1, 0, 2.4, 1, 1, -.5, 1, 0, -1, 2.5, 1, -.5, .5, -1, -.5, .5, -2.3, -1, -1, -2, -1.4, -2.5)
  (lb_20s_global1 <- occuMulti(detFormulas_height, occFormulas_global1, lion_bob_20s_umf, silent = TRUE)) #fails even with different optimizers/maxit
  (lb_20s_global2 <- occuMulti(detFormulas_height, occFormulas_global2, lion_bob_20s_umf, silent = TRUE)) #fails even with different optimizers/maxit
  lb_20s_occ_fld <- fitList(lb_20s_null1, lb_20s_null2, lb_20s_hab, lb_20s_preygroups, lb_20s_preydiversity, 
                            lb_20s_anthro) #, lb_20s_global1, lb_20s_global2 
  #' Model selection
  modSel(lb_20s_occ_fld)
  summary(lb_20s_null1) # lb_20s_anthro technically best but coeffs are questionable & null1 (deltaAIC 0.55)
  lb_20s_top <- lb_20s_null1 #lb_20s_anthro technically best but coeffs are questionable & null1 is next best with P/U
  
  ####  Lion - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_bob_21s_umf)  # trail ( P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lb_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_bob_21s_umf, silent = TRUE))
  (lb_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_bob_21s_umf, silent = TRUE)) #fails with P-only data
  (lb_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, lion_bob_21s_umf, silent = TRUE)) #fails with P/U data even with different optimizer/maxit 
  (lb_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, lion_bob_21s_umf, silent = TRUE)) 
  (lb_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, lion_bob_21s_umf, silent = TRUE)) #fails even with different optimizer/maxit
  lb_21s_occ_fld <- fitList(lb_21s_null1, lb_21s_null2, lb_21s_hab, lb_21s_preygroups, #lb_21s_preydiversity, 
                            lb_21s_anthro, lb_21s_global1)   #lb_21s_global2 
  #' Model selection
  modSel(lb_21s_occ_fld)
  summary(lb_21s_hab)     # but nothing is significant...
  lb_21s_top <- lb_21s_hab
  
  
  ####  Lion - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_coy_20s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_coy_20s_umf, silent = TRUE))
  (lc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_coy_20s_umf, silent = TRUE)) #fails with P/U data even with different optimizers, maxit, and starting values
  (lc_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, lion_coy_20s_umf, silent = TRUE)) 
  (lc_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, lion_coy_20s_umf, silent = TRUE)) 
  (lc_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, lion_coy_20s_umf, silent = TRUE))
  lc_20s_occ_fld <- fitList(lc_20s_null1, lc_20s_null2, lc_20s_hab, lc_20s_preygroups, lc_20s_preydiversity, 
                            lc_20s_anthro, lc_20s_global1, lc_20s_global2) 
  modSel(lc_20s_occ_fld)
  summary(lc_20s_global2)
  lc_20s_top <- lc_20s_global2
  
  ####  Lion - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = lion_coy_21s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (lc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, lion_coy_21s_umf, silent = TRUE))
  (lc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, lion_coy_21s_umf, silent = TRUE)) 
  (lc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, lion_coy_21s_umf, silent = TRUE)) # questionable coeffs with P-only data
  (lc_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, lion_coy_21s_umf, silent = TRUE)) # very questionable coeffs with P-only data
  (lc_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, lion_coy_21s_umf, silent = TRUE)) # very questionable coeffs with P-only data
  (lc_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, lion_coy_21s_umf, silent = TRUE, control = list(maxit = 10000))) # very questionable coeffs with P-only data
  lc_21s_occ_fld <- fitList(lc_21s_null1, lc_21s_null2, lc_21s_hab, lc_21s_preygroups, lc_21s_preydiversity, 
                            lc_21s_anthro, lc_21s_global1, lc_21s_global2)    
  #' Model selection
  modSel(lc_21s_occ_fld)
  summary(lc_21s_global1) 
  lc_21s_top <- lc_21s_global1
  
  
  ####  Bear - Bobcat Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_bob_20s_umf)  # effort (P-only) setup (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (brb_20s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, bear_bob_20s_umf, silent = TRUE))
  (brb_20s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, bear_bob_20s_umf, silent = TRUE)) 
  (brb_20s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, bear_bob_20s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) 
  (brb_20s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, bear_bob_20s_umf, silent = TRUE)) 
  (brb_20s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, bear_bob_20s_umf, silent = TRUE)) #fails with P-only data
  brb_20s_occ_fld <- fitList(brb_20s_null1, brb_20s_null2, brb_20s_hab, brb_20s_preygroups, brb_20s_preydiversity, 
                             brb_20s_anthro, brb_20s_global1) #, brb_20s_global2
  #' Model selection
  modSel(brb_20s_occ_fld)
  summary(brb_20s_hab)    
  brb_20s_top <- brb_20s_hab
  
  ####  Bear - Bobcat Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_bob_21s_umf)  # effort (P-only) trail (P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (brb_21s_null1 <- occuMulti(detFormulas_effort, occFormulas_null1, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_null2 <- occuMulti(detFormulas_effort, occFormulas_null2, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_hab <- occuMulti(detFormulas_effort, occFormulas_hab, bear_bob_21s_umf, silent = TRUE))
  (brb_21s_preygroups <- occuMulti(detFormulas_effort, occFormulas_preygroups, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_preydiversity <- occuMulti(detFormulas_effort, occFormulas_preydiversity, bear_bob_21s_umf, silent = TRUE, method = "L-BFGS-B")) 
  (brb_21s_anthro <- occuMulti(detFormulas_effort, occFormulas_anthro, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_global1 <- occuMulti(detFormulas_effort, occFormulas_global1, bear_bob_21s_umf, silent = TRUE)) 
  (brb_21s_global2 <- occuMulti(detFormulas_effort, occFormulas_global2, bear_bob_21s_umf, silent = TRUE, method = "L-BFGS-B", control = list(maxit = 10000))) 
  brb_21s_occ_fld <- fitList(brb_21s_null1, brb_21s_null2, brb_21s_hab, brb_21s_preygroups, brb_21s_preydiversity, 
                             brb_21s_anthro, brb_21s_global1, brb_21s_global2)   
  #' Model selection
  modSel(brb_21s_occ_fld)
  summary(brb_21s_hab) 
  brb_21s_top <- brb_21s_hab #brb_21s_preygroups with P/U data
  
  
  ####  Bear - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_coy_20s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (brc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bear_coy_20s_umf, silent = TRUE))
  (brc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bear_coy_20s_umf, silent = TRUE)) #fails with P/U data even with different optimizer/maxit
  (brc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bear_coy_20s_umf, silent = TRUE)) #fails with P/U data event with different optimizer/maxit
  (brc_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, bear_coy_20s_umf, silent = TRUE)) # very large questionable coeffs with P/U data , starts = c(2.5, 2, 2, -1, -.5, 2.5, -1, 2, 1, -.5, -1, -1, 2.5, 1, -2,-1,-1,1,1,-2,1,-1.5,-.5,-.5,1)
  (brc_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, bear_coy_20s_umf, silent = TRUE)) #fails even with different optimizer/maxit
  (brc_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, bear_coy_20s_umf, silent = TRUE)) #fails with P/U data even with different optimizer/maxit
  brc_20s_occ_fld <- fitList(brc_20s_null1, brc_20s_null2, brc_20s_hab, brc_20s_preygroups, brc_20s_preydiversity, 
                             brc_20s_anthro, brc_20s_global2) #, brc_20s_global1
  #' Model selection
  modSel(brc_20s_occ_fld)
  summary(brc_20s_anthro) 
  brc_20s_top <- brc_20s_anthro
  
  ####  Bear - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bear_coy_21s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (brc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bear_coy_21s_umf, silent = TRUE))
  (brc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bear_coy_21s_umf, silent = TRUE)) # fails with P-only data 
  (brc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bear_coy_21s_umf, silent = TRUE)) # fails with P-only data
  (brc_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, bear_coy_21s_umf, silent = TRUE)) # questionable coeffs
  (brc_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, bear_coy_21s_umf, silent = TRUE)) # fails with P-only data; questionable coeffs with P/U data
  (brc_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, bear_coy_21s_umf, silent = TRUE)) # fails with P-only data
  brc_21s_occ_fld <- fitList(brc_21s_null1, brc_21s_null2, brc_21s_hab, #brc_21s_preygroups, brc_21s_preydiversity, 
                             brc_21s_anthro) #, brc_21s_global1, brc_21s_global2
  #' Model selection
  modSel(brc_21s_occ_fld)
  summary(brc_21s_anthro)
  brc_21s_top <- brc_21s_anthro# brc_21s_global1 with P/U but questionable coeffs
  
  
  ####  Bobcat - Coyote Summer 2020  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bob_coy_20s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (bc_20s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bob_coy_20s_umf, silent = TRUE))
  (bc_20s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bob_coy_20s_umf, silent = TRUE)) 
  (bc_20s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bob_coy_20s_umf, silent = TRUE)) 
  (bc_20s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, bob_coy_20s_umf, silent = TRUE)) # fails even with different optimizer/maxit
  (bc_20s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, bob_coy_20s_umf, silent = TRUE)) # fails even with different optimizer/maxit
  (bc_20s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, bob_coy_20s_umf, silent = TRUE)) # fails even with different optimizer/maxit
  bc_20s_occ_fld <- fitList(bc_20s_null1, bc_20s_null2, bc_20s_hab, bc_20s_preygroups, bc_20s_preydiversity) 
                            #bc_20s_anthro, bc_20s_global1, bc_20s_global2) 
  #' Model selection
  modSel(bc_20s_occ_fld)
  summary(bc_20s_preydiversity)  
  bc_20s_top <- bc_20s_preydiversity
  
  ####  Bobcat - Coyote Summer 2021  ####
  #'  Review detection sub-models
  choose_det_submod(det_submod, umf = bob_coy_21s_umf)  # trail (P-only & P/U)
  
  #'  Run occupancy models using best supported detection sub-model
  (bc_21s_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_hab <- occuMulti(detFormulas_trail, occFormulas_hab, bob_coy_21s_umf, silent = TRUE))
  (bc_21s_preygroups <- occuMulti(detFormulas_trail, occFormulas_preygroups, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_preydiversity <- occuMulti(detFormulas_trail, occFormulas_preydiversity, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_anthro <- occuMulti(detFormulas_trail, occFormulas_anthro, bob_coy_21s_umf, silent = TRUE)) #fails even with different optimizer/maxit
  (bc_21s_global1 <- occuMulti(detFormulas_trail, occFormulas_global1, bob_coy_21s_umf, silent = TRUE)) 
  (bc_21s_global2 <- occuMulti(detFormulas_trail, occFormulas_global2, bob_coy_21s_umf, silent = TRUE)) 
  bc_21s_occ_fld <- fitList(bc_21s_null1, bc_21s_null2, bc_21s_hab, bc_21s_preygroups, bc_21s_preydiversity, 
                            bc_21s_global1, bc_21s_global2) #bc_21s_anthro,   
  #' Model selection
  modSel(bc_21s_occ_fld)
  summary(bc_21s_global1)   
  bc_21s_top <- bc_21s_global1
  
  
  #' Save model outputs in one giant R image
  save.image(file = "./Outputs/MultiSpp_OccMod_Outputs/MultiSpp_CoOcc_Models_2spp_PredatorCamsOnly.RData")
  
  
  #'  ------------------
  ####  Summary tables  ####
  #'  ------------------
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model, appends species/season info,
  #'  re-formats data into wide format for better comparison across species
  #'  NOTE: ignore warning about computing `Parameter = ifelse(...)`
  
  #'  Function to save occupancy results
  rounddig <- 2
  occ_out <- function(mod, spp1, spp2, season) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Season = rep(season, nrow(.)),
        Estimate = round(Estimate, rounddig),
        SE = round(SE, rounddig),
        Pval = round(`P(>|z|)`, 3)
      ) %>%
      dplyr::select(-`P(>|z|)`) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species1, .before = Parameter) %>%
      relocate(Species2, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%  
      dplyr::select(-z) %>%
      mutate(
        # SE = round(SE, rounddig),
        SE = paste0("(", SE, ")"),
        #'  Generalize species names (necessary for wide format)
        Parameter = ifelse(grepl(Species1, Parameter), str_replace(Parameter, Species1, "Species 1"), Parameter),
        Parameter = ifelse(grepl(Species2, Parameter), str_replace(Parameter, Species2, "Species 2"), Parameter)
      ) 
    
    return(out)
  }
  
  #'  Run each model through function
  occ_wolfbear_20s <- occ_out(wbr_20s_top, "wolf", "bear", "Summer 2020")
  occ_wolfbear_21s <- occ_out(wbr_21s_top, "wolf", "bear", "Summer 2021")
  occ_wolfbob_20s <- occ_out(wb_20s_top, "wolf", "bobcat", "Summer 2020")
  occ_wolfbob_21s <- occ_out(wb_21s_top, "wolf", "bobcat", "Summer 2021")
  occ_wolfcoy_20s <- occ_out(wc_20s_top, "wolf", "coyote", "Summer 2020")
  occ_wolfcoy_21s <- occ_out(wc_21s_top, "wolf", "coyote", "Summer 2021")
  occ_wolflion_20s <- occ_out(wl_20s_top, "wolf", "lion", "Summer 2020")
  occ_wolflion_21s <- occ_out(wl_21s_top, "wolf", "lion", "Summer 2021")
  occ_lionbear_20s <- occ_out(lbr_20s_top, "lion", "bear", "Summer 2020")
  occ_lionbear_21s <- occ_out(lbr_21s_top, "lion", "bear", "Summer 2021")
  occ_lionbob_20s <- occ_out(lb_20s_top, "lion", "bobcat", "Summer 2020")
  occ_lionbob_21s <- occ_out(lb_21s_top, "lion", "bobcat", "Summer 2021")
  occ_lioncoy_20s <- occ_out(lc_20s_top, "lion", "coyote", "Summer 2020")
  occ_lioncoy_21s <- occ_out(lc_21s_top, "lion", "coyote", "Summer 2021")
  occ_bearbob_20s <- occ_out(brb_20s_top, "bear", "bobcat", "Summer 2020")
  occ_bearbob_21s <- occ_out(brb_21s_top, "bear", "bobcat", "Summer 2021")
  occ_bearcoy_20s <- occ_out(brc_20s_top, "bear", "coyote", "Summer 2020")
  occ_bearcoy_21s <- occ_out(brc_21s_top, "bear", "coyote", "Summer 2021")
  occ_bobcoy_20s <- occ_out(bc_20s_top, "bobcat", "coyote", "Summer 2020")
  occ_bobcoy_21s <- occ_out(bc_21s_top, "bobcat", "coyote", "Summer 2021")
  
  #'  Merge into larger data frames based on best supported model
  #'  Null models
  occ_results_null <- occ_lionbob_20s #lb_20s_null1
  #'  Habitat models
  # occ_results_habitat <- rbind(occ_wolfbear_20s, occ_lionbear_21s, occ_lionbob_21s, occ_bearbob_20s) #wbr_20s_hab, lbr_21s_hab, lb_21s_hab, brb_20s_hab
  occ_results_habitat <- rbind(occ_wolfbear_20s, occ_wolflion_20s, occ_lionbear_21s, occ_lionbob_21s, occ_bearbob_20s, occ_bearbob_21s) #wbr_20s_hab, wl_20s_hab, lbr_21s_hab, lb_21s_hab, brb_20s_hab, brb_21s_hab
  #'  Prey relative abundance by functional group
  # occ_results_preygroup <- occ_bearbob_21s #brb_21s_preygroups
  #'  Prey relative abundance by species
  # occ_results_preydiversity <- occ_bobcoy_20s #bc_20s_preydiversity
  occ_results_preydiversity <- rbind(occ_wolfbob_21s, occ_bobcoy_20s) #wb_21s_preydiversity, bc_20s_preydiversity
  #'  Proxies for anthropogenic risk and disturbance
  # occ_results_anthro <- rbind(occ_wolfbear_21s, occ_wolfcoy_20s, occ_wolflion_21s, occ_lionbob_20s, occ_bearcoy_20s) #wbr_21s_anthro, wc_20s_anthro, wl_21s_anthro, lb_20s_anthro, brc_20s_anthro
  occ_results_anthro <- rbind(occ_wolfbear_21s, occ_wolfcoy_20s, occ_wolflion_21s, occ_lionbear_20s, occ_bearcoy_20s, occ_bearcoy_21s) #wbr_21s_anthro, wc_20s_anthro, wl_21s_anthro, lbr_20s_anthro, brc_20s_anthro, brc_21s_anthro
  #'  Global model 1
  # occ_results_global1 <- rbind(occ_wolfbob_21s, occ_lioncoy_21s, occ_bearcoy_21s, occ_bobcoy_21s) #wb_21s_global1, lc_21s_global1, brc_21s_global1, bc_21s_global1
  occ_results_global1 <- rbind(occ_wolfbob_20s, occ_wolfcoy_21s, occ_lioncoy_21s, occ_bobcoy_21s) #wb_20s_global1, wc_21s_global1, lc_21s_global1, bc_21s_global1
  #'  Global model 2
  # occ_results_global2 <- rbind(occ_wolfbob_20s, occ_wolfcoy_21s, occ_wolflion_20s, occ_lionbear_20s, occ_lioncoy_20s) #wb_20s_global2, wc_21s_global2, wl_20s_global2, lbr_20s_global2, lc_20s_global2
  occ_results_global2 <- rbind(occ_lioncoy_20s) #lc_20s_global2
  #'  Merge all results together
  # occ_results_allmodels <- rbind(occ_wolfbear_20s, occ_lionbear_21s, occ_lionbob_21s, 
  #                                occ_bearbob_20s,occ_bearbob_21s, occ_bobcoy_20s, 
  #                                occ_wolfbear_21s, occ_wolfcoy_20s, occ_wolflion_21s, 
  #                                occ_lionbob_20s, occ_bearcoy_20s, occ_wolfbob_21s, 
  #                                occ_lioncoy_21s, occ_bearcoy_21s, occ_bobcoy_21s,
  #                                occ_wolfbob_20s, occ_wolfcoy_21s, occ_wolflion_20s, 
  #                                occ_lionbear_20s, occ_lioncoy_20s) 
  occ_results_allmodels <- rbind(occ_lionbob_20s, occ_wolfbear_20s, occ_wolflion_20s, 
                                 occ_lionbear_21s, occ_lionbob_21s, occ_bearbob_20s, 
                                 occ_bearbob_21s, occ_wolfbob_21s, occ_bobcoy_20s, 
                                 occ_wolfbear_21s, occ_wolfcoy_20s, occ_wolflion_21s, 
                                 occ_lionbear_20s, occ_bearcoy_20s, occ_bearcoy_21s, 
                                 occ_wolfbob_20s, occ_wolfcoy_21s, occ_lioncoy_21s, 
                                 occ_bobcoy_21s, occ_lioncoy_20s)
  
  
  ####  Spread results into wide format and rename  ####
  #'  Habitat filtering hypothesis results
  results_occmod_wide_habitat <- occ_results_habitat %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  #'  Relative abundance of prey functional groups hypothesis results
  results_occmod_wide_preygroup <- occ_results_preygroup %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] Nlagomorph", .after = "[Species 1] Nsmall_deer") %>%
    relocate("[Species 2] Nlagomorph", .after = "[Species 2] Nsmall_deer") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlagomorph") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Nbig_deer", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] Nsmall_deer", .after = "[Species 1:Species 2] Nbig_deer") %>%
    relocate("[Species 1:Species 2] Nlagomorph", .after = "[Species 1:Species 2] Nsmall_deer") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Nbig_deer", c("[Species 1] Large deer activity (SE)", "[Species 1] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nbig_deer", c("[Species 2] Large deer activity (SE)", "[Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nsmall_deer", c("[Species 1] Small deer activity (SE)", "[Species 1] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nsmall_deer", c("[Species 2] Small deer activity (SE)", "[Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlagomorph", c("[Species 1] Lagomorph activity (SE)", "[Species 1] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlagomorph", c("[Species 2] Lagomorph activity (SE)", "[Species 2] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nbig_deer", c("[Species 1:Species 2] Large deer activity (SE)", "[Species 1:Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nsmall_deer", c("[Species 1:Species 2] Small deer activity (SE)", "[Species 1:Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlagomorph", c("[Species 1:Species 2] Lagomorph activity (SE)", "[Species 1:Species 2] Lagomorph activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  Relative abundance of prey species hypothesis results
  results_occmod_wide_preydiversity <- occ_results_preydiversity %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] Nlagomorph", .after = "[Species 1] Nwtd") %>%
    relocate("[Species 2] Nlagomorph", .after = "[Species 2] Nwtd") %>%
    relocate("[Species 1] Nmd", .after = "[Species 1] Nmoose") %>%
    relocate("[Species 2] Nmd", .after = "[Species 2] Nmoose") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlagomorph") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Nelk", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] Nmoose", .after = "[Species 1:Species 2] Nelk") %>%
    relocate("[Species 1:Species 2] Nmd", .after = "[Species 1:Species 2] Nmoose") %>%
    relocate("[Species 1:Species 2] Nwtd", .after = "[Species 1:Species 2] Nmd") %>%
    relocate("[Species 1:Species 2] Nlagomorph", .after = "[Species 1:Species 2] Nwtd") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Nelk", c("[Species 1] Elk activity (SE)", "[Species 1] Elk activity Pval"), sep = "_") %>%
    separate("[Species 2] Nelk", c("[Species 2] Elk activity (SE)", "[Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmoose", c("[Species 1] Moose activity (SE)", "[Species 1] Moose activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmoose", c("[Species 2] Moose activity (SE)", "[Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmd", c("[Species 1] Mule deer activity (SE)", "[Species 1] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmd", c("[Species 2] Mule deer activity (SE)", "[Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nwtd", c("[Species 1] White-tailed deer activity (SE)", "[Species 1] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nwtd", c("[Species 2] White-tailed deer activity (SE)", "[Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlagomorph", c("[Species 1] Lagomorph activity (SE)", "[Species 1] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlagomorph", c("[Species 2] Lagomorph activity (SE)", "[Species 2] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nelk", c("[Species 1:Species 2] Elk activity (SE)", "[Species 1:Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmoose", c("[Species 1:Species 2] Moose activity (SE)", "[Species 1:Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmd", c("[Species 1:Species 2] Mule deer activity (SE)", "[Species 1:Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nwtd", c("[Species 1:Species 2] White-tailed deer activity (SE)", "[Species 1:Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlagomorph", c("[Species 1:Species 2] Lagomorph activity (SE)", "[Species 1:Species 2] Lagomorph activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  Anthropogenic hypothesis results
  results_occmod_wide_anthro <- occ_results_anthro %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlivestock") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Dist2Burbs", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] logNearestRd", .after = "[Species 1:Species 2] Dist2Burbs") %>%
    relocate("[Species 1:Species 2] Nhuman", .after = "[Species 1:Species 2] logNearestRd") %>%
    relocate("[Species 1:Species 2] Nlivestock", .after = "[Species 1:Species 2] Nhuman") %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Dist2Burbs", c("[Species 1] Distance to Suburban (SE)", "[Species 1] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 2] Dist2Burbs", c("[Species 2] Distance to Suburban (SE)", "[Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1] logNearestRd", c("[Species 1] log(Nearest Road) (SE)", "[Species 1] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 2] logNearestRd", c("[Species 2] log(Nearest Road) (SE)", "[Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1] Nhuman", c("[Species 1] Human activity (SE)", "[Species 1] Human activity Pval"), sep = "_") %>%
    separate("[Species 2] Nhuman", c("[Species 2] Human activity (SE)", "[Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlivestock", c("[Species 1] Livestock activity (SE)", "[Species 1] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlivestock", c("[Species 2] Livestock activity (SE)", "[Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Dist2Burbs", c("[Species 1:Species 2] Distance to Suburban (SE)", "[Species 1:Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] logNearestRd", c("[Species 1:Species 2] log(Nearest Road) (SE)", "[Species 1:Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nhuman", c("[Species 1:Species 2] Human activity (SE)", "[Species 1:Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlivestock", c("[Species 1:Species 2] Livestock activity (SE)", "[Species 1:Species 2] Livestock activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  Global model #1 with relative abundance of prey functional groups
  results_occmod_wide_global1 <- occ_results_global1 %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] Nbig_deer", .after = "[Species 1] Nlivestock") %>%
    relocate("[Species 2] Nbig_deer", .after = "[Species 2] Nlivestock") %>%
    relocate("[Species 1] Nlagomorph", .after = "[Species 1] Nsmall_deer") %>%
    relocate("[Species 2] Nlagomorph", .after = "[Species 2] Nsmall_deer") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlagomorph") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Dist2Burbs", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] logNearestRd", .after = "[Species 1:Species 2] Dist2Burbs") %>%
    relocate("[Species 1:Species 2] Nhuman", .after = "[Species 1:Species 2] logNearestRd") %>%
    relocate("[Species 1:Species 2] Nlivestock", .after = "[Species 1:Species 2] Nhuman") %>%
    relocate("[Species 1:Species 2] Nbig_deer", .after = "[Species 1:Species 2] Nlivestock") %>%
    relocate("[Species 1:Species 2] Nsmall_deer", .after = "[Species 1:Species 2] Nbig_deer") %>%
    relocate("[Species 1:Species 2] Nlagomorph", .after = "[Species 1:Species 2] Nsmall_deer") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Dist2Burbs", c("[Species 1] Distance to Suburban (SE)", "[Species 1] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 2] Dist2Burbs", c("[Species 2] Distance to Suburban (SE)", "[Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1] logNearestRd", c("[Species 1] log(Nearest Road) (SE)", "[Species 1] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 2] logNearestRd", c("[Species 2] log(Nearest Road) (SE)", "[Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1] Nhuman", c("[Species 1] Human activity (SE)", "[Species 1] Human activity Pval"), sep = "_") %>%
    separate("[Species 2] Nhuman", c("[Species 2] Human activity (SE)", "[Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlivestock", c("[Species 1] Livestock activity (SE)", "[Species 1] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlivestock", c("[Species 2] Livestock activity (SE)", "[Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1] Nbig_deer", c("[Species 1] Large deer activity (SE)", "[Species 1] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nbig_deer", c("[Species 2] Large deer activity (SE)", "[Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nsmall_deer", c("[Species 1] Small deer activity (SE)", "[Species 1] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nsmall_deer", c("[Species 2] Small deer activity (SE)", "[Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlagomorph", c("[Species 1] Lagomorph activity (SE)", "[Species 1] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlagomorph", c("[Species 2] Lagomorph activity (SE)", "[Species 2] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Dist2Burbs", c("[Species 1:Species 2] Distance to Suburban (SE)", "[Species 1:Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] logNearestRd", c("[Species 1:Species 2] log(Nearest Road) (SE)", "[Species 1:Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nhuman", c("[Species 1:Species 2] Human activity (SE)", "[Species 1:Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlivestock", c("[Species 1:Species 2] Livestock activity (SE)", "[Species 1:Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nbig_deer", c("[Species 1:Species 2] Large deer activity (SE)", "[Species 1:Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nsmall_deer", c("[Species 1:Species 2] Small deer activity (SE)", "[Species 1:Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlagomorph", c("[Species 1:Species 2] Lagomorph activity (SE)", "[Species 1:Species 2] Lagomorph activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  Global model #2 with relative abundance of individual prey species
  results_occmod_wide_global2 <- occ_results_global2 %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] Nlagomorph", .after = "[Species 1] Nwtd") %>%
    relocate("[Species 2] Nlagomorph", .after = "[Species 2] Nwtd") %>%
    relocate("[Species 1] Nmd", .after = "[Species 1] Nmoose") %>%
    relocate("[Species 2] Nmd", .after = "[Species 2] Nmoose") %>%
    relocate("[Species 1] Nelk", .after = "[Species 1] Nlivestock") %>%
    relocate("[Species 2] Nelk", .after = "[Species 2] Nlivestock") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlagomorph") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Dist2Burbs", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] logNearestRd", .after = "[Species 1:Species 2] Dist2Burbs") %>%
    relocate("[Species 1:Species 2] Nhuman", .after = "[Species 1:Species 2] logNearestRd") %>%
    relocate("[Species 1:Species 2] Nlivestock", .after = "[Species 1:Species 2] Nhuman") %>%
    relocate("[Species 1:Species 2] Nelk", .after = "[Species 1:Species 2] Nlivestock") %>%
    relocate("[Species 1:Species 2] Nmoose", .after = "[Species 1:Species 2] Nelk") %>%
    relocate("[Species 1:Species 2] Nmd", .after = "[Species 1:Species 2] Nmoose") %>%
    relocate("[Species 1:Species 2] Nwtd", .after = "[Species 1:Species 2] Nmd") %>%
    relocate("[Species 1:Species 2] Nlagomorph", .after = "[Species 1:Species 2] Nwtd") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Dist2Burbs", c("[Species 1] Distance to Suburban (SE)", "[Species 1] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 2] Dist2Burbs", c("[Species 2] Distance to Suburban (SE)", "[Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1] logNearestRd", c("[Species 1] log(Nearest Road) (SE)", "[Species 1] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 2] logNearestRd", c("[Species 2] log(Nearest Road) (SE)", "[Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1] Nhuman", c("[Species 1] Human activity (SE)", "[Species 1] Human activity Pval"), sep = "_") %>%
    separate("[Species 2] Nhuman", c("[Species 2] Human activity (SE)", "[Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlivestock", c("[Species 1] Livestock activity (SE)", "[Species 1] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlivestock", c("[Species 2] Livestock activity (SE)", "[Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1] Nelk", c("[Species 1] Elk activity (SE)", "[Species 1] Elk activity Pval"), sep = "_") %>%
    separate("[Species 2] Nelk", c("[Species 2] Elk activity (SE)", "[Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmoose", c("[Species 1] Moose activity (SE)", "[Species 1] Moose activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmoose", c("[Species 2] Moose activity (SE)", "[Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmd", c("[Species 1] Mule deer activity (SE)", "[Species 1] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmd", c("[Species 2] Mule deer activity (SE)", "[Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nwtd", c("[Species 1] White-tailed deer activity (SE)", "[Species 1] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nwtd", c("[Species 2] White-tailed deer activity (SE)", "[Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlagomorph", c("[Species 1] Lagomorph activity (SE)", "[Species 1] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlagomorph", c("[Species 2] Lagomorph activity (SE)", "[Species 2] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Dist2Burbs", c("[Species 1:Species 2] Distance to Suburban (SE)", "[Species 1:Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] logNearestRd", c("[Species 1:Species 2] log(Nearest Road) (SE)", "[Species 1:Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nhuman", c("[Species 1:Species 2] Human activity (SE)", "[Species 1:Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlivestock", c("[Species 1:Species 2] Livestock activity (SE)", "[Species 1:Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nelk", c("[Species 1:Species 2] Elk activity (SE)", "[Species 1:Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmoose", c("[Species 1:Species 2] Moose activity (SE)", "[Species 1:Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmd", c("[Species 1:Species 2] Mule deer activity (SE)", "[Species 1:Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nwtd", c("[Species 1:Species 2] White-tailed deer activity (SE)", "[Species 1:Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlagomorph", c("[Species 1:Species 2] Lagomorph activity (SE)", "[Species 1:Species 2] Lagomorph activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  OH MY GOD ALL OF THEM COMBINED INTO ONE CRAZY PLOT
  results_occmod_wide_allmodels <- occ_results_allmodels %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1] Elev", .after = "[Species 1] (Intercept)") %>%
    relocate("[Species 2] Elev", .after = "[Species 2] (Intercept)") %>%
    relocate("[Species 1] PercForest", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] PercForest", .after = "[Species 2] Elev") %>%
    relocate("[Species 1] Nbig_deer", .after = "[Species 1] Nlivestock") %>%
    relocate("[Species 2] Nbig_deer", .after = "[Species 2] Nlivestock") %>%
    relocate("[Species 1] Nsmall_deer", .after = "[Species 1] Nbig_deer") %>%
    relocate("[Species 2] Nsmall_deer", .after = "[Species 2] Nbig_deer") %>%
    relocate("[Species 1] Nlagomorph", .after = "[Species 1] Nwtd") %>%
    relocate("[Species 2] Nlagomorph", .after = "[Species 2] Nwtd") %>%
    relocate("[Species 1] Nmd", .after = "[Species 1] Nmoose") %>%
    relocate("[Species 2] Nmd", .after = "[Species 2] Nmoose") %>%
    relocate("[Species 1] Nelk", .after = "[Species 1] Nsmall_deer") %>%
    relocate("[Species 2] Nelk", .after = "[Species 2] Nsmall_deer") %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Nlagomorph") %>%
    relocate("[Species 1:Species 2] Elev", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] PercForest", .after = "[Species 1:Species 2] Elev") %>%
    relocate("[Species 1:Species 2] Dist2Burbs", .after = "[Species 1:Species 2] PercForest") %>%
    relocate("[Species 1:Species 2] logNearestRd", .after = "[Species 1:Species 2] Dist2Burbs") %>%
    relocate("[Species 1:Species 2] Nhuman", .after = "[Species 1:Species 2] logNearestRd") %>%
    relocate("[Species 1:Species 2] Nlivestock", .after = "[Species 1:Species 2] Nhuman") %>%
    relocate("[Species 1:Species 2] Nbig_deer", .after = "[Species 1:Species 2] Nlivestock") %>%
    relocate("[Species 1:Species 2] Nsmall_deer", .after = "[Species 1:Species 2] Nbig_deer") %>%
    relocate("[Species 1:Species 2] Nelk", .after = "[Species 1:Species 2] Nsmall_deer") %>%
    relocate("[Species 1:Species 2] Nmoose", .after = "[Species 1:Species 2] Nelk") %>%
    relocate("[Species 1:Species 2] Nmd", .after = "[Species 1:Species 2] Nmoose") %>%
    relocate("[Species 1:Species 2] Nwtd", .after = "[Species 1:Species 2] Nmd") %>%
    relocate("[Species 1:Species 2] Nlagomorph", .after = "[Species 1:Species 2] Nwtd") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] Percent Forest (SE)", "[Species 1] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] Percent Forest (SE)", "[Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1] Dist2Burbs", c("[Species 1] Distance to Suburban (SE)", "[Species 1] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 2] Dist2Burbs", c("[Species 2] Distance to Suburban (SE)", "[Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1] logNearestRd", c("[Species 1] log(Nearest Road) (SE)", "[Species 1] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 2] logNearestRd", c("[Species 2] log(Nearest Road) (SE)", "[Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1] Nhuman", c("[Species 1] Human activity (SE)", "[Species 1] Human activity Pval"), sep = "_") %>%
    separate("[Species 2] Nhuman", c("[Species 2] Human activity (SE)", "[Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlivestock", c("[Species 1] Livestock activity (SE)", "[Species 1] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlivestock", c("[Species 2] Livestock activity (SE)", "[Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1] Nbig_deer", c("[Species 1] Large deer activity (SE)", "[Species 1] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nbig_deer", c("[Species 2] Large deer activity (SE)", "[Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nsmall_deer", c("[Species 1] Small deer activity (SE)", "[Species 1] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nsmall_deer", c("[Species 2] Small deer activity (SE)", "[Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nelk", c("[Species 1] Elk activity (SE)", "[Species 1] Elk activity Pval"), sep = "_") %>%
    separate("[Species 2] Nelk", c("[Species 2] Elk activity (SE)", "[Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmoose", c("[Species 1] Moose activity (SE)", "[Species 1] Moose activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmoose", c("[Species 2] Moose activity (SE)", "[Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1] Nmd", c("[Species 1] Mule deer activity (SE)", "[Species 1] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nmd", c("[Species 2] Mule deer activity (SE)", "[Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nwtd", c("[Species 1] White-tailed deer activity (SE)", "[Species 1] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 2] Nwtd", c("[Species 2] White-tailed deer activity (SE)", "[Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1] Nlagomorph", c("[Species 1] Lagomorph activity (SE)", "[Species 1] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 2] Nlagomorph", c("[Species 2] Lagomorph activity (SE)", "[Species 2] Lagomorph activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Elev", c("[Species 1:Species 2] Elev (SE)", "[Species 1:Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] PercForest", c("[Species 1:Species 2] Percent Forest (SE)", "[Species 1:Species 2] Percent Forest Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Dist2Burbs", c("[Species 1:Species 2] Distance to Suburban (SE)", "[Species 1:Species 2] Distance to Suburban Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] logNearestRd", c("[Species 1:Species 2] log(Nearest Road) (SE)", "[Species 1:Species 2] log(Nearest Road) Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nhuman", c("[Species 1:Species 2] Human activity (SE)", "[Species 1:Species 2] Human activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlivestock", c("[Species 1:Species 2] Livestock activity (SE)", "[Species 1:Species 2] Livestock activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nbig_deer", c("[Species 1:Species 2] Large deer activity (SE)", "[Species 1:Species 2] Large deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nsmall_deer", c("[Species 1:Species 2] Small deer activity (SE)", "[Species 1:Species 2] Small deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nelk", c("[Species 1:Species 2] Elk activity (SE)", "[Species 1:Species 2] Elk activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmoose", c("[Species 1:Species 2] Moose activity (SE)", "[Species 1:Species 2] Moose activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nmd", c("[Species 1:Species 2] Mule deer activity (SE)", "[Species 1:Species 2] Mule deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nwtd", c("[Species 1:Species 2] White-tailed deer activity (SE)", "[Species 1:Species 2] White-tailed deer activity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Nlagomorph", c("[Species 1:Species 2] Lagomorph activity (SE)", "[Species 1:Species 2] Lagomorph activity Pval"), sep = "_") %>%
    arrange(match(Species1, c("bear", "bobcat", "coyote", "lion", "wolf")))
  
  
  #'  Save!
  # write.csv(results_occmod_wide_habitat, paste0("./Outputs/Tables/CoOcc_OccProb_Habitat_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_preygroup, paste0("./Outputs/Tables/CoOcc_OccProb_PreyGroup_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_preydiversity, paste0("./Outputs/Tables/CoOcc_OccProb_PreyDiversity_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_anthro, paste0("./Outputs/Tables/CoOcc_OccProb_Anthro_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_global1, paste0("./Outputs/Tables/CoOcc_OccProb_GlobalPreyGroup_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_global2, paste0("./Outputs/Tables/CoOcc_OccProb_GlobalPreyDiversity_wide_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_allmodels, paste0("./Outputs/Tables/CoOcc_OccProb_AllModels_wide_", Sys.Date(), ".csv"))
  
  write.csv(results_occmod_wide_habitat, paste0("./Outputs/Tables/CoOcc_OccProb_Habitat_wide_Ponly_", Sys.Date(), ".csv"))
  # write.csv(results_occmod_wide_preygroup, paste0("./Outputs/Tables/CoOcc_OccProb_PreyGroup_wide_Ponly_", Sys.Date(), ".csv"))
  write.csv(results_occmod_wide_preydiversity, paste0("./Outputs/Tables/CoOcc_OccProb_PreyDiversity_wide_Ponly_", Sys.Date(), ".csv"))
  write.csv(results_occmod_wide_anthro, paste0("./Outputs/Tables/CoOcc_OccProb_Anthro_wide_Ponly_", Sys.Date(), ".csv"))
  write.csv(results_occmod_wide_global1, paste0("./Outputs/Tables/CoOcc_OccProb_GlobalPreyGroup_wide_Ponly_", Sys.Date(), ".csv"))
  write.csv(results_occmod_wide_global2, paste0("./Outputs/Tables/CoOcc_OccProb_GlobalPreyDiversity_wide_Ponly_", Sys.Date(), ".csv"))
  write.csv(results_occmod_wide_allmodels, paste0("./Outputs/Tables/CoOcc_OccProb_AllModels_wide_Ponly_", Sys.Date(), ".csv"))

  
