  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations when t-1 --> t-1: bottom-up interference FINAL
  #'  Sept 2026
  #'  ----------------------------------------------------------------
  #'  List inputs for active regressions to iterate over for each d-separation 
  #'  test when variable A t-1 --> variable B t-1.
  #'  
  #'  Input for each iteration are: 
  #'  spp = the response variable in the independence claim (must be the same
  #'    how the spp is represented in each spp.latent estimate)
  #'  covariate_array = the input covariate data used as the focal explanatory 
  #'    variable in the independence claim
  #'  
  #'  This list of inputs differs from the input structure in d_Sep_active_regressions_bottomup.R
  #'  because we are not building an entirely new model with each iteration here. 
  #'  Instead, we are pulling the posterior of the latent variable of interest from
  #'  the original fitted SEM and using this as the "observed" data in an updated
  #'  version of the SEM. In other words, the spp.latent output from the original
  #'  formulation of the SEM is used as the spp.hat input data for y in the updated
  #'  SEM for that species, specifically. This was necessary under the situation 
  #'  when both the response and the explanatory variables in an independence claim
  #'  are for the t-1 time step since the original formulation of the SEM only has 
  #'  regressions for response variables for time step t. This is a somewhat hacky
  #'  way to test these independence claims but is the easiest approach given how
  #'  the SEMs are coded in JAGS. As a result, these d-Sep tests essentially assume
  #'  the latent estimates (aka the response variable of interest) are measured
  #'  without error. The uncertainty in the latent posterior estimates are therefore
  #'  not propagated in these d-Sep tests. But this allows us to at least assess
  #'  if there's any evidence of correlation between the variables and then fully
  #'  propagate uncertainty in the final updated version of the SEM.
  #'  
  #'  This script also grabs the input covariate data (e.g., forest) needed
  #'  to test the independence claim for each iteration. The fit_aux_claim() function
  #'  in the d_Sep_test_and_FisherC_for_Bayesian_SEM.R then refits the SEM as it
  #'  was originally parameterized but with the addition of the added covariate.
  #'  ----------------------------------------------------------------
  
  dSep_iterations_bottomup_inter_tmin1_only <- list(
    ########
    # Regression #1
    # bs_bottomup_inter_final[[1]] # "wsi.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    # ########
    # # Regression #2
    # # bs_bottomup_inter_final[[2]] # "wsi.tmin1"    "forest.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),     # exog on exog
    ########
    # Regression #3
    # bs_bottomup_inter_final[[3]] # "wsi.tmin1"   "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    ########
    # Regression #4
    # bs_bottomup_inter_final[[4]] # "wsi.tmin1"  "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    ########
    # Regression #5
    # bs_bottomup_inter_final[[5]] # "wsi.tmin1"  "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    ########
    # Regression #6
    # bs_bottomup_inter_final[[6]] # "wsi.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    ########
    # Regression #8
    # bs_bottomup_inter_final[[8]] # "wsi.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wsi),
    ########
    # Regression #12                # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT FOREST
    # bs_bottomup_inter_final[[12]] # "coy.tmin1"    "forest.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #13
    # bs_bottomup_inter_final[[13]] # "coy.tmin1"   "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_bottomup_inter_final$coy.hat),
    ########
    # Regression #15
    # bs_bottomup_inter_final[[15]] # "coy.tmin1"  "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_bottomup_inter_final$coy.hat),
    ########
    # Regression #16
    # bs_bottomup_inter_final[[16]] # "coy.tmin1"  "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_bottomup_inter_final$coy.hat),
    ########
    # Regression #17
    # bs_bottomup_inter_final[[17]] # "coy.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$coy.hat),
    ########
    # Regression #19
    # bs_bottomup_inter_final[[19]] # "coy.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$coy.hat),
    ########
    # Regression #24
    # bs_bottomup_inter_final[[24]] # "forest.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #25
    # bs_bottomup_inter_final[[25]] # "forest.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #26
    # bs_bottomup_inter_final[[26]] # "forest.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #27
    # bs_bottomup_inter_final[[27]] # "forest.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #29
    # bs_bottomup_inter_final[[29]] # "forest.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$forest),
    ########
    # Regression #32
    # bs_bottomup_inter_final[[32]] # "moose.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_bottomup_inter_final$moose.hat),
    ########
    # Regression #33
    # bs_bottomup_inter_final[[33]] # "moose.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_bottomup_inter_final$moose.hat),
    ########
    # Regression #34
    # bs_bottomup_inter_final[[34]] # "moose.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$moose.hat),
    ########
    # Regression #37
    # bs_bottomup_inter_final[[37]] # "moose.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$moose.hat),
    ########
    # Regression #50
    # bs_bottomup_inter_final[[50]] # "wolf.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wolf.hat),
    ########
    # Regression #51
    # bs_bottomup_inter_final[[51]] # "wolf.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wolf.hat),
    ########
    # Regression #53
    # bs_bottomup_inter_final[[53]] # "wolf.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wolf.hat),
    ########
    # Regression #56
    # bs_bottomup_inter_final[[56]] # "bear.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_bottomup_inter_final$bear.hat),
    ########
    # Regression #58
    # bs_bottomup_inter_final[[58]] # "bear.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$bear.hat),
    ########
    # Regression #61
    # bs_bottomup_inter_final[[61]] # "wtd.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_bottomup_inter_final$wtd.hat)
  )
  
