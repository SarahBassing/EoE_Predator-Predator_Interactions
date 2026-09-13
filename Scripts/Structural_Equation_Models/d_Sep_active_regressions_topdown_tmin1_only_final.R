  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations when t-1 --> t-1: top-down FINAL
  #'  Aug 2026
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
  #'  This list of inputs differs from the input structure in d_Sep_active_regressions_topdown.R
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
  #'  This script also grabs the input covariate data (e.g., deerHarv) needed
  #'  to test the independence claim for each iteration. The fit_aux_claim() function
  #'  in the d_Sep_test_and_FisherC_for_Bayesian_SEM.R then refits the SEM as it
  #'  was originally parameterized but with the addition of the added covariate.
  #'  ----------------------------------------------------------------
  
  dSep_iterations_topdown_tmin1_only_final <- list(
    # Regression #2
    # bs_topdown_final[[2]] # "deerHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #3
    # bs_topdown_final[[3]] # "deerHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #5
    # bs_topdown_final[[5]] # "deerHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #6
    # bs_topdown_final[[6]] # "deerHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #8
    # bs_topdown_final[[8]] # "deerHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #10
    # bs_topdown_final[[10]] # "deerHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #13
    # bs_topdown_final[[13]] # "deerHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$deerHarv),
    # Regression #18
    # bs_topdown_final[[18]] # "elkHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #19
    # bs_topdown_final[[19]] # "elkHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #21
    # bs_topdown_final[[21]] # "elkHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #22
    # bs_topdown_final[[22]] # "elkHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #24
    # bs_topdown_final[[24]] # "elkHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #26
    # bs_topdown_final[[26]] # "elkHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #30
    # bs_topdown_final[[30]] # "elkHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$elkHarv),
    # Regression #34
    # bs_topdown_final[[34]] # "lion.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #35         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[35]] # "lion.tmin1"     "bearHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #36
    # bs_topdown_final[[36]] # "lion.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #37
    # bs_topdown_final[[37]] # "lion.tmin1"  "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #38         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[38]] # "lion.tmin1"     "wolfHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #39
    # bs_topdown_final[[39]] # "lion.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #41
    # bs_topdown_final[[41]] # "lion.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #44
    # bs_topdown_final[[44]] # "lion.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$lion.hat),
    # Regression #46         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[46]] # "lion.tmin1"     "lionHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #48         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[48]] # "coy.tmin1"      "bearHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #49
    # bs_topdown_final[[49]] # "coy.tmin1"  "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$coy.hat),
    # Regression #50
    # bs_topdown_final[[50]] # "coy.tmin1"   "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$coy.hat),
    # Regression #51         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[51]] # "coy.tmin1"      "wolfHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #52
    # bs_topdown_final[[52]] # "coy.tmin1"  "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$coy.hat),
    # Regression #54
    # bs_topdown_final[[54]] # "coy.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$coy.hat),
    # Regression #56
    # bs_topdown_final[[56]] # "coy.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$coy.hat),
    # Regression #59         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[59]] # "coy.tmin1"      "lionHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #61
    # bs_topdown_final[[61]] # "bearHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #62
    # bs_topdown_final[[62]] # "bearHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #64
    # bs_topdown_final[[64]] # "bearHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #66
    # bs_topdown_final[[66]] # "bearHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #69
    # bs_topdown_final[[69]] # "bearHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #74
    # bs_topdown_final[[74]] # "bear.tmin1"  "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$bear.hat),
    # Regression #75         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[75]] # "bear.tmin1"     "wolfHarv.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #76
    # bs_topdown_final[[76]] # "bear.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$bear.hat),
    # Regression #78
    # bs_topdown_final[[78]] # "bear.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$bear.hat),
    # Regression #81
    # bs_topdown_final[[81]] # "bear.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$bear.hat),
    # Regression #83         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[83]] # "bear.tmin1"     "lionHarv.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #85         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[85]] # "moose.tmin1"    "wolfHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #86
    # bs_topdown_final[[86]] # "moose.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$moose.hat),
    # Regression #87
    # bs_topdown_final[[87]] # "moose.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$moose.hat),
    # Regression #91
    # bs_topdown_final[[91]] # "moose.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$moose.hat),
    # Regression #93         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[93]] # "moose.tmin1"    "lionHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #95
    # bs_topdown_final[[95]] # "wolfHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #97
    # bs_topdown_final[[97]] # "wolfHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #101
    # bs_topdown_final[[101]] # "wolfHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$wolfHarv),
    # Regression #105
    # bs_topdown_final[[105]] # "wolf.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$wolf.hat),
    # Regression #108
    # bs_topdown_final[[108]] # "wolf.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$wolf.hat),
    # Regression #109         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[109]] # "wolf.tmin1"     "lionHarv.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #121
    # bs_topdown_final[[121]] # "wtd.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$wtd.hat),
    # Regression #123         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[123]] # "wtd.tmin1"      "lionHarv.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv),
    # Regression #139         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[139]] # "elk.tmin1"      "lionHarv.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_final$lionHarv)
  )
  