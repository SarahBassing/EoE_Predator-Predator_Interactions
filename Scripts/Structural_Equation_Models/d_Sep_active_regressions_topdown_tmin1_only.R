  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations when t-1 --> t-1: top-down
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
  
  dSep_iterations_topdown_tmin1_only <- list(
    # Regression #1
    # bs_topdown[[1]] # "deerHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #2
    # bs_topdown[[2]] # "deerHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #4
    # bs_topdown[[4]] # "deerHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #5
    # bs_topdown[[5]] # "deerHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #6
    # bs_topdown[[6]] # "deerHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #9
    # bs_topdown[[9]] # "deerHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #12
    # bs_topdown[[12]] # "deerHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$deerHarv),
    # Regression #18
    # bs_topdown[[18]] # "wtd.tmin1"   "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    #########
    # Regression #19   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[19]] # "wtd.tmin1"     "elkHarv.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    #########
    # Regression #20
    # bs_topdown[[20]] # "wtd.tmin1"  "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    # Regression #21
    # bs_topdown[[21]] # "wtd.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    # Regression #22
    # bs_topdown[[22]] # "wtd.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    # ##########
    # Regression #24   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[24]] # "wtd.tmin1"      "bearHarv.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ##########
    # Regression #25
    # bs_topdown[[25]] # "wtd.tmin1"  "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    # ##########
    # Regression #27   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[27]] # "wtd.tmin1"      "wolfHarv.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ##########
    # Regression #28
    # bs_topdown[[28]] # "wtd.tmin1"  "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$wtd.hat),
    # ##########
    # Regression #32   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[32]] # "wtd.tmin1"      "lionHarv.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # Regression #34   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[34]] # "moose.tmin1"   "elkHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # ##########
    # Regression #35
    # bs_topdown[[35]] # "moose.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$moose.hat),
    # Regression #37
    # bs_topdown[[37]] # "moose.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$moose.hat),
    # Regression #38
    # bs_topdown[[38]] # "moose.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$moose.hat),
    # ##########
    # Regression #40   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[40]] # "moose.tmin1"    "bearHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ##########
    # Regression #41
    # bs_topdown[[41]] # "moose.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$moose.hat),
    # ##########
    # Regression #43   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[43]] # "moose.tmin1"    "wolfHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ##########
    # Regression #44
    # bs_topdown[[44]] # "moose.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$moose.hat),
    # ##########
    # Regression #47   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[47]] # "moose.tmin1"    "lionHarv.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # ##########
    # Regression #49
    # bs_topdown[[49]] # "elkHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # Regression #51
    # bs_topdown[[51]] # "elkHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # Regression #52
    # bs_topdown[[52]] # "elkHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # Regression #55
    # bs_topdown[[55]]  # "elkHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # Regression #58
    # bs_topdown[[58]]  # "elkHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$elkHarv),
    # Regression #63
    # bs_topdown[[63]]  # "lion.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$lion.hat),
    # Regression #64
    # bs_topdown[[64]]  # "lion.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$lion.hat),
    # ##########
    # Regression #66   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[66]] # "lion.tmin1"     "bearHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ##########
    # Regression #67
    # bs_topdown[[67]]  # "lion.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$lion.hat),
    # ##########
    # Regression #69   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[69]] # "lion.tmin1"     "wolfHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ##########
    # Regression #70
    # bs_topdown[[70]]  # "lion.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$lion.hat),
    # ##########
    # Regression #73    # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[73]]  # "lion.tmin1"     "lionHarv.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # ##########
    # Regression #88
    # bs_topdown[[88]] # "elk.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$elk.hat),
    # ##########
    # Regression #90   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[90]] # "elk.tmin1"      "bearHarv.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ##########
    # Regression #91
    # bs_topdown[[91]] # "elk.tmin1"  "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$elk.hat),
    # ##########
    # Regression #93   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[93]] # "elk.tmin1"  "wolfHarv.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ##########
    # Regression #94
    # bs_topdown[[94]] # "elk.tmin1"  "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$elk.hat),
    # ##########
    # Regression #97   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[97]] # "elk.tmin1"  "lionHarv.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # Regression #99   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[99]] # "coy.tmin1"      "bearHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ##########
    # Regression #100
    # bs_topdown[[100]] # "coy.tmin1"      "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$coy.hat),
    # ##########
    # Regression #102   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[102]] # "coy.tmin1"      "wolfHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ##########
    # Regression #103
    # bs_topdown[[103]] # "coy.tmin1"      "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$coy.hat),
    # ##########
    # Regression #107   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[107]] # "coy.tmin1"      "lionHarv.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # ##########
    # Regression #119
    # bs_topdown[[119]] # "bearHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # Regression #121
    # bs_topdown[[121]] # "bearHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$bearHarv),
    # ############
    # Regression #127   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[127]] # "bear.tmin1"     "wolfHarv.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # ############
    # Regression #128
    # bs_topdown[[128]] # "bear.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$bear.hat),
    # ############
    # Regression #132   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[132]] # "bear.tmin1" "lionHarv.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown$lionHarv),
    # ############
    # Regression #141
    # bs_topdown[[141]] # "wolfHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$wolfHarv),
    # #########
    # Regression #146   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown[[146]] # "wolf.tmin1"     "lionHarv.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown$lionHarv)
    # #########
  )
  