  #'  ------------------------------------------------------------------------
  #'  Inputs for for d-Sep iterations when t-1 --> t-1: top-down interference
  #'  Aug 2026
  #'  ------------------------------------------------------------------------
  #'  List inputs for active regressions to iterate over for each d-separation 
  #'  test when variable A t-1 --> variable B t-1.
  #'  
  #'  Input for each iteration are: 
  #'  spp = the response variable in the independence claim (must be the same
  #'    how the spp is represented in each spp.latent estimate)
  #'  covariate_array = the input covariate data used as the focal explanatory 
  #'    variable in the independence claim
  #'  
  #'  This list of inputs differs from the input structure in d_Sep_active_regressions_topdown_inter.R
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
  #'  ------------------------------------------------------------------------
  
  dSep_iterations_topdown_inter_tmin1_only <- list(
    # Regression #1
    # bs_topdown_inter_harv[[1]] # "deerHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #2
    # bs_topdown_inter_harv[[2]] # "deerHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #4
    # bs_topdown_inter_harv[[4]] # "deerHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #5
    # bs_topdown_inter_harv[[5]] # "deerHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #6
    # bs_topdown_inter_harv[[6]] # "deerHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #8
    # bs_topdown_inter_harv[[8]] # "deerHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #11
    # bs_topdown_inter_harv[[11]] # "deerHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$deerHarv),
    # Regression #18
    # bs_topdown_inter_harv[[18]] # "wtd.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #19
    # bs_topdown_inter_harv[[19]] # "wtd.tmin1" "elkHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #20
    # bs_topdown_inter_harv[[20]] # "wtd.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #21
    # bs_topdown_inter_harv[[21]] # "wtd.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #22
    # bs_topdown_inter_harv[[22]] # "wtd.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #23
    # bs_topdown_inter_harv[[23]] # "wtd.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #24
    # bs_topdown_inter_harv[[24]] # "wtd.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #25
    # bs_topdown_inter_harv[[25]] # "wtd.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #26
    # bs_topdown_inter_harv[[26]] # "wtd.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #27
    # bs_topdown_inter_harv[[27]] # "wtd.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
    # Regression #34
    # bs_topdown_inter_harv[[34]] # "moose.tmin1" "elkHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #35
    # bs_topdown_inter_harv[[35]] # "moose.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
    # Regression #36
    # bs_topdown_inter_harv[[36]] # "moose.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
    # Regression #38
    # bs_topdown_inter_harv[[38]] # "moose.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
    # Regression #39
    # bs_topdown_inter_harv[[39]] # "moose.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #40
    # bs_topdown_inter_harv[[40]] # "moose.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
    # Regression #41
    # bs_topdown_inter_harv[[41]] # "moose.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #42
    # bs_topdown_inter_harv[[42]] # "moose.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #43
    # bs_topdown_inter_harv[[43]] # "moose.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
    # Regression #49
    # bs_topdown_inter_harv[[49]] # "elkHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #50
    # bs_topdown_inter_harv[[50]] # "elkHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #52
    # bs_topdown_inter_harv[[52]] # "elkHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #54
    # bs_topdown_inter_harv[[54]] # "elkHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #57
    # bs_topdown_inter_harv[[57]] # "elkHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$elkHarv),
    # Regression #63
    # bs_topdown_inter_harv[[63]] # "elk.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
    # Regression #65
    # bs_topdown_inter_harv[[65]] # "elk.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
    # Regression #66
    # bs_topdown_inter_harv[[66]] # "elk.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #67
    # bs_topdown_inter_harv[[67]] # "elk.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
    # Regression #68
    # bs_topdown_inter_harv[[68]] # "elk.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #69
    # bs_topdown_inter_harv[[69]] # "elk.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #70
    # bs_topdown_inter_harv[[70]] # "elk.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
    # Regression #76
    # bs_topdown_inter_harv[[76]] # "lion.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
    # Regression #77
    # bs_topdown_inter_harv[[77]] # "lion.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #78
    # bs_topdown_inter_harv[[78]] # "lion.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
    # Regression #79
    # bs_topdown_inter_harv[[79]] # "lion.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #80
    # bs_topdown_inter_harv[[80]] # "lion.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #81
    # bs_topdown_inter_harv[[81]] # "lion.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
    # Regression #98
    # bs_topdown_inter_harv[[98]] # "coy.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #99
    # bs_topdown_inter_harv[[99]] # "coy.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$coy.hat),
    # Regression #100
    # bs_topdown_inter_harv[[100]] # "coy.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #101
    # bs_topdown_inter_harv[[101]] # "coy.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #102
    # bs_topdown_inter_harv[[102]] # "coy.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$coy.hat),
    # Regression #108
    # bs_topdown_inter_harv[[108]] # "bearHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #111
    # bs_topdown_inter_harv[[111]] # "bearHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
    # Regression #117
    # bs_topdown_inter_harv[[117]] # "bear.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #118
    # bs_topdown_inter_harv[[118]] # "bear.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
    # Regression #119
    # bs_topdown_inter_harv[[119]] # "bear.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$bear.hat),
    # Regression #126
    # bs_topdown_inter_harv[[126]] # "wolfHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
    # Regression #132
    # bs_topdown_inter_harv[[132]] # "lionHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv))
    
    
    
    
    
  #   # Regression #1
  #   # bs_topdown_inter[[1]] # "wtd.tmin1"   "moose.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # Regression #2
  #   # bs_topdown_inter[[2]] # "wtd.tmin1" "elk.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # Regression #3
  #   # bs_topdown_inter[[3]] # "wtd.tmin1"  "lion.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # Regression #4
  #   # bs_topdown_inter[[4]] # "wtd.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # #########
  #   # Regression #5         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[5]] # "wtd.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # Regression #6
  #   # bs_topdown_inter[[6]] # "wtd.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # #########
  #   # Regression #7         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[7]] # "wtd.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #8         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[8]] # "wtd.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #9
  #   # bs_topdown_inter[[9]] # "wtd.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$wtd.hat),
  #   # Regression #16
  #   # bs_topdown_inter[[16]] # "moose.tmin1" "elk.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
  #   # Regression #17
  #   # bs_topdown_inter[[17]] # "moose.tmin1" "lion.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
  #   # Regression #19
  #   # bs_topdown_inter[[19]] # "moose.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
  #   # #########
  #   # Regression #20         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[20]] # "moose.tmin1"    "bearHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # Regression #21
  #   # bs_topdown_inter[[21]] # "moose.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
  #   # #########
  #   # Regression #22         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[22]] # "moose.tmin1"    "wolfHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #23         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[23]] # "moose.tmin1"    "lionHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #24
  #   # bs_topdown_inter[[24]] # "moose.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$moose.hat),
  #   # Regression #32
  #   # bs_topdown_inter[[32]] # "elk.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
  #   # #########
  #   # Regression #33         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[33]] # "elk.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # Regression #34
  #   # bs_topdown_inter[[34]] # "elk.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
  #   # #########
  #   # Regression #35         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[35]] # "elk.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #36         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[36]] # "elk.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #37
  #   # bs_topdown_inter[[37]] # "elk.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$elk.hat),
  #   # # Regression #43
  #   # bs_topdown_inter[[43]] # "lion.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
  #   # #########
  #   # Regression #44         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[44]] # "lion.tmin1"     "bearHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # Regression #45
  #   # bs_topdown_inter[[45]] # "lion.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
  #   # #########
  #   # Regression #46         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[46]] # "lion.tmin1"     "wolfHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #47         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[47]] # "lion.tmin1"     "lionHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #48
  #   # bs_topdown_inter[[48]] # "lion.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$lion.hat),
  #   # #########
  #   # Regression #65         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[65]] # "coy.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # Regression #66
  #   # bs_topdown_inter[[66]] # "coy.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$coy.hat),
  #   # #########
  #   # Regression #67         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[67]] # "coy.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #68         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[68]] # "coy.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #69
  #   # bs_topdown_inter[[69]] # "coy.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$coy.hat),
  #   # Regression #75
  #   # bs_topdown_inter[[75]] # "bearHarv.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # # Regression #76
  #   # # bs_topdown_inter[[76]] # "bearHarv.tmin1" "wolfHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # # Regression #77
  #   # # bs_topdown_inter[[77]] # "bearHarv.tmin1" "lionHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # Regression #78
  #   # bs_topdown_inter[[78]] # "bearHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$bearHarv),
  #   # #########
  #   # Regression #84         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[84]] # "bear.tmin1"     "wolfHarv.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # #########
  #   # Regression #85         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter[[85]] # "bear.tmin1"     "lionHarv.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv),
  #   # Regression #86
  #   # bs_topdown_inter[[86]] # "bear.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$bear.hat),
  #   # # Regression #92
  #   # # bs_topdown_inter[[92]] # "wolfHarv.tmin1" "lionHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # Regression #93
  #   # bs_topdown_inter[[93]] # "wolfHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$wolfHarv),
  #   # Regression #99
  #   # bs_topdown_inter[[99]] # "lionHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter$lionHarv)
  # )
  