  #'  ------------------------------------------------------------------------
  #'  Inputs for for d-Sep iterations when t-1 --> t-1: top-down interference FINAL
  #'  Sept 2026
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
  #'  This list of inputs differs from the input structure in d_Sep_active_regressions_topdown_inter_final.R
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
  
  dSep_iterations_topdown_inter_tmin1_only_final <- list(
    # Regression #1
    # bs_topdown_inter_final[[1]] # "deerHarv.tmin1" "wtd.tmin1"
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #2
    # bs_topdown_inter_final[[2]] # "deerHarv.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #4
    # bs_topdown_inter_final[[4]] # "deerHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #5
    # bs_topdown_inter_final[[5]] # "deerHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #6
    # bs_topdown_inter_final[[6]] # "deerHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #9
    # bs_topdown_inter_final[[9]] # "deerHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #11
    # bs_topdown_inter_final[[11]] # "deerHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #18
    # bs_topdown_inter_final[[18]] # "wtd.tmin1" "moose.tmin1"
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #19
    # bs_topdown_inter_final[[19]] # "wtd.tmin1" "elkHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #20
    # bs_topdown_inter_final[[20]] # "wtd.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #21
    # bs_topdown_inter_final[[21]] # "wtd.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #22
    # bs_topdown_inter_final[[22]] # "wtd.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #23
    # bs_topdown_inter_final[[23]] # "wtd.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #24
    # bs_topdown_inter_final[[24]] # "wtd.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #25
    # bs_topdown_inter_final[[25]] # "wtd.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #26
    # bs_topdown_inter_final[[26]] # "wtd.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #27
    # bs_topdown_inter_final[[27]] # "wtd.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
    # Regression #34
    # bs_topdown_inter_final[[34]] # "moose.tmin1" "elkHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #35
    # bs_topdown_inter_final[[35]] # "moose.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
    # Regression #36
    # bs_topdown_inter_final[[36]] # "moose.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
    # Regression #37
    # bs_topdown_inter_final[[37]] # "moose.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
    # Regression #38
    # bs_topdown_inter_final[[38]] # "moose.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #39
    # bs_topdown_inter_final[[39]] # "moose.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #40
    # bs_topdown_inter_final[[40]] # "moose.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
    # Regression #42
    # bs_topdown_inter_final[[42]] # "moose.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #43
    # bs_topdown_inter_final[[43]] # "moose.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
    # Regression #49
    # bs_topdown_inter_final[[49]] # "elkHarv.tmin1" "elk.tmin1"
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #50
    # bs_topdown_inter_final[[50]] # "elkHarv.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #51
    # bs_topdown_inter_final[[51]] # "elkHarv.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #54
    # bs_topdown_inter_final[[54]] # "elkHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #57
    # bs_topdown_inter_final[[57]] # "elkHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #63
    # bs_topdown_inter_final[[63]] # "elk.tmin1" "lion.tmin1"
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
    # Regression #64
    # bs_topdown_inter_final[[64]] # "elk.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
    # Regression #65
    # bs_topdown_inter_final[[65]] # "elk.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #66
    # bs_topdown_inter_final[[66]] # "elk.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #67
    # bs_topdown_inter_final[[67]] # "elk.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
    # Regression #69
    # bs_topdown_inter_final[[69]] # "elk.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #70
    # bs_topdown_inter_final[[70]] # "elk.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
    # Regression #76
    # bs_topdown_inter_final[[76]] # "lion.tmin1" "coy.tmin1"
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
    # Regression #77
    # bs_topdown_inter_final[[77]] # "lion.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #78
    # bs_topdown_inter_final[[78]] # "lion.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #79
    # bs_topdown_inter_final[[79]] # "lion.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
    # Regression #80
    # bs_topdown_inter_final[[80]] # "lion.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #81
    # bs_topdown_inter_final[[81]] # "lion.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
    # Regression #86
    # bs_topdown_inter_final[[86]] # "coy.tmin1" "bearHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #87
    # bs_topdown_inter_final[[87]] # "coy.tmin1" "wolfHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #88
    # bs_topdown_inter_final[[88]] # "coy.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$coy.hat),
    # Regression #89
    # bs_topdown_inter_final[[89]] # "coy.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #90
    # bs_topdown_inter_final[[90]] # "coy.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$coy.hat),
    # Regression #97
    # bs_topdown_inter_final[[97]] # "bearHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #100
    # bs_topdown_inter_final[[100]] # "bearHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #106
    # bs_topdown_inter_final[[106]] # "wolfHarv.tmin1" "bear.tmin1"
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #109
    # bs_topdown_inter_final[[109]] # "wolfHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
    # Regression #115
    # bs_topdown_inter_final[[115]] # "bear.tmin1" "lionHarv.tmin1"
    # NOTE: FLIPPED so the response is the species, not harvest
    list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
    # Regression #116
    # bs_topdown_inter_final[[116]] # "bear.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$bear.hat),
    # Regression #128
    # bs_topdown_inter_final[[128]] # "lionHarv.tmin1" "wolf.tmin1"
    list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv)
  )
    
    
    
    
    
    
  #   # Regression #1
  #   # bs_topdown_inter_final[[1]] # "wtd.tmin1"   "moose.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #2
  #   # bs_topdown_inter_final[[2]] # "wtd.tmin1" "elk.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #3
  #   # bs_topdown_inter_final[[3]] # "wtd.tmin1"  "lion.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #4
  #   # bs_topdown_inter_final[[4]] # "wtd.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #5               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[5]] # "wtd.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #6
  #   # bs_topdown_inter_final[[6]] # "wtd.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #7               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[7]] # "wtd.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #8               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[8]] # "wtd.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("wtd"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #9
  #   # bs_topdown_inter_final[[9]] # "wtd.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$wtd.hat),
  #   # Regression #16
  #   # bs_topdown_inter_final[[16]] # "moose.tmin1" "elk.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
  #   # Regression #17
  #   # bs_topdown_inter_final[[17]] # "moose.tmin1" "lion.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
  #   # Regression #18
  #   # bs_topdown_inter_final[[18]] # "moose.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
  #   # Regression #19               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[19]] # "moose.tmin1"    "bearHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #20
  #   # bs_topdown_inter_final[[20]] # "moose.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
  #   # Regression #22               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[22]] # "moose.tmin1"    "wolfHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #23               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[23]] # "moose.tmin1"    "lionHarv.tmin1"
  #   list(spp = c("moose"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #24
  #   # bs_topdown_inter_final[[24]] # "moose.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$moose.hat),
  #   # Regression #30
  #   # bs_topdown_inter_final[[30]] # "elk.tmin1"  "lion.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
  #   # Regression #31
  #   # bs_topdown_inter_final[[31]] # "elk.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #32               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[32]] # "elk.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #33
  #   # bs_topdown_inter_final[[33]] # "elk.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
  #   # Regression #35               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[35]] # "elk.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #36               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[36]] # "elk.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("elk"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #37
  #   # bs_topdown_inter_final[[37]] # "elk.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$elk.hat),
  #   # Regression #43
  #   # bs_topdown_inter_final[[43]] # "lion.tmin1" "coy.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
  #   # Regression #44               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[44]] # "lion.tmin1"     "bearHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #45
  #   # bs_topdown_inter_final[[45]] # "lion.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
  #   # Regression #46               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[46]] # "lion.tmin1"     "wolfHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #47               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[47]] # "lion.tmin1"     "lionHarv.tmin1"
  #   list(spp = c("lion"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #48
  #   # bs_topdown_inter_final[[48]] # "lion.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$lion.hat),
  #   # Regression #53               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[53]] # "coy.tmin1"      "bearHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #54
  #   # bs_topdown_inter_final[[54]] # "coy.tmin1"  "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$coy.hat),
  #   # Regression #55               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[55]] # "coy.tmin1"      "wolfHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #56               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[56]] # "coy.tmin1"      "lionHarv.tmin1"
  #   list(spp = c("coy"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #57
  #   # bs_topdown_inter_final[[57]] # "coy.tmin1"  "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$coy.hat),
  #   # Regression #63
  #   # bs_topdown_inter_final[[63]] # "bearHarv.tmin1" "bear.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #67
  #   # bs_topdown_inter_final[[67]] # "bearHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #73               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[73]] # "bear.tmin1"     "wolfHarv.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #74               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # bs_topdown_inter_final[[74]] # "bear.tmin1"     "lionHarv.tmin1"
  #   list(spp = c("bear"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv),
  #   # Regression #75
  #   # bs_topdown_inter_final[[75]] # "bear.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$bear.hat),
  #   # Regression #91
  #   # bs_topdown_inter_final[[91]] # "wolfHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$wolfHarv),
  #   # Regression #97
  #   # bs_topdown_inter_final[[97]] # "lionHarv.tmin1" "wolf.tmin1"
  #   list(spp = c("wolf"), covariate_array = data_JAGS_bundle_topdown_inter_final$lionHarv)
  # )
