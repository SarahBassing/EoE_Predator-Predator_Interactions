  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down exploitative FINAL
  #'  Sep 2026
  #'  ----------------------------------------------------------------
  #'  List active regressions to iterate over for each d-separation test. 
  #'  Necessary objects for each active regression include:
  #'  
  #'  -- dSep_test: integer that indicates which regression (1:7) within the SEM 
  #'     gets the custom regression where: 
  #'            1 = cougar 
  #'            2 = wolf 
  #'            3 = black bear 
  #'            4 = coyote 
  #'            5 = elk 
  #'            6 = moose
  #'            7 = wtd 
  #'  -- covariates: vector that include the covariates included in the SEM 
  #'     a priori, plus the one needed for the conditional independence test.
  #'     NOTE: the last covariate in this vector should always be the explanatory
  #'     variable you are testing the independence claim on. This will make life
  #'     easier when using the p.rope() and bayes_pvalue() functions in 
  #'     d_Sep_test_and_FishersC_for_Bayesian_SEM.R
  #'     NOTE: in cases where independence claim is testing whether variable A at 
  #'     time t affects variable B at time t-1 (which is not possible given our
  #'     current understanding of space and time), I switched which variable was
  #'     the response vs explanatory variable. This effectively tests whether 
  #'     there is a reasonable correlation between the variables. Although the 
  #'     exact independence claim cannot be included in the SEM (due to how time
  #'     works), it allows us to include the independence claim when calculating
  #'     GoF. This is important for accurately assessing whether the hypothesized 
  #'     causal model fits the data.
  #'  -- spp: vector containing the suffix of each parameter name to be appended
  #'     to "beta" in the model template. Must include a period before the term's
  #'     name to work.
  #'  -- indices: vector containing the index value for each parameter. For most
  #'     regressions, each parameter's name is unique so the index is 1. But for some,
  #'     the parameter name is used twice in the model representing different coefficients
  #'     (e.g., beta.harvest[1] and beta.harvest[2] correspond to the effects of 
  #'     wolfHarv and lionHarv in the same regression). Double check indexing for 
  #'     each regression to ensure appropriate coefficient estimates are saved and 
  #'     used in ROPE / Bayesian p-value d-Sep tests.
  #'  -- lags: vector containing character strings that indicate whether the year
  #'     index should be lagged by one year. "y-1" applies a lag to the variable;
  #'     "y" does not apply a lag.
  #'  
  #'  Some independence claims cannot be assess with this script given how the
  #'  SEM regressions are written. Those are included here for consistency but 
  #'  are commented out and assessed in other scripts: 
  #'  d_Sep_active_regression_topdown_tmin1_only.R  and 
  #'  d_Sep_active_regression_topdown_exog_only.R
  #'  ---------------------------------------------------------------- 
  
  dSep_iterations_topdown_final <- list(
    # # Regression #1
    # # bs_topdown_final[[1]] # "deerHarv.tmin1" "elkHarv.tmin1" 
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #2
    # # bs_topdown_final[[2]] # "deerHarv.tmin1" "lion.tmin1"    
    # list(dSep_test = 1, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # #########
    # # Regression #3
    # # bs_topdown_final[[3]] # "deerHarv.tmin1" "coy.tmin1" 
    # list(dSep_test = 4, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #4
    # # bs_topdown_final[[4]] # "deerHarv.tmin1" "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #5
    # # bs_topdown_final[[5]] # "deerHarv.tmin1" "bear.tmin1"  
    # list(dSep_test = 3, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # #########
    # # Regression #6
    # # bs_topdown_final[[6]] # "deerHarv.tmin1" "moose.tmin1" 
    # list(dSep_test = 6, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #7
    # # bs_topdown_final[[7]] # "deerHarv.tmin1" "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #8
    # # bs_topdown_final[[8]] # "deerHarv.tmin1" "wolf.tmin1"    
    # list(dSep_test = 2, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #9
    # bs_topdown_final[[9]] # "deerHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"  
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "deerHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #10
    # # bs_topdown_final[[10]] # "deerHarv.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #11
    # bs_topdown_final[[11]] # "deerHarv.tmin1" "coy.t"          "coy.tmin1"      "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "deerHarv"), spp = c(".coy", ".wtd", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #12
    # bs_topdown_final[[12]] # "deerHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "deerHarv"), spp = c(".harvest", ".bear", ".wtd", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #13
    # # bs_topdown_final[[13]] # "deerHarv.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #14
    # bs_topdown_final[[14]] # "deerHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "deerHarv"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #15
    # bs_topdown_final[[15]] # "deerHarv.tmin1" "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "deerHarv"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".harvest"), indices = c(1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # # Regression #16
    # # bs_topdown_final[[16]] # "deerHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # Regression #17
    # bs_topdown_final[[17]] # "deerHarv.tmin1" "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "deerHarv"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".harvest"), indices = c(1,2,1,1,2), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #18
    # # bs_topdown_final[[18]] # "elkHarv.tmin1" "lion.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #19
    # # bs_topdown_final[[19]] # "elkHarv.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # # Regression #20
    # # bs_topdown_final[[20]] # "elkHarv.tmin1"  "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #21
    # # bs_topdown_final[[21]] # "elkHarv.tmin1" "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #22
    # # bs_topdown_final[[22]] # "elkHarv.tmin1" "moose.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # # Regression #23
    # # bs_topdown_final[[23]] # "elkHarv.tmin1"  "wolfHarv.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #24
    # # bs_topdown_final[[24]] # "elkHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #25
    # bs_topdown_final[[25]] # "elkHarv.tmin1" "moose.t"       "moose.tmin1"   "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elkHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #26
    # # bs_topdown_final[[26]] # "elkHarv.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #27
    # bs_topdown_final[[27]] # "elkHarv.tmin1" "coy.t"         "coy.tmin1"     "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "elkHarv"), spp = c(".coy", ".wtd", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #28
    # bs_topdown_final[[28]] # "elkHarv.tmin1"  "wtd.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "elkHarv"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #29
    # bs_topdown_final[[29]] # "elkHarv.tmin1"  "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "elkHarv"), spp = c(".harvest", ".bear", ".wtd", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #30
    # # bs_topdown_final[[30]] # "elkHarv.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #31
    # bs_topdown_final[[31]] # "elkHarv.tmin1"  "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "elkHarv"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".harvest"), indices = c(1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # # Regression #32
    # # bs_topdown_final[[32]] # "elkHarv.tmin1"  "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # Regression #33
    # bs_topdown_final[[33]] # "elkHarv.tmin1"  "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "elkHarv"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".harvest"), indices = c(1,2,1,1,2), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #34
    # # bs_topdown_final[[34]] # "lion.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #35         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[35]] # "lion.tmin1"     "bearHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #36
    # # bs_topdown_final[[36]] # "lion.tmin1" "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #37
    # # bs_topdown_final[[37]] # "lion.tmin1"  "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #38         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[38]] # "lion.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #39
    # # bs_topdown_final[[39]] # "lion.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #40
    # bs_topdown_final[[40]] # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #41
    # # bs_topdown_final[[41]] # "lion.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #42
    # bs_topdown_final[[42]] # "lion.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "lion.latent"), spp = c(".coy", ".wtd", ".lion"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #43
    # bs_topdown_final[[43]] # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "lion.latent"), spp = c(".harvest", ".bear", ".wtd", ".lion"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #44
    # # bs_topdown_final[[44]] # "lion.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #45
    # bs_topdown_final[[45]] # "lion.tmin1"     "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "lion.latent"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".lion"), indices = c(1,1,1,1,2,1,1), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # #########
    # # Regression #46         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[46]] # "lion.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #47
    # bs_topdown_final[[47]] # "lion.tmin1"     "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "lion.latent"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".lion"), indices = c(1,2,1,1,1), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #48         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[48]] # "coy.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # # Regression #49
    # # bs_topdown_final[[49]] # "coy.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # # Regression #50
    # # bs_topdown_final[[50]] # "coy.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #51         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[51]] # "coy.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #52
    # # bs_topdown_final[[52]] # "coy.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #53
    # bs_topdown_final[[53]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #54
    # # bs_topdown_final[[54]] # "coy.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #55
    # bs_topdown_final[[55]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "coy.latent"), spp = c(".harvest", ".bear", ".wtd", ".coy"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #56
    # # bs_topdown_final[[56]] # "coy.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #57
    # bs_topdown_final[[57]] # "coy.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "bear.tmin1"    "wolf.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "coy.latent"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".coy"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #58
    # bs_topdown_final[[58]] # "coy.tmin1"      "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "coy.latent"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".coy"), indices = c(1,1,1,1,2,1,1), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # #########
    # # Regression #59         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[59]] # "coy.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #60
    # bs_topdown_final[[60]] # "coy.tmin1"      "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "coy.latent"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".coy"), indices = c(1,2,1,1,1), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #61
    # # bs_topdown_final[[61]] # "bearHarv.tmin1" "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #62
    # # bs_topdown_final[[62]] # "bearHarv.tmin1" "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # # Regression #63
    # # bs_topdown_final[[63]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # #########
    # # Regression #64
    # # bs_topdown_final[[64]] # "bearHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #65
    # bs_topdown_final[[65]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #66
    # # bs_topdown_final[[66]] # "bearHarv.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #67
    # bs_topdown_final[[67]] # "bearHarv.tmin1" "coy.t"          "coy.tmin1"      "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "bearHarv"), spp = c(".coy", ".wtd", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #68
    # bs_topdown_final[[68]] # "bearHarv.tmin1" "wtd.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "bearHarv"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #69
    # # bs_topdown_final[[69]] # "bearHarv.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #70
    # bs_topdown_final[[70]] # "bearHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "bearHarv"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #71
    # bs_topdown_final[[71]] # "bearHarv.tmin1" "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "bearHarv"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".harvest"), indices = c(1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # # Regression #72
    # # bs_topdown_final[[72]] # "bearHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # Regression #73
    # bs_topdown_final[[73]] # "bearHarv.tmin1" "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "bearHarv"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".harvest"), indices = c(1,2,1,1,2), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #74
    # # bs_topdown_final[[74]] # "bear.tmin1"  "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #75         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[75]] # "bear.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #76
    # # bs_topdown_final[[76]] # "bear.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #77
    # bs_topdown_final[[77]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #78
    # # bs_topdown_final[[78]] # "bear.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #79
    # bs_topdown_final[[79]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "bear.latent"), spp = c(".coy", ".wtd", ".bear"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #80
    # bs_topdown_final[[80]] # "bear.tmin1"     "wtd.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "bear.latent"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".bear"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #81
    # # bs_topdown_final[[81]] # "bear.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #82
    # bs_topdown_final[[82]] # "bear.tmin1"     "wolf.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "bear.latent"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".bear"), indices = c(1,1,1,1,2,1,1), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # #########
    # # Regression #83         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[83]] # "bear.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #84
    # bs_topdown_final[[84]] # "bear.tmin1"     "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "bear.latent"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".bear"), indices = c(1,2,1,1,1), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #85         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[85]] # "moose.tmin1"    "wolfHarv.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #86
    # # bs_topdown_final[[86]] # "moose.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #87
    # # bs_topdown_final[[87]] # "moose.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #88
    # bs_topdown_final[[88]] # "moose.tmin1" "coy.t"       "coy.tmin1"   "wtd.tmin1"
    list(dSep_test = 3, covariates = c("coy.latent", "wtd.latent", "moose.latent"), spp = c(".coy", ".wtd", ".moose"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #89
    # bs_topdown_final[[89]] # "moose.tmin1"    "wtd.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "moose.latent"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".moose"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #90
    # bs_topdown_final[[90]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "moose.latent"), spp = c(".harvest", ".bear", ".wtd", ".moose"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #91
    # # bs_topdown_final[[91]] # "moose.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #92
    # bs_topdown_final[[92]] # "moose.tmin1"   "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "bear.tmin1"    "wolf.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "moose.latent"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".moose"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #93         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[93]] # "moose.tmin1"    "lionHarv.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #94
    # bs_topdown_final[[94]] # "moose.tmin1"    "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "moose.latent"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".moose"), indices = c(1,2,1,1,1), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #95
    # # bs_topdown_final[[95]] # "wolfHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #96
    # bs_topdown_final[[96]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #97
    # # bs_topdown_final[[97]] # "wolfHarv.tmin1" "wtd.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #98
    # bs_topdown_final[[98]] # "wolfHarv.tmin1" "coy.t"          "coy.tmin1"      "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wolfHarv"), spp = c(".coy", ".wtd", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #99
    # bs_topdown_final[[99]] # "wolfHarv.tmin1" "wtd.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "wolfHarv"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #100
    # bs_topdown_final[[100]] # "wolfHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".wtd", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #101
    # # bs_topdown_final[[101]] # "wolfHarv.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #102
    # bs_topdown_final[[102]] # "wolfHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "wolfHarv"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #103
    # # bs_topdown_final[[103]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),   # exog on exog
    # Regression #104
    # bs_topdown_final[[104]] # "wolfHarv.tmin1" "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "wolfHarv"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".harvest"), indices = c(1,2,1,1,2), lags = c("y-1","y","y-1","y-1","y-1")),
    # #########
    # # Regression #105
    # # bs_topdown_final[[105]] # "wolf.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #106
    # bs_topdown_final[[106]] # "wolf.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wolf.latent"), spp = c(".coy", ".wtd", ".wolf"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #107
    # bs_topdown_final[[107]] # "wolf.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "wolf.latent"), spp = c(".harvest", ".bear", ".wtd", ".wolf"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #108
    # # bs_topdown_final[[108]] # "wolf.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #109         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[109]] # "wolf.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #110
    # bs_topdown_final[[110]] # "wolf.tmin1"     "lion.t"         "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "wolf.latent"), spp = c(".wtd", ".wtd", ".elk", ".harvest", ".wolf"), indices = c(1,2,1,1,1), lags = c("y-1","y","y-1","y-1","y-1")),
    # Regression #111         # NOTE: NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_final[[111]] # "moose.t"     "wtd.tmin1"   "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #112
    # bs_topdown_final[[112]] # "moose.t"     "coy.t"       "moose.tmin1" "wolf.tmin1"  "coy.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("moose.latent", "wolf.latent", "coy.latent", "wtd.latent", "moose.latent"), spp = c(".moose", ".wolf", ".coy", ".wtd", ".moose"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #113
    # bs_topdown_final[[113]] # "moose.t"        "wtd.t"          "moose.tmin1"    "wolf.tmin1"     "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wtd.tmin1"
    list(dSep_test = 7, covariates = c("moose.latent", "wolf.latent", "deerHarv", "lion.latent", "coy.latent", "wtd.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".coy", ".wtd", ".moose"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #114
    # bs_topdown_final[[114]] # "moose.t"        "bear.t"         "moose.tmin1"    "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("moose.latent", "wolf.latent", "bearHarv", "bear.latent", "wtd.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".wtd", ".moose"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y","y")),
    # Regression #115         # NOTE: NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_final[[115]] # "moose.t"     "elk.tmin1"   "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #116
    # bs_topdown_final[[116]] # "moose.t"       "elk.t"         "moose.tmin1"   "wolf.tmin1"    "elkHarv.tmin1" "lion.tmin1"    "bear.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elkHarv", "lion.latent", "bear.latent", "elk.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".bear", ".elk", ".moose"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #117
    # bs_topdown_final[[117]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1" "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "wtd.latent", "wtd.latent", "elk.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".wtd", ".wtd", ".elk", ".moose"), indices = c(1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y")),
    # Regression #118         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[118]] # "moose.t"        "lionHarv.tmin1" "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #119
    # bs_topdown_final[[119]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wolf", ".wtd", ".wtd", ".elk", ".harvest", ".moose"), indices = c(1,1,1,2,1,1,2), lags = c("y-1","y-1","y-1","y","y-1","y-1","y")),
    # Regression #120
    # bs_topdown_final[[120]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wtd", ".wtd"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #121
    # # bs_topdown_final[[121]] # "wtd.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #122
    # bs_topdown_final[[122]] # "wtd.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "bear.tmin1"    "wolf.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "wtd.latent"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".wtd"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #123         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[123]] # "wtd.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #124
    # bs_topdown_final[[124]] # "coy.t"          "wtd.t"          "coy.tmin1"      "wtd.tmin1"      "deerHarv.tmin1" "lion.tmin1"     "wolf.tmin1"
    list(dSep_test = 7, covariates = c("coy.latent", "wtd.latent", "deerHarv", "lion.latent", "wolf.latent", "coy.latent"), spp = c(".coy", ".wtd", ".harvest", ".lion", ".wolf", ".coy"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #125
    # bs_topdown_final[[125]] # "coy.t"          "bear.t"         "coy.tmin1"      "wtd.tmin1"      "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("coy.latent", "wtd.latent", "bearHarv", "bear.latent", "wtd.latent", "coy.latent"), spp = c(".coy", ".wtd", ".harvest", ".bear", ".wtd", ".coy"), indices = c(1,1,1,1,2,2), lags = c("y-1","y-1","y-1","y-1","y","y")),
    # Regression #126         # NOTE: NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_final[[126]] # "coy.t"     "elk.tmin1" "coy.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "elk.latent"), spp = c(".coy", ".wtd", ".elk"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #127
    # bs_topdown_final[[127]] # "coy.t"         "elk.t"         "coy.tmin1"     "wtd.tmin1"     "elkHarv.tmin1" "lion.tmin1"    "bear.tmin1"    "wolf.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("coy.latent", "wtd.latent", "elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "coy.latent"), spp = c(".coy", ".wtd", ".harvest", ".lion", ".bear", ".wolf", ".elk", ".coy"), indices = c(1,1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #128
    # bs_topdown_final[[128]] # "coy.t"          "wolf.t"         "coy.tmin1"      "wtd.tmin1"      "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.t"   "elk.tmin1"
    list(dSep_test = 2, covariates = c("coy.latent", "wtd.latent", "moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "elk.latent", "coy.latent"), spp = c(".coy", ".wtd", ".moose", ".harvest", ".wolf", ".wtd", ".elk", ".coy"), indices = c(1,1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y","y-1","y")),
    # Regression #129         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[129]] # "coy.t"          "lionHarv.tmin1" "coy.tmin1"      "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "lionHarv"), spp = c(".coy", ".wtd", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #130
    # bs_topdown_final[[130]] # "coy.t"          "lion.t"         "coy.tmin1"      "wtd.tmin1"      "wtd.t"          "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".elk", ".harvest", ".coy"), indices = c(1,1,2,1,1,2), lags = c("y-1","y-1","y","y-1","y-1","y")),
    # Regression #131         # NOTE: NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_final[[131]] # "wtd.t"          "elk.tmin1"      "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "elk.latent"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".elk"), indices = c(1,1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #132
    # bs_topdown_final[[132]] # "wtd.t"          "elk.t"          "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"      "elkHarv.tmin1"     "bear.tmin1"     "elk.tmin1"
    list(dSep_test = 5, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "elkHarv", "bear.latent", "elk.latent", "wtd.latent"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".harvest", ".bear", ".elk", ".wtd"), indices = c(1,1,1,1,1,2,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #133        # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[133]] # "wtd.t"          "lionHarv.tmin1" "deerHarv.tmin1" "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wtd.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent", "lionHarv"), spp = c(".harvest", ".lion", ".coy", ".wolf", ".wtd", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #134         # NOTE: NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_final[[134]] # "bear.t"         "elk.tmin1"      "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "elk.latent"), spp = c(".harvest", ".bear", ".wtd", ".elk"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #135
    # bs_topdown_final[[135]] # "bear.t"         "elk.t"          "bearHarv.tmin1" "bear.tmin1"     "wtd.t"          "elkHarv.tmin1"  "lion.tmin1"     "wolf.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("bearHarv", "bear.latent", "wtd.latent", "elkHarv", "lion.latent", "wolf.latent", "elk.latent", "bear.latent"), spp = c(".harvest", ".bear", ".wtd", ".harvest", ".lion", ".wolf", ".elk", ".bear"), indices = c(1,1,1,2,1,1,1,2), lags = c("y-1","y-1","y","y-1","y-1","y-1","y-1","y")),
    # Regression #136
    # bs_topdown_final[[136]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"          "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"    "wtd.tmin1"      "elk.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wtd.latent", "moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "elk.latent", "bear.latent"), spp = c(".harvest", ".bear", ".wtd", ".moose", ".harvest", ".wolf", ".wtd", ".elk", ".bear"), indices = c(1,1,1,1,2,1,2,1,2), lags = c("y-1","y-1","y","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #137         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[137]] # "bear.t"         "lionHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"     "wtd.t"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent", "lionHarv"), spp = c(".harvest", ".bear", ".wtd", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # Regression #138
    # bs_topdown_final[[138]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "wtd.t"          "wtd.tmin1"      "elk.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wtd", ".wtd", ".elk", ".harvest", ".bear"), indices = c(1,1,1,2,1,2,2), lags = c("y-1","y-1","y","y-1","y-1","y-1","y")),
    # #########
    # # Regression #139         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_final[[139]] # "elk.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #140
    # bs_topdown_final[[140]] # "elk.t"          "wolf.t"         "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"      "moose.tmin1"   "wolfHarv.tmin1" "wtd.tmin1"      "wtd.t"
    list(dSep_test = 2, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "moose.latent", "wolfHarv", "wtd.latent", "wtd.latent", "elk.latent"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".moose", ".harvest", ".wtd", ".wtd", ".elk"), indices = c(1,1,1,1,1,1,2,1,2,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y-1","y-1","y","y")),
    # Regression #141         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[141]] # "elk.t"          "lionHarv.tmin1" "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "lionHarv"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".harvest"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #142
    # bs_topdown_final[[142]] # "elk.t"          "lion.t"         "elkHarv.tmin1"  "lion.tmin1"     "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"      "wtd.tmin1"  "wtd.t"          "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elkHarv", "lion.latent", "bear.latent", "wolf.latent", "elk.latent", "wtd.latent", "wtd.latent", "lionHarv", "elk.latent"), spp = c(".harvest", ".lion", ".bear", ".wolf", ".elk", ".wtd", ".wtd", ".harvest", ".elk"), indices = c(1,1,1,1,1,1,2,2,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y","y-1","y")),
    # Regression #143         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # bs_topdown_final[[143]] # "wolf.t"         "lionHarv.tmin1" "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "lionHarv"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".harvest"), indices = c(1,1,1,1,2,1,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1")),
    # Regression #144
    # bs_topdown_final[[144]] # "wolf.t"         "lion.t"         "moose.tmin1"    "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "wtd.t"          "elk.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolfHarv", "wolf.latent", "wtd.latent", "wtd.latent", "elk.latent", "lionHarv", "wolf.latent"), spp = c(".moose", ".harvest", ".wolf", ".wtd", ".wtd", ".elk", ".harvest", ".wolf"), indices = c(1,1,1,1,2,1,2,2), lags = c("y-1","y-1","y-1","y-1","y","y-1","y-1","y"))
    )
