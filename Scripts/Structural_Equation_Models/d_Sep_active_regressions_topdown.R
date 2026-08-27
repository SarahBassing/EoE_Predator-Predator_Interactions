  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down
  #'  Aug 2026
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
  
  dSep_iterations_topdown <- list(
    
    # # Regression #1
    # # bs_topdown[[1]] # "deerHarv.tmin1" "wtd.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #2
    # # bs_topdown[[2]] # "deerHarv.tmin1" "moose.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #3
    # # bs_topdown[[3]] # "deerHarv.tmin1" "elkHarv.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #4
    # # bs_topdown[[4]] # "deerHarv.tmin1" "lion.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #5
    # # bs_topdown[[5]] # "deerHarv.tmin1" "elk.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #6
    # # bs_topdown[[6]] # "deerHarv.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #7
    # bs_topdown[[7]] "deerHarv.tmin1" "coy.t"          "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "deerHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #8
    # # bs_topdown[[8]] # "deerHarv.tmin1" "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #9
    # # bs_topdown[[9]] # "deerHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #10
    # bs_topdown[[10]] # "deerHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "deerHarv"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # # Regression #11
    # # bs_topdown[[11]] # "deerHarv.tmin1" "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #12
    # # bs_topdown[[12]] # "deerHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #13
    # bs_topdown[[13]] #  "deerHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "deerHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #14
    # bs_topdown[[14]] #  "deerHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "deerHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #15
    # bs_topdown[[15]] # "deerHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "deerHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # # Regression #16
    # # bs_topdown[[16]] # "deerHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("deerHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # Regression #17
    # bs_topdown[[17]] # "deerHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "deerHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # # Regression #18
    # # bs_topdown[[18]] # "wtd.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # #########
    # # Regression #19   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[19]] # "wtd.tmin1"     "elkHarv.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # #########
    # # Regression #20
    # # bs_topdown[[20]] # "wtd.tmin1"  "lion.tmin1"
    # list(dSep_test = NA, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # # Regression #21
    # # bs_topdown[[21]] # "wtd.tmin1" "elk.tmin1"
    # list(dSep_test = NA, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # # Regression #22
    # # bs_topdown[[22]] # "wtd.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # Regression #23
    # bs_topdown[[23]] # "wtd.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # # Regression #24   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[24]] # "wtd.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # # Regression #25
    # # bs_topdown[[25]] # "wtd.tmin1"  "bear.tmin1"
    # list(dSep_test = NA, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # Regression #26
    # bs_topdown[[26]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #27   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[27]] # "wtd.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # # Regression #28
    # # bs_topdown[[28]] # "wtd.tmin1"  "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
    # Regression #29
    # bs_topdown[[29]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #30
    # bs_topdown[[30]] # "wtd.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #31
    # bs_topdown[[31]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #32   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[32]] # "wtd.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # Regression #33
    # bs_topdown[[33]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wtd.latent"), spp = c(".harvest", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # # Regression #34   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[34]] # "moose.tmin1"   "elkHarv.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # # Regression #35
    # # bs_topdown[[35]] # "moose.tmin1" "lion.tmin1"
    # list(dSep_test = NA, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")), 
    # Regression #36
    # bs_topdown[[36]] # "moose.tmin1"    "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #37
    # # bs_topdown[[37]] # "moose.tmin1" "elk.tmin1"  
    # list(dSep_test = NA, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")), 
    # # Regression #38
    # # bs_topdown[[38]] # "moose.tmin1" "coy.tmin1" 
    # list(dSep_test = NA, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")), 
    # Regression #39
    # bs_topdown[[39]] #  "moose.tmin1" "coy.t"       "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "moose.latent"), spp = c(".coy", ".moose"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # # Regression #40   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[40]] # "moose.tmin1"    "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
    # ##########
    # # Regression #41
    # # bs_topdown[[41]] # "moose.tmin1" "bear.tmin1" 
    # list(dSep_test = NA, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")), 
    # Regression #42
    # bs_topdown[[42]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "moose.latent"), spp = c(".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #43   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[43]] # "moose.tmin1"    "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # # Regression #44
    # # bs_topdown[[44]] # "moose.tmin1" "wolf.tmin1" 
    # list(dSep_test = NA, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
    # # Regression #45
    # bs_topdown[[45]] # "moose.tmin1"   "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #46
    # bs_topdown[[46]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #47   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[47]] # "moose.tmin1"    "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########
    # Regression #48
    # bs_topdown[[48]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "moose.latent"), spp = c(".harvest",  ".moose"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #49
    # # bs_topdown[[49]] # "elkHarv.tmin1" "lion.tmin1" 
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #50
    # bs_topdown[[50]] #  "elkHarv.tmin1"  "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #51
    # # bs_topdown[[51]] # "elkHarv.tmin1" "elk.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # # Regression #52
    # # bs_topdown[[52]] # "elkHarv.tmin1" "coy.tmin1" 
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #53
    # bs_topdown[[53]] # "elkHarv.tmin1" "coy.t"         "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "elkHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #54
    # # bs_topdown[[54]] # "elkHarv.tmin1"  "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #55 
    # # bs_topdown[[55]]  # "elkHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #56
    # bs_topdown[[56]]  # "elkHarv.tmin1"  "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "elkHarv"), spp = c(".harvest", ".bear", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # # Regression #57 
    # # bs_topdown[[57]]  # "elkHarv.tmin1"  "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #58 
    # # bs_topdown[[58]]  # "elkHarv.tmin1" "wolf.tmin1"  
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # Regression #59
    # bs_topdown[[59]]  # "elkHarv.tmin1" "moose.t"       "moose.tmin1"   "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elkHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #60
    # bs_topdown[[60]]  # "elkHarv.tmin1"  "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elkHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # # Regression #61
    # # bs_topdown[[61]]  # "elkHarv.tmin1"  "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("elkHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
    # Regression #62
    # bs_topdown[[62]]  # "elkHarv.tmin1"  "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv",  "elkHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # # Regression #63
    # # bs_topdown[[63]]  # "lion.tmin1" "elk.tmin1"
    # list(dSep_test = NA, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
    # # Regression #64
    # # bs_topdown[[64]]  # "lion.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
    # Regression #65
    # bs_topdown[[65]]  # "lion.tmin1" "coy.t"      "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "lion.latent"), spp = c(".coy", ".lion"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########    
    # # Regression #66   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[66]] # "lion.tmin1"     "bearHarv.tmin1"
    # list(dSep_test = 1, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########    
    # # Regression #67
    # # bs_topdown[[67]]  # "lion.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
    # Regression #68
    # bs_topdown[[68]]  # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "lion.latent"), spp = c(".harvest", ".bear", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########    
    # # Regression #69   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[69]] # "lion.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
    # ##########    
    # # Regression #70
    # # bs_topdown[[70]]  # "lion.tmin1" "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
    # Regression #71
    # bs_topdown[[71]]  # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #72
    # bs_topdown[[72]]  # "lion.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########    
    # # Regression #73    # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[73]]  # "lion.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
    # ##########    
    # Regression #74
    # bs_topdown[[74]]  # "lion.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "lion.latent"), spp = c(".harvest", ".lion"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #75    # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[75]]  # "wtd.t"          "elk.tmin1"      "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elk.latent"), spp = c(".harvest", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #76   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[76]] # "wtd.t"          "coy.tmin1"      "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "coy.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # ##########
    # Regression #77
    # bs_topdown[[77]]  # "wtd.t"          "coy.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("deerHarv", "wtd.latent", "lion.latent", "coy.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #78   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[78]] # "wtd.t"          "bearHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bearHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #79   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[79]] # "wtd.t"          "bear.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bear.latent"), spp = c(".harvest", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # ##########
    # Regression #80
    # bs_topdown[[80]] # "wtd.t"          "bear.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bearHarv", "bear.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".bear", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #81   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[81]] # "wtd.t"          "wolfHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolfHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #82   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[82]] # "wtd.t"          "wolf.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolf.latent"), spp = c(".harvest", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # ##########
    # Regression #83
    # bs_topdown[[83]] # "wtd.t"          "moose.t"        "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #84
    # bs_topdown[[84]] # "wtd.t"          "elk.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #85
    # bs_topdown[[85]] # "wtd.t"          "wolf.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #86   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[86]] #  "wtd.t"          "lionHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "lionHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # ##########
    # Regression #87
    # bs_topdown[[87]] # "wtd.t"          "lion.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("deerHarv", "wtd.latent", "lion.latent", "lionHarv", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wtd"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # # Regression #88
    # # bs_topdown[[88]] # "elk.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c("elk.latent"), spp = c(".elk"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #89
    # bs_topdown[[89]] # "elk.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "elk.latent"), spp = c(".coy", ".elk"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # # Regression #90   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[90]] # "elk.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # # Regression #91
    # # bs_topdown[[91]] # "elk.tmin1"  "bear.tmin1"
    # list(dSep_test = NA, covariates = c("elk.latent"), spp = c(".elk"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #92
    # bs_topdown[[92]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "elk.latent"), spp = c(".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #93   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[93]] # "elk.tmin1"  "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # # Regression #94
    # # bs_topdown[[94]] # "elk.tmin1"  "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("elk.latent"), spp = c(".elk"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #95
    # bs_topdown[[95]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #96
    # bs_topdown[[96]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #97   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[97]] # "elk.tmin1"  "lionHarv.tmin1"
    # list(dSep_test = 5, covariates = c("lionHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # Regression #98
    # bs_topdown[[98]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "elk.latent"), spp = c(".harvest", ".elk"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # # Regression #99   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[99]] # "coy.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # # Regression #100
    # # bs_topdown[[100]] # "coy.tmin1"      "bear.tmin1"
    # list(dSep_test = NA, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #101
    # bs_topdown[[101]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "coy.latent"), spp = c(".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #102   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[102]] # "coy.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # # Regression #103
    # # bs_topdown[[103]] # "coy.tmin1"      "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #104
    # bs_topdown[[104]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #105
    # bs_topdown[[105]] # "coy.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #106
    # bs_topdown[[106]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # # Regression #107   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[107]] # "coy.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ##########
    # Regression #108
    # bs_topdown[[108]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "coy.latent"), spp = c(".harvest", ".coy"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #109   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[109]] # "coy.t"          "bearHarv.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "bearHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #110   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[110]] # "coy.t"      "bear.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "bear.latent"), spp = c(".coy", ".bear"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #111
    # bs_topdown[[111]] # "coy.t"          "bear.t"         "coy.tmin1"      "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("coy.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".coy", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # ##########
    # Regression #112   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[112]] # "coy.t"          "wolfHarv.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolfHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #113   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[113]] # "coy.t"      "wolf.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent"), spp = c(".coy", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #114
    # bs_topdown[[114]] # "coy.t"       "moose.t"     "coy.tmin1"   "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("coy.latent", "moose.latent", "wolf.latent", "coy.latent"), spp = c(".coy", ".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #115
    # bs_topdown[[115]] #  "coy.t"         "elk.t"         "coy.tmin1"     "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("coy.latent", "elkHarv", "lion.latent", "elk.latent", "wolf.latent", "coy.latent"), spp = c(".coy", ".harvest", ".lion", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #116
    # bs_topdown[[116]] # "coy.t"          "wolf.t"         "coy.tmin1"      "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("coy.latent", "wolfHarv", "wolf.latent", "coy.latent"), spp = c(".coy", ".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # ##########
    # Regression #117   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[117]] # "coy.t"          "lionHarv.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "lionHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #118
    # bs_topdown[[118]] # "coy.t"          "lion.t"         "coy.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("coy.latent", "lionHarv", "coy.latent"), spp = c(".coy", ".harvest", ".coy"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y")),
    # # Regression #119
    # # bs_topdown[[119]] # "bearHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # # Regression #120
    # # # bs_topdown[[120]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    # # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),   # harvest effects harvest...
    # # Regression #121
    # # bs_topdown[[121]] # "bearHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #122
    # bs_topdown[[122]] #  "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #123
    # bs_topdown[[123]] #  "bearHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bearHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #124
    # bs_topdown[[124]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # # Regression #125
    # # bs_topdown[[125]] # "bearHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("bearHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),   # harvest effects harvest...
    # Regression #126
    # bs_topdown[[126]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "bearHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # ############
    # # Regression #127   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[127]] # "bear.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = 3, covariates = c("wolfHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ############
    # # Regression #128
    # # bs_topdown[[128]] # "bear.tmin1" "wolf.tmin1"
    # list(dSep_test = NA, covariates = c("bear.latent"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #129
    # bs_topdown[[129]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #130
    # bs_topdown[[130]] # "bear.tmin1"    "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #131
    # bs_topdown[[131]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ############
    # # Regression #132   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[132]] # "bear.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # ############
    # Regression #133
    # bs_topdown[[133]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "bear.latent"), spp = c(".harvest", ".bear"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ##########
    # Regression #134   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[134]] # "bear.t"         "wolfHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #135   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[135]] # "bear.t"         "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent"), spp = c(".harvest", ".bear", ".wolf"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # Regression #136
    # bs_topdown[[136]] # "bear.t"         "moose.t"        "bearHarv.tmin1" "bear.tmin1"     "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("bearHarv", "bear.latent", "moose.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #137
    # bs_topdown[[137]] # "bear.t"         "elk.t"          "bearHarv.tmin1" "bear.tmin1"     "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("bearHarv", "bear.latent", "elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".lion", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,2,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #138
    # bs_topdown[[138]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #139   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[139]] # "bear.t"         "lionHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "lionHarv"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # ##########
    # Regression #140
    # bs_topdown[[140]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # # Regression #141
    # # bs_topdown[[141]] # "wolfHarv.tmin1" "wolf.tmin1" 
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #142
    # bs_topdown[[142]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #143
    # bs_topdown[[143]] # "wolfHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # # Regression #144
    # # bs_topdown[[144]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("wolfHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),   # harvest effects harvest...
    # Regression #145
    # bs_topdown[[145]] #  "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolfHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # #########
    # # Regression #146   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown[[146]] # "wolf.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c("lionHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    # #########
    # Regression #147
    # bs_topdown[[147]] #  "wolf.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #148
    # bs_topdown[[148]] # "moose.t"       "elk.t"         "moose.tmin1"   "wolf.tmin1"    "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elkHarv", "lion.latent", "elk.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #149
    # bs_topdown[[149]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # ##########
    # Regression #150   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[150]] # "moose.t"        "lionHarv.tmin1" "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # ##########
    # Regression #151
    # bs_topdown[[151]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wtd", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #152
    # bs_topdown[[152]] # "elk.t"          "wolf.t"         "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wolfHarv", "elk.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #153   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[153]] # "elk.t"          "lionHarv.tmin1" "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # ##########
    # Regression #154
    # bs_topdown[[154]] # "elk.t"          "lion.t"         "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ##########
    # Regression #155 #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[155]] # "wolf.t"         "lionHarv.tmin1" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lionHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # ##########
    # Regression #156
    # bs_topdown[[156]] # "wolf.t"         "lion.t"         "wolfHarv.tmin1" "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wolfHarv", "wolf.latent", "lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf", ".harvest", ".wolf"), indices = as.integer(c(1,1,2,2)), lags = c("y-1","y-1","y-1","y"))
  )
  