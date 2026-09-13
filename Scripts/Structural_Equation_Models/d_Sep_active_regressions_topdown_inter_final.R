  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down interference FINAL
  #'  Sept 2026
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
  #'  d_Sep_active_regression_topdown_inter_tmin1_only.R  and 
  #'  d_Sep_active_regression_topdown_inter_exog_only.R
  #'  ---------------------------------------------------------------- 
  
  dSep_iterations_topdown_inter_final <- list(
    # #########
    # # Regression #1
    # # bs_topdown_inter_final[[1]] # "wtd.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #2
    # # bs_topdown_inter_final[[2]] # "wtd.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #3
    # # bs_topdown_inter_final[[3]] # "wtd.tmin1"  "lion.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #4
    # # bs_topdown_inter_final[[4]] # "wtd.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #5               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[5]] # "wtd.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #6
    # # bs_topdown_inter_final[[6]] # "wtd.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #7               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[7]] # "wtd.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #8               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[8]] # "wtd.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #9
    # # bs_topdown_inter_final[[9]] # "wtd.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #10
    # bs_topdown_inter_final[[10]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #11
    # bs_topdown_inter_final[[11]] # "wtd.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wtd.latent"), spp = c(".elk", ".lion", ".wolf", ".wtd"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #12
    # bs_topdown_inter_final[[12]] # "wtd.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".lion", ".coy", ".wolf", ".wtd"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #13
    # bs_topdown_inter_final[[13]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"   
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wolf", ".wtd"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #14
    # bs_topdown_inter_final[[14]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"   
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".bear", ".wtd"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #15
    # bs_topdown_inter_final[[15]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #16
    # # bs_topdown_inter_final[[16]] # "moose.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #17
    # # bs_topdown_inter_final[[17]] # "moose.tmin1" "lion.tmin1" 
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #18
    # # bs_topdown_inter_final[[18]] # "moose.tmin1" "coy.tmin1"  
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #19               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[19]] # "moose.tmin1"    "bearHarv.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #20
    # # bs_topdown_inter_final[[20]] # "moose.tmin1" "bear.tmin1" 
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #21
    # bs_topdown_inter_final[[21]] # "moose.tmin1" "wtd.t"       "wtd.tmin1"   "lion.tmin1"  "coy.tmin1"   "bear.tmin1" 
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "moose.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".moose"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #22               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[22]] # "moose.tmin1"    "wolfHarv.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #23               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[23]] # "moose.tmin1"    "lionHarv.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #24
    # # bs_topdown_inter_final[[24]] # "moose.tmin1" "wolf.tmin1" 
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #25
    # bs_topdown_inter_final[[25]] # "moose.tmin1" "elk.t"       "elk.tmin1"   "lion.tmin1"  "wolf.tmin1" 
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "moose.latent"), spp = c(".elk", ".lion", ".wolf", ".moose"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #26
    # bs_topdown_inter_final[[26]] # "moose.tmin1" "coy.t"       "lion.tmin1"  "coy.tmin1"   "wolf.tmin1" 
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "moose.latent"), spp = c(".lion", ".coy", ".wolf", ".moose"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #27
    # bs_topdown_inter_final[[27]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"   
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".bear", ".wolf", ".moose"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #28
    # bs_topdown_inter_final[[28]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"     
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".bear", ".moose"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #29
    # bs_topdown_inter_final[[29]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1" "wolf.tmin1"  
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #30
    # # bs_topdown_inter_final[[30]] # "elk.tmin1"  "lion.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #31
    # # bs_topdown_inter_final[[31]] # "elk.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #32               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[32]] # "elk.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #33
    # # bs_topdown_inter_final[[33]] # "elk.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #34
    # bs_topdown_inter_final[[34]] # "elk.tmin1"  "wtd.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "bear.tmin1"
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "elk.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".elk"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #35               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[35]] # "elk.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #36               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[36]] # "elk.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #37
    # # bs_topdown_inter_final[[37]] # "elk.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #38
    # bs_topdown_inter_final[[38]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #39
    # bs_topdown_inter_final[[39]] # "elk.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "elk.latent"), spp = c(".lion", ".coy", ".wolf", ".elk"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #40
    # bs_topdown_inter_final[[40]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"  
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".harvest", ".bear", ".wolf", ".elk"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #41
    # bs_topdown_inter_final[[41]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"    
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".bear", ".elk"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #42
    # bs_topdown_inter_final[[42]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #43
    # # bs_topdown_inter_final[[43]] # "lion.tmin1" "coy.tmin1" 
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #44               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[44]] # "lion.tmin1"     "bearHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #45
    # # bs_topdown_inter_final[[45]] # "lion.tmin1" "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #46               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[46]] # "lion.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #47               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[47]] # "lion.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = 1, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #48
    # # bs_topdown_inter_final[[48]] # "lion.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #49
    # bs_topdown_inter_final[[49]] # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #50
    # bs_topdown_inter_final[[50]] # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"    
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lion.latent"), spp = c(".harvest", ".bear", ".wolf", ".lion"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #51
    # bs_topdown_inter_final[[51]] # "lion.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".bear", ".lion"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #52
    # bs_topdown_inter_final[[52]] # "lion.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #53               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[53]] # "coy.tmin1"      "bearHarv.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #54
    # # bs_topdown_inter_final[[54]] # "coy.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #55               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[55]] # "coy.tmin1"      "wolfHarv.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #56               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[56]] # "coy.tmin1"      "lionHarv.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #57
    # # bs_topdown_inter_final[[57]] # "coy.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #58
    # bs_topdown_inter_final[[58]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #59
    # bs_topdown_inter_final[[59]] # "coy.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent"), spp = c(".elk", ".lion", ".wolf", ".coy"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #60
    # bs_topdown_inter_final[[60]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"  
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".bear", ".wolf", ".coy"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #61
    # bs_topdown_inter_final[[61]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t" 
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".bear", ".coy"), indices = c(1,1,1,1), lags = c("y-1","y-1","y","y-1")),
    # Regression #62
    # bs_topdown_inter_final[[62]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #63
    # # bs_topdown_inter_final[[63]] # "bearHarv.tmin1" "bear.tmin1"    
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #64
    # bs_topdown_inter_final[[64]] # "bearHarv.tmin1" "wtd.t"          "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1"
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "bearHarv"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #65
    # # bs_topdown_inter_final[[65]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),     # exog on exog
    # #########
    # # Regression #66
    # # bs_topdown_inter_final[[66]] # "bearHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),     # exog on exog
    # #########
    # # Regression #67
    # # bs_topdown_inter_final[[67]] # "bearHarv.tmin1" "wolf.tmin1"    
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #68
    # bs_topdown_inter_final[[68]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"    
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #69
    # bs_topdown_inter_final[[69]] # "bearHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"   
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #70
    # bs_topdown_inter_final[[70]] # "bearHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1" 
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #71
    # bs_topdown_inter_final[[71]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"        
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".bear", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # Regression #72
    # bs_topdown_inter_final[[72]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = c(1,1,2), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #73               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[73]] # "bear.tmin1"     "wolfHarv.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #74               # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
    # # bs_topdown_inter_final[[74]] # "bear.tmin1"     "lionHarv.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = c(), lags = c()),
    # #########
    # # Regression #75
    # # bs_topdown_inter_final[[75]] # "bear.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #76
    # bs_topdown_inter_final[[76]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf",  ".bear"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #77
    # bs_topdown_inter_final[[77]] # "bear.tmin1" "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bear.latent"), spp = c(".elk", ".lion", ".wolf", ".bear"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #78
    # bs_topdown_inter_final[[78]] # "bear.tmin1" "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bear.latent"), spp = c(".lion", ".coy", ".wolf", ".bear"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #79
    # bs_topdown_inter_final[[79]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"     
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear", ".bear"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # Regression #80
    # bs_topdown_inter_final[[80]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #81               #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_inter_final[[81]] # "wtd.t"          "wolfHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1" 
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "wolfHarv"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #82               #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_inter_final[[82]] # "wtd.t"          "lionHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1"   
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "lionHarv"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #83               #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown_inter_final[[83]] # "wtd.t"      "wolf.tmin1" "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "bear.tmin1"
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "wolf.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".wolf"), indices = c(1,1,1,1,1), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #84
    # bs_topdown_inter_final[[84]] # "wtd.t"       "moose.t"     "wtd.tmin1"   "lion.tmin1"  "coy.tmin1"   "bear.tmin1"  "moose.tmin1" "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".moose", ".wolf", ".wtd"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #85
    # bs_topdown_inter_final[[85]] # "wtd.t"      "elk.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "bear.tmin1" "elk.tmin1"  "wolf.tmin1"
    list(dSep_test = 5, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".elk", ".wolf", ".wtd"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #86
    # bs_topdown_inter_final[[86]] # "wtd.t"      "coy.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".wolf", ".wtd"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #87
    # bs_topdown_inter_final[[87]] # "wtd.t"          "bear.t"         "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1"     "bearHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "bearHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest", ".wolf", ".wtd"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #88
    # bs_topdown_inter_final[[88]] # "wtd.t"          "wolf.t"         "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"
    list(dSep_test = 2, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "wolfHarv", "wolf.latent", "bear.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest", ".wolf", ".bear", ".wtd"), indices = c(1,1,1,1,1,1,2,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y","y")),
    # Regression #89
    # bs_topdown_inter_final[[89]] # "wtd.t"          "lion.t"         "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"      "bear.tmin1"     "lionHarv.tmin1" "wolf.tmin1"  
    list(dSep_test = 1, covariates = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent", "lionHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".bear", ".harvest", ".wolf", ".wtd"), indices = c(1,1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # #########
    # # Regression #90
    # # bs_topdown_inter_final[[90]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(), lags = c()),     # exog on exog
    # #########
    # # Regression #91
    # # bs_topdown_inter_final[[91]] # "wolfHarv.tmin1" "wolf.tmin1"    
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #92
    # bs_topdown_inter_final[[92]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"  
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #93
    # bs_topdown_inter_final[[93]] # "wolfHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"  
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #94
    # bs_topdown_inter_final[[94]] # "wolfHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"   
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #95
    # bs_topdown_inter_final[[95]] # "wolfHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"   
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #96
    # bs_topdown_inter_final[[96]] # "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"   
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = c(1,1,2), lags = c("y-1","y-1","y-1")),
    # #########
    # # Regression #97
    # # bs_topdown_inter_final[[97]] # "lionHarv.tmin1" "wolf.tmin1"    
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = c(), lags = c()),
    # Regression #98
    # bs_topdown_inter_final[[98]] # "lionHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"    
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = c(1,1,1), lags = c("y-1","y-1","y-1")),
    # Regression #99
    # bs_topdown_inter_final[[99]] # "lionHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"  
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #100
    # bs_topdown_inter_final[[100]] # "lionHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"    
    list(dSep_test = 3, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = c(1,1,1,1), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #101
    # bs_topdown_inter_final[[101]] # "lionHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"   
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #102
    # bs_topdown_inter_final[[102]] # "lionHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"     
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "lionHarv"), spp = c(".harvest", ".wolf", ".bear", ".harvest"), indices = c(1,1,1,2), lags = c("y-1","y-1","y","y-1")),
    # Regression #103
    # bs_topdown_inter_final[[103]] # "moose.t"     "elk.t"       "moose.tmin1" "wolf.tmin1"  "elk.tmin1"   "lion.tmin1" 
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elk.latent", "lion.latent", "moose.latent"), spp = c(".moose", ".wolf", ".elk", ".lion", ".moose"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #104
    # bs_topdown_inter_final[[104]] # "moose.t"     "coy.t"       "moose.tmin1" "wolf.tmin1"  "lion.tmin1"  "coy.tmin1"  
    list(dSep_test = 4, covariates = c("moose.latent", "wolf.latent", "lion.latent", "coy.latent", "moose.latent"), spp = c(".moose", ".wolf", ".lion", ".coy", ".moose"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #105
    # bs_topdown_inter_final[[105]] # "moose.t"        "bear.t"         "moose.tmin1"    "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"   
    list(dSep_test = 3, covariates = c("moose.latent", "wolf.latent", "bearHarv", "bear.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #106
    # bs_topdown_inter_final[[106]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1" "bear.t"       
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "bear.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y","y")),
    # Regression #107
    # bs_topdown_inter_final[[107]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = c(1,1,1,2), lags = c("y-1","y-1","y-1","y")),
    # Regression #108
    # bs_topdown_inter_final[[108]] # "elk.t"      "coy.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1" "coy.tmin1" 
    list(dSep_test = 4, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".coy", ".elk"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #109
    # bs_topdown_inter_final[[109]] # "elk.t"          "bear.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"    
    list(dSep_test = 3, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv", "bear.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".bear", ".elk"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #110
    # bs_topdown_inter_final[[110]] # "elk.t"          "wolf.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1" "bear.t"    
    list(dSep_test = 2, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv", "bear.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".bear", ".elk"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y","y")),
    # Regression #111
    # bs_topdown_inter_final[[111]] # "elk.t"          "lion.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #112
    # bs_topdown_inter_final[[112]] # "coy.t"          "bear.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"    
    list(dSep_test = 3, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".bear", ".coy"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #113
    # bs_topdown_inter_final[[113]] # "coy.t"          "wolf.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wolfHarv.tmin1" "bear.t"     
    list(dSep_test = 2, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv", "bear.latent", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".bear", ".coy"), indices = c(1,1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y","y")),
    # Regression #114
    # bs_topdown_inter_final[[114]] # "coy.t"          "lion.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = c(1,1,1,1,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #115
    # bs_topdown_inter_final[[115]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = c(1,1,1,2,2), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #116
    # bs_topdown_inter_final[[116]] # "wolf.t"         "lion.t"         "wolfHarv.tmin1" "wolf.tmin1"     "bear.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("wolfHarv", "wolf.latent", "bear.latent", "lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf", ".bear", ".harvest", ".wolf"), indices = c(1,1,1,2,2), lags = c("y-1","y-1","y","y-1","y"))
  )
  
  

  
  
  
  
  
  
  
  
  