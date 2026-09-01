  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: bottom-up interference
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
  #'  d_Sep_active_regression_bottomup_inter_tmin1_only.R  and 
  #'  d_Sep_active_regression_bottomup_inter_exog_only.R
  #'  ---------------------------------------------------------------- 
  
  dSep_iterations_bottomup_inter <- list(
    # # Regression #1
    # # bs_bottomup_inter[[1]] # "wsi.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #2
    # # bs_bottomup_inter[[2]] # "wsi.tmin1"    "forest.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c(1)), lags = c("y-1")),    # enviro effects enviro...
    # # Regression #3
    # # bs_bottomup_inter[[3]] # "wsi.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #4
    # # bs_bottomup_inter[[4]] # "wsi.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #5
    # # bs_bottomup_inter[[5]] # "wsi.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #6
    # # bs_bottomup_inter[[6]] # "wsi.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #7
    # bs_bottomup_inter[[7]] # "wsi.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wolf.latent", "wsi"), spp = c(".forest", ".bear", ".wolf", ".wsi"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #8
    # bs_bottomup_inter[[8]] # "wsi.tmin1"   "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "elk.latent", "wolf.latent", "wsi"), spp = c(".moose", ".elk", ".wolf", ".wsi"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #9
    # # bs_bottomup_inter[[9]] # "wsi.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c("wsi"), spp = c(".wsi"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #10
    # bs_bottomup_inter[[10]] # "wsi.tmin1"  "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent", "wtd.latent", "wsi"), spp = c(".coy", ".wolf", ".wtd", ".wsi"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #4
    # bs_bottomup_inter[[4]] # "wsi.tmin1"  "lion.t"     "wtd.tmin1"  
    list(dSep_test = 1, covariates = c("wtd.latent", "wsi"), spp = c(".wtd", ".wsi"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # ########
    # # Regression #12          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT FOREST
    # # bs_bottomup_inter[[12]] # "coy.tmin1"    "forest.tmin1"
    # list(dSep_test = 4, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #13
    # # bs_bottomup_inter[[13]] # "coy.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #14
    # # bs_bottomup_inter[[14]] # "coy.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #15
    # bs_bottomup_inter[[15]] # "coy.tmin1"    "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "coy.latent"), spp = c(".wsi", ".forest", ".moose", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #16
    # # bs_bottomup_inter[[16]] # "coy.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #17
    # bs_bottomup_inter[[17]] # "coy.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "elk.latent", "coy.latent"), spp = c(".wsi", ".forest", ".elk", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #18
    # # bs_bottomup_inter[[18]] # "coy.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #19
    # bs_bottomup_inter[[19]] # "coy.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wolf.latent", "coy.latent"), spp = c(".forest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #20
    # bs_bottomup_inter[[20]] # "coy.tmin1"   "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "elk.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #21
    # # bs_bottomup_inter[[21]] # "coy.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #22
    # bs_bottomup_inter[[22]] # "coy.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "wtd.latent", "coy.latent"), spp = c(".wsi", ".forest", ".wtd", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #23
    # bs_bottomup_inter[[23]] # "coy.tmin1"  "lion.t"     "wtd.tmin1"  
    list(dSep_test = 1, covariates = c("wtd.latent", "coy.latent"), spp = c(".wtd", ".coy"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #24
    # # bs_bottomup_inter[[24]] # "forest.tmin1" "bear.tmin1"  
    # list(dSep_test = 3, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #25
    # # bs_bottomup_inter[[25]] # "forest.tmin1" "moose.tmin1" 
    # list(dSep_test = 7, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #26
    # # bs_bottomup_inter[[26]] # "forest.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # # Regression #27
    # # bs_bottomup_inter[[27]] # "forest.tmin1" "wolf.tmin1" 
    # list(dSep_test = 2, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #28
    # bs_bottomup_inter[[28]] # "forest.tmin1" "wolf.t"       "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "elk.latent", "wolf.latent", "forest"), spp = c(".moose", ".elk", ".wolf", ".forest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #29
    # # bs_bottomup_inter[[29]] # "forest.tmin1" "wtd.tmin1"   
    # list(dSep_test = 7, covariates = c("forest"), spp = c(".forest"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #30
    # bs_bottomup_inter[[30]] # "forest.tmin1" "coy.t"        "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent", "wtd.latent", "forest"), spp = c(".coy", ".wolf", ".wtd", ".forest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #31
    # bs_bottomup_inter[[31]] # "forest.tmin1" "lion.t"       "wtd.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "forest"), spp = c(".wtd", ".forest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #32
    # # bs_bottomup_inter[[32]] # "bear.tmin1"  "moose.tmin1"
    # list(dSep_test = 6, covariates = c("bear.latent"), spp = c(".bear"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #33
    # bs_bottomup_inter[[33]] # "bear.tmin1"   "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "bear.latent"), spp = c(".wsi", ".forest", ".moose", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #34
    # # bs_bottomup_inter[[34]] # "bear.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c("bear.latent"), spp = c(".bear"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #35
    # bs_bottomup_inter[[35]] # "bear.tmin1"   "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "elk.latent", "bear.latent"), spp = c(".wsi", ".forest", ".elk", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #36
    # # bs_bottomup_inter[[36]] # "bear.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c("bear.latent"), spp = c(".bear"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #37
    # bs_bottomup_inter[[37]] # "bear.tmin1"  "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "elk.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #38
    # # bs_bottomup_inter[[38]] # "bear.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c("bear.latent"), spp = c(".bear"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #39
    # bs_bottomup_inter[[39]] # "bear.tmin1"   "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "wtd.latent", "bear.latent"), spp = c(".wsi", ".forest", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #40
    # bs_bottomup_inter[[40]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent", "wtd.latent", "bear.latent"), spp = c(".coy", ".wolf", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #41
    # bs_bottomup_inter[[41]] # "bear.tmin1" "lion.t"     "wtd.tmin1" 
    list(dSep_test = 1, covariates = c("wtd.latent", "bear.latent"), spp = c(".wtd", ".bear"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # # Regression #42
    # # bs_bottomup_inter[[42]] # "moose.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c("moose.latent"), spp = c(".moose"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #43
    # bs_bottomup_inter[[43]] # "moose.tmin1"  "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".elk", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #44
    # # bs_bottomup_inter[[44]] # "moose.tmin1" "wolf.tmin1" 
    # list(dSep_test = 2, covariates = c("moose.latent"), spp = c(".moose"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #45
    # bs_bottomup_inter[[45]] # "moose.tmin1"  "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".forest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #46
    # # bs_bottomup_inter[[46]] # "moose.tmin1" "wtd.tmin1"  
    # list(dSep_test = 7, covariates = c("moose.latent"), spp = c(".moose"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #47
    # bs_bottomup_inter[[47]] # "moose.tmin1"  "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #48
    # bs_bottomup_inter[[48]] # "moose.tmin1" "coy.t"       "coy.tmin1"   "wolf.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent",  "wtd.latent", "moose.latent"), spp = c(".coy", ".wolf", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #49
    # bs_bottomup_inter[[49]] # "moose.tmin1" "lion.t"      "wtd.tmin1"  
    list(dSep_test = 1, covariates = c("wtd.latent", "moose.latent"), spp = c(".wtd", ".moose"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ########
    # Regression #50          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[50]] # "moose.t"      "elk.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "elk.latent"), spp = c(".wsi", ".forest", ".moose",".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #51
    # bs_bottomup_inter[[51]] # "moose.t"      "elk.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "moose.latent", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # ########
    # Regression #52          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[52]] # "moose.t"      "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "wolf.latent"), spp = c(".wsi", ".forest", ".moose", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #53
    # bs_bottomup_inter[[53]] # "moose.t"      "bear.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wsi", "forest", "moose.latent", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #54
    # bs_bottomup_inter[[54]] # "moose.t"      "wolf.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wsi", "forest", "moose.latent", "elk.latent", "wolf.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".elk", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ########
    # Regression #55          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[55]] # "moose.t"      "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "wtd.latent"), spp = c(".wsi", ".forest", ".moose", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #56
    # bs_bottomup_inter[[56]] # "moose.t"      "wtd.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "moose.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #57
    # bs_bottomup_inter[[57]] # "moose.t"      "coy.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("wsi", "forest", "moose.latent", "coy.latent", "wolf.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".coy", ".wolf", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #58
    # bs_bottomup_inter[[58]] # "moose.t"      "lion.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"   
    list(dSep_test = 1, covariates = c("wsi", "forest", "moose.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # # Regression #59
    # # bs_bottomup_inter[[59]] # "elk.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c("elk.latent"), spp = c(".elk"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #60
    # bs_bottomup_inter[[60]] # "elk.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".forest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # # Regression #61
    # # bs_bottomup_inter[[61]] # "elk.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c("elk.latent"), spp = c(".elk"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #62
    # bs_bottomup_inter[[62]] # "elk.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "wtd.latent", "elk.latent"), spp = c(".wsi", ".forest", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #63
    # bs_bottomup_inter[[63]] # "elk.tmin1"  "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wolf.latent", "wtd.latent", "elk.latent"), spp = c(".coy", ".wolf", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #64
    # bs_bottomup_inter[[64]] # "elk.tmin1" "lion.t"    "wtd.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "elk.latent"), spp = c(".wtd", ".elk"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # ########
    # Regression #65          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[65]] # "elk.t"        "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "elk.latent", "wolf.latent"), spp = c(".wsi", ".forest", ".elk", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #66
    # bs_bottomup_inter[[66]] # "elk.t"        "bear.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wsi", "forest", "elk.latent", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".wsi", ".forest", ".elk", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #67
    # bs_bottomup_inter[[67]] # "elk.t"        "wolf.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "moose.tmin1"  "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wsi", "forest", "elk.latent", "moose.latent", "wolf.latent", "elk.latent"), spp = c(".wsi", ".forest", ".elk", ".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ########
    # Regression #68          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[68]] # "elk.t"        "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "elk.latent", "wtd.latent"), spp = c(".wsi", ".forest", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #69
    # bs_bottomup_inter[[69]] # "elk.t"        "wtd.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi", "forest", "elk.latent", "wtd.latent", "elk.latent"), spp = c(".wsi", ".forest", ".elk", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #70
    # bs_bottomup_inter[[70]] # "elk.t"        "coy.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("wsi", "forest", "elk.latent", "coy.latent", "wolf.latent", "wtd.latent", "elk.latent"), spp = c(".wsi", ".forest", ".elk", ".coy", ".wolf", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #71
    # bs_bottomup_inter[[45]] # "elk.t"        "lion.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "wtd.tmin1"   
    list(dSep_test = 1, covariates = c("wsi", "forest", "elk.latent", "wtd.latent", "elk.latent"), spp = c(".wsi", ".forest", ".elk", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # # Regression #72
    # # bs_bottomup_inter[[72]] # "wolf.tmin1" "wtd.tmin1" 
    # list(dSep_test = 7, covariates = c("wolf.latent"), spp = c(".wolf"), indices = as.integer(c(1)), lags = c("y-1")),
    # Regression #73
    # bs_bottomup_inter[[73]] # "wolf.tmin1"   "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"  
    list(dSep_test = 7, covariates = c("wsi", "forest", "wtd.latent", "wolf.latent"), spp = c(".wsi", ".forest", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #74
    # bs_bottomup_inter[[74]] # "wolf.tmin1" "lion.t"     "wtd.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "wolf.latent"), spp = c(".wtd", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #75
    # bs_bottomup_inter[[75]] # "bear.t"       "wolf.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "moose.tmin1"  "elk.tmin1" 
    list(dSep_test = 2, covariates = c("forest", "bear.latent", "wolf.latent", "moose.latent", "elk.latent", "bear.latent"), spp = c(".forest", ".bear", ".wolf", ".moose", ".elk", ".bear"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # ########
    # Regression #76          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[76]] # "bear.t"       "wtd.tmin1"    "forest.tmin1" "bear.tmin1"   "wolf.tmin1"  
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".forest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #77
    # bs_bottomup_inter[[77]] # "bear.t"       "wtd.t"        "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "wsi.tmin1"    "wtd.tmin1"
    list(dSep_test = 7, covariates = c("forest", "bear.latent", "wolf.latent", "wsi", "wtd.latent", "bear.latent"), spp = c(".forest", ".bear", ".wolf", ".wsi", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #78
    # bs_bottomup_inter[[78]] # "bear.t"       "coy.t"        "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "coy.tmin1"    "wtd.tmin1"
    list(dSep_test = 4, covariates = c("forest", "bear.latent", "wolf.latent", "coy.latent", "wtd.latent", "bear.latent"), spp = c(".forest", ".bear", ".wolf", ".coy", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #79
    # bs_bottomup_inter[[79]] # "bear.t"       "lion.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "wtd.tmin1"    
    list(dSep_test = 1, covariates = c("forest", "bear.latent", "wolf.latent", "wtd.latent", "bear.latent"), spp = c(".forest", ".bear", ".wolf", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #80          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_bottomup_inter[[54]] # "wolf.t"      "wtd.tmin1"   "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #81
    # bs_bottomup_inter[[81]] # "wolf.t"       "wtd.t"        "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("moose.latent", "elk.latent", "wolf.latent", "wsi", "forest", "wtd.latent", "wolf.latent"), spp = c(".moose", ".elk", ".wolf", ".wsi", ".forest", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #82
    # bs_bottomup_inter[[82]] # "wolf.t"      "coy.t"       "moose.tmin1" "elk.tmin1"   "wolf.tmin1"  "coy.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("moose.latent", "elk.latent", "wolf.latent", "coy.latent", "wtd.latent", "wolf.latent"), spp = c(".moose", ".elk", ".wolf", ".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #83
    # bs_bottomup_inter[[83]] # "wolf.t"      "lion.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"  "wtd.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "elk.latent", "wolf.latent", "wtd.latent", "lion.latent", "wolf.latent"), spp = c(".moose", ".elk", ".wolf", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #84
    # bs_bottomup_inter[[84]] # "wtd.t"        "coy.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "coy.tmin1"    "wolf.tmin1"
    list(dSep_test = 4, covariates = c("wsi", "forest", "wtd.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".wsi", ".forest", ".wtd", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #85
    # bs_bottomup_inter[[85]] # "wtd.t"        "lion.t"       "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "wolf.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("wsi", "forest", "wtd.latent", "wolf.latent", "lion.latent", "wtd.latent"), spp = c(".wsi", ".forest", ".wtd", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #86
    # bs_bottomup_inter[[63]] # "coy.t"      "lion.t"     "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 1, covariates = c("coy.latent", "wolf.latent", "wtd.latent", "lion.latent", "coy.latent"), spp = c(".coy", ".wolf", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y"))
  )
  
  