  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down interference
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
  #'  d_Sep_active_regression_topdown_inter_tmin1_only.R  and 
  #'  d_Sep_active_regression_topdown_inter_exog_only.R
  #'  ---------------------------------------------------------------- 

  dSep_iterations_topdown_int <- list(
    # Regression #12
    # bs_topdown_inter[[12]] # "deerHarv.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "deerHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #13
    # bs_topdown_inter[[13]] # "deerHarv.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "deerHarv"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #14
    # bs_topdown_inter[[14]] # "deerHarv.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "deerHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #15
    # bs_topdown_inter[[15]] # "deerHarv.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "deerHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #16
    # bs_topdown_inter[[16]] # "deerHarv.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "deerHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #17
    # bs_topdown_inter[[17]] # "deerHarv.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "deerHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #28
    # bs_topdown_inter[[28]] # "wtd.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #29
    # bs_topdown_inter[[29]] # "wtd.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #30
    # bs_topdown_inter[[30]] # "wtd.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #31
    # bs_topdown_inter[[31]] # "wtd.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #32
    # bs_topdown_inter[[32]] # "wtd.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #33
    # bs_topdown_inter[[33]] # "wtd.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #37
    # bs_topdown_inter[[37]] # "moose.tmin1" "wtd.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #44
    # bs_topdown_inter[[44]] # "moose.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #45
    # bs_topdown_inter[[45]] # "moose.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "moose.latent"), spp = c(".lion", ".coy", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #46
    # bs_topdown_inter[[46]] # "moose.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #47
    # bs_topdown_inter[[47]] # "moose.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #48
    # bs_topdown_inter[[48]] # "moose.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #51
    # bs_topdown_inter[[51]] # "elkHarv.tmin1" "wtd.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #58
    # bs_topdown_inter[[58]] # "elkHarv.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elkHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #59
    # bs_topdown_inter[[59]] # "elkHarv.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "elkHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #60
    # bs_topdown_inter[[60]] # "elkHarv.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "elkHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #61
    # bs_topdown_inter[[61]] # "elkHarv.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elkHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #62
    # bs_topdown_inter[[62]] # "elkHarv.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "elkHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #64
    # bs_topdown_inter[[64]] # "elk.tmin1" "wtd.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elk.latent"), spp = c(".harvest", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #71
    # bs_topdown_inter[[71]] # "elk.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #72
    # bs_topdown_inter[[72]] # "elk.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "elk.latent"), spp = c(".lion", ".coy", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #73
    # bs_topdown_inter[[73]] # "elk.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".harvest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #74
    # bs_topdown_inter[[74]] # "elk.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #75
    # bs_topdown_inter[[75]] # "elk.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #82
    # bs_topdown_inter[[82]] # "lion.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #83
    # bs_topdown_inter[[83]] # "lion.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lion.latent"), spp = c(".harvest", ".bear", ".wolf", ".lion"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #84
    # bs_topdown_inter[[84]] # "lion.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #85
    # bs_topdown_inter[[85]] # "lion.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #86
    # bs_topdown_inter[[86]] # "wtd.t" "coy.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "coy.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #87
    # bs_topdown_inter[[87]] # "wtd.t" "bearHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bearHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #88
    # bs_topdown_inter[[88]] # "wtd.t" "bear.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bear.latent"), spp = c(".harvest", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #89
    # bs_topdown_inter[[89]] # "wtd.t" "wolfHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolfHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #90
    # bs_topdown_inter[[90]] # "wtd.t" "lionHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "lionHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #91
    # bs_topdown_inter[[91]] # "wtd.t" "wolf.tmin1" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1"
    # NOTE: FLIPPED so the .t variable (wtd.t) is the response
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolf.latent"), spp = c(".harvest", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #92
    # bs_topdown_inter[[92]] # "wtd.t" "moose.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "moose.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = moose.t
    list(dSep_test = 6, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #93
    # bs_topdown_inter[[93]] # "wtd.t" "elk.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "elkHarv.tmin1" "elk.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = elk.t
    list(dSep_test = 5, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #94
    # bs_topdown_inter[[94]] # "wtd.t" "coy.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = coy.t
    list(dSep_test = 4, covariates = c("deerHarv", "wtd.latent", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #95
    # bs_topdown_inter[[95]] # "wtd.t" "bear.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = bear.t
    list(dSep_test = 3, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #96
    # bs_topdown_inter[[96]] # "wtd.t" "wolf.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "wolfHarv.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = wolf.t
    list(dSep_test = 2, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #97
    # bs_topdown_inter[[97]] # "wtd.t" "lion.t" "deerHarv.tmin1" "wtd.tmin1" "lion.tmin1" "lionHarv.tmin1" "wolf.tmin1"
    # NOTE: both .t -> tested term "wtd.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("deerHarv", "wtd.latent", "lion.latent", "lionHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #103
    # bs_topdown_inter[[103]] # "coy.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #104
    # bs_topdown_inter[[104]] # "coy.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #105
    # bs_topdown_inter[[105]] # "coy.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #106
    # bs_topdown_inter[[106]] # "coy.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #107
    # bs_topdown_inter[[107]] # "coy.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #112
    # bs_topdown_inter[[112]] # "bearHarv.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #113
    # bs_topdown_inter[[113]] # "bearHarv.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "bearHarv"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #114
    # bs_topdown_inter[[114]] # "bearHarv.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #115
    # bs_topdown_inter[[115]] # "bearHarv.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #116
    # bs_topdown_inter[[116]] # "bearHarv.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #120
    # bs_topdown_inter[[120]] # "bear.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #121
    # bs_topdown_inter[[121]] # "bear.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #122
    # bs_topdown_inter[[122]] # "bear.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bear.latent"), spp = c(".lion", ".coy", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #123
    # bs_topdown_inter[[123]] # "bear.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #124
    # bs_topdown_inter[[124]] # "bear.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #127
    # bs_topdown_inter[[127]] # "wolfHarv.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #128
    # bs_topdown_inter[[128]] # "wolfHarv.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #129
    # bs_topdown_inter[[129]] # "wolfHarv.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #130
    # bs_topdown_inter[[130]] # "wolfHarv.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #131
    # bs_topdown_inter[[131]] # "wolfHarv.tmin1" "lion.t" "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #133
    # bs_topdown_inter[[133]] # "lionHarv.tmin1" "moose.t" "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1", "y-1", "y-1")),
    # Regression #134
    # bs_topdown_inter[[134]] # "lionHarv.tmin1" "elk.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1")),
    # Regression #135
    # bs_topdown_inter[[135]] # "lionHarv.tmin1" "coy.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #136
    # bs_topdown_inter[[136]] # "lionHarv.tmin1" "bear.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1")),
    # Regression #137
    # bs_topdown_inter[[137]] # "lionHarv.tmin1" "wolf.t" "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lionHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1", "y-1", "y-1")),
    # Regression #138
    # NOTE: both .t -> tested term "moose.t" enters as CONTEMPORANEOUS (lag "y"), outcome = elk.t
    # bs_topdown_inter[[138]] # "moose.t" "elk.t" "moose.tmin1" "wolf.tmin1" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1"
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elkHarv", "elk.latent", "lion.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #139
    # bs_topdown_inter[[139]] # "moose.t" "coy.t" "moose.tmin1" "wolf.tmin1" "lion.tmin1" "coy.tmin1"
    # NOTE: both .t -> tested term "moose.t" enters as CONTEMPORANEOUS (lag "y"), outcome = coy.t
    list(dSep_test = 4, covariates = c("moose.latent", "wolf.latent", "lion.latent", "coy.latent", "moose.latent"), spp = c(".moose", ".wolf", ".lion", ".coy", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #140
    # bs_topdown_inter[[140]] # "moose.t" "bear.t" "moose.tmin1" "wolf.tmin1" "bearHarv.tmin1" "bear.tmin1"
    # NOTE: both .t -> tested term "moose.t" enters as CONTEMPORANEOUS (lag "y"), outcome = bear.t
    list(dSep_test = 3, covariates = c("moose.latent", "wolf.latent", "bearHarv", "bear.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #141
    # bs_topdown_inter[[141]] # "moose.t" "wolf.t" "moose.tmin1" "wolf.tmin1" "wolfHarv.tmin1"
    # NOTE: both .t -> tested term "moose.t" enters as CONTEMPORANEOUS (lag "y"), outcome = wolf.t
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y")),
    # Regression #142
    # bs_topdown_inter[[142]] # "moose.t" "lion.t" "moose.tmin1" "wolf.tmin1" "lionHarv.tmin1"
    # NOTE: both .t -> tested term "moose.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y")),
    # Regression #143
    # bs_topdown_inter[[143]] # "elk.t" "coy.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1" "coy.tmin1"
    # NOTE: both .t -> tested term "elk.t" enters as CONTEMPORANEOUS (lag "y"), outcome = coy.t
    list(dSep_test = 4, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "coy.latent", "elk.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".coy", ".elk"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #144
    # bs_topdown_inter[[144]] # "elk.t" "bear.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1" "bearHarv.tmin1" "bear.tmin1"
    # NOTE: both .t -> tested term "elk.t" enters as CONTEMPORANEOUS (lag "y"), outcome = bear.t
    list(dSep_test = 3, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "bearHarv", "bear.latent", "elk.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1,1,2,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #145
    # bs_topdown_inter[[145]] # "elk.t" "wolf.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1" "wolfHarv.tmin1"
    # NOTE: both .t -> tested term "elk.t" enters as CONTEMPORANEOUS (lag "y"), outcome = wolf.t
    list(dSep_test = 2, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "wolfHarv", "elk.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #146
    # bs_topdown_inter[[146]] # "elk.t" "lion.t" "elkHarv.tmin1" "elk.tmin1" "lion.tmin1" "wolf.tmin1" "lionHarv.tmin1"
    # NOTE: both .t -> tested term "elk.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("elkHarv", "elk.latent", "lion.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".harvest", ".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #147
    # bs_topdown_inter[[147]] # "coy.t" "bear.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1" "bearHarv.tmin1" "bear.tmin1"
    # NOTE: both .t -> tested term "coy.t" enters as CONTEMPORANEOUS (lag "y"), outcome = bear.t
    list(dSep_test = 3, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #148
    # bs_topdown_inter[[148]] # "coy.t" "wolf.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1" "wolfHarv.tmin1"
    # NOTE: both .t -> tested term "coy.t" enters as CONTEMPORANEOUS (lag "y"), outcome = wolf.t
    list(dSep_test = 2, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #149
    # bs_topdown_inter[[149]] # "coy.t" "lion.t" "lion.tmin1" "coy.tmin1" "wolf.tmin1" "lionHarv.tmin1"
    # NOTE: both .t -> tested term "coy.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #150
    # bs_topdown_inter[[150]] # "bear.t" "wolf.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1" "wolfHarv.tmin1"
    # NOTE: both .t -> tested term "bear.t" enters as CONTEMPORANEOUS (lag "y"), outcome = wolf.t
    list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #151
    # bs_topdown_inter[[151]] # "bear.t" "lion.t" "bearHarv.tmin1" "bear.tmin1" "wolf.tmin1" "lionHarv.tmin1"
    # NOTE: both .t -> tested term "bear.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1", "y-1", "y-1", "y-1", "y")),
    # Regression #152
    # bs_topdown_inter[[152]] # "wolf.t" "lion.t" "wolfHarv.tmin1" "wolf.tmin1" "lionHarv.tmin1"
    # NOTE: both .t -> tested term "wolf.t" enters as CONTEMPORANEOUS (lag "y"), outcome = lion.t
    list(dSep_test = 1, covariates = c("wolfHarv", "wolf.latent", "lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf", ".harvest", ".wolf"), indices = as.integer(c(1,1,2,2)), lags = c("y-1", "y-1", "y-1", "y")))
  
    
    
    
    
    
    
    
    
    
    
    
    
    
  #   # # Regression #1
  #   # # bs_topdown_inter[[1]] # "wtd.tmin1"   "moose.tmin1"
  #   # list(dSep_test = 5, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # # Regression #2
  #   # # bs_topdown_inter[[2]] # "wtd.tmin1" "elk.tmin1"
  #   # list(dSep_test = 4, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # # Regression #3
  #   # # bs_topdown_inter[[3]] # "wtd.tmin1"  "lion.tmin1"
  #   # list(dSep_test = 1, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # # Regression #4
  #   # # bs_topdown_inter[[4]] # "wtd.tmin1" "coy.tmin1"
  #   # list(dSep_test = 3, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #5         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[5]] # "wtd.tmin1"      "bearHarv.tmin1"
  #   # list(dSep_test = 7, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #6
  #   # # bs_topdown_inter[[6]] # "wtd.tmin1"  "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #7         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[7]] # "wtd.tmin1"      "wolfHarv.tmin1"
  #   # list(dSep_test = 7, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #8         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[8]] # "wtd.tmin1"      "lionHarv.tmin1"
  #   # list(dSep_test = 7, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #9
  #   # # bs_topdown_inter[[9]] # "wtd.tmin1"  "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("wtd.latent"), spp = c(".wtd"), indices = c(1), lags = c("y-1")),
  #   # Regression #10
  #   # bs_topdown_inter[[10]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #11
  #   # bs_topdown_inter[[11]] # "wtd.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wtd.latent"), spp = c(".elk", ".lion", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #12
  #   # bs_topdown_inter[[12]] # "wtd.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #13
  #   # bs_topdown_inter[[13]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #14
  #   # bs_topdown_inter[[14]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #15
  #   # bs_topdown_inter[[15]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #16
  #   # # bs_topdown_inter[[16]] # "moose.tmin1" "elk.tmin1"
  #   # list(dSep_test = 5, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
  #   # # Regression #17
  #   # # bs_topdown_inter[[17]] # "moose.tmin1" "lion.tmin1"
  #   # list(dSep_test = 1, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
  #   # Regression #18
  #   # bs_topdown_inter[[18]] # "moose.tmin1" "wtd.t"       "wtd.tmin1"   "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "moose.latent"), spp = c(".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #19
  #   # # bs_topdown_inter[[19]] # "moose.tmin1" "coy.tmin1"
  #   # list(dSep_test = 3, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #20         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[20]] # "moose.tmin1"    "bearHarv.tmin1"
  #   # list(dSep_test = 6, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #21
  #   # # bs_topdown_inter[[21]] # "moose.tmin1" "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("moose.latent"), spp = c(".moose"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #22         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[22]] # "moose.tmin1"    "wolfHarv.tmin1"
  #   # list(dSep_test = 6, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # #########
  #   # # Regression #23
  #   # # bs_topdown_inter[[23]] # "moose.tmin1"    "lionHarv.tmin1"
  #   # list(dSep_test = 6, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #24
  #   # # bs_topdown_inter[[24]] # "moose.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("wolf.latent"), spp = c(".wolf"), indices = c(1), lags = c("y-1")),
  #   # Regression #25
  #   # bs_topdown_inter[[25]] # "moose.tmin1" "elk.t"       "elk.tmin1"   "lion.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 5, covaraites = c("elk.latent", "lion.latent", "wolf.latent", "moose.latent"), spp = c(".elk", ".lion", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #26
  #   # bs_topdown_inter[[26]] # "moose.tmin1" "coy.t"       "lion.tmin1"  "coy.tmin1"   "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "moose.latent"), spp = c(".lion", ".coy", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #27
  #   # bs_topdown_inter[[27]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #28
  #   # bs_topdown_inter[[28]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #29
  #   # bs_topdown_inter[[29]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #30
  #   # # bs_topdown_inter[[30]] # "elk.tmin1"  "lion.tmin1"
  #   # list(dSep_test = 1, covariates = c("elk.latent"), spp = c(".elk"), indices = c(1), lags = c("y-1")),
  #   # Regression #31
  #   # bs_topdown_inter[[31]] # "elk.tmin1"  "wtd.t"      "wtd.tmin1"  "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "elk.latent"), spp = c(".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #32
  #   # # bs_topdown_inter[[32]] # "elk.tmin1" "coy.tmin1"
  #   # list(dSep_test = 3, covariates = c("elk.latent"), spp = c(".elk"), indices = c(1), lags = c("y-1")),
  #   # # Regression #33         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[33]] # "elk.tmin1"      "bearHarv.tmin1"
  #   # list(dSep_test = 5, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #34
  #   # # bs_topdown_inter[[34]] # "elk.tmin1"  "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("elk.latent"), spp = c(".elk"), indices = c(1), lags = c("y-1")),
  #   # # Regression #35         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[35]] # "elk.tmin1"      "wolfHarv.tmin1"
  #   # list(dSep_test = 5, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #36         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[36]] # "elk.tmin1"      "lionHarv.tmin1"
  #   # list(dSep_test = 5, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #37
  #   # # bs_topdown_inter[[37]] # "elk.tmin1"  "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("elk.latent"), spp = c(".elk"), indices = c(1), lags = c("y-1")),
  #   # Regression #38
  #   # bs_topdown_inter[[38]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #39
  #   # bs_topdown_inter[[39]] "elk.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "elk.latent"), spp = c(".lion", ".coy", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #40
  #   # bs_topdown_inter[[40]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".harvest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #41
  #   # bs_topdown_inter[[41]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #42
  #   # bs_topdown_inter[[42]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #43
  #   # # bs_topdown_inter[[43]] # "lion.tmin1" "coy.tmin1"
  #   # list(dSep_test = 4, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
  #   # # Regression #44         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[44]] # "lion.tmin1"     "bearHarv.tmin1"
  #   # list(dSep_test = 1, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #45
  #   # # bs_topdown_inter[[45]] # "lion.tmin1" "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
  #   # # Regression #46         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[46]] # "lion.tmin1"     "wolfHarv.tmin1"
  #   # list(dSep_test = 1, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #47         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[47]] # "lion.tmin1"     "lionHarv.tmin1"
  #   # list(dSep_test = 1, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #48
  #   # # bs_topdown_inter[[48]] # "lion.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("lion.latent"), spp = c(".lion"), indices = c(1), lags = c("y-1")),
  #   # Regression #49
  #   # bs_topdown_inter[[49]] # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #50
  #   # bs_topdown_inter[[50]] # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lion.latent"), spp = c(".harvest", ".bear", ".wolf", ".lion"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #51
  #   # bs_topdown_inter[[51]] # "lion.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #52
  #   # bs_topdown_inter[[52]] # "lion.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #53         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[23]] # "wtd.t"      "coy.tmin1"  "wtd.tmin1"  "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "coy.latent"), spp = c(".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #54         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[54]] # "wtd.t"          "bearHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "bearHarv"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #55         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[55]] # "wtd.t"      "bear.tmin1" "wtd.tmin1"  "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "bear.latent"), spp = c(".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #56         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[56]] # "wtd.t"          "wolfHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "wolfHarv"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #57         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[57]] # "wtd.t"          "lionHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "lionHarv"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #58         #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
  #   # bs_topdown_inter[[28]] # "wtd.t"      "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
  #   list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "wolf.latent"), spp = c(".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #59
  #   # bs_topdown_inter[[59]] # "wtd.t"       "moose.t"     "wtd.tmin1"   "lion.tmin1"  "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("wtd.latent", "lion.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #60
  #   # bs_topdown_inter[[60]] # "wtd.t"      "elk.t"      "wtd.tmin1"  "lion.tmin1" "elk.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("wtd.latent", "lion.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regerssion #61
  #   # bs_topdown_inter[[61]] # "wtd.t"      "coy.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("wtd.latent", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #62
  #   # bs_topdown_inter[[62]] # "wtd.t"          "bear.t"         "wtd.tmin1"      "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("wtd.latent", "lion.latent", "bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
  #   # Regression #63
  #   # bs_topdown_inter[[63]] # "wtd.t"          "wolf.t"         "wtd.tmin1"      "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wtd.latent", "lion.latent", "wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #64
  #   # bs_topdown_inter[[64]] # "wtd.t"          "lion.t"         "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("wtd.latent", "lion.latent", "lionHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # # Regression #65         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[65]] # "coy.tmin1"      "bearHarv.tmin1"
  #   # list(dSep_test = 4, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #66
  #   # # bs_topdown_inter[[66]] # "coy.tmin1"  "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("coy.latent"), spp = c(".coy"), indices = c(1), lags = c("y-1")),
  #   # # Regression #67         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[67]] # "coy.tmin1"      "wolfHarv.tmin1"
  #   # list(dSep_test = 4, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #68         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[68]] # "coy.tmin1"      "lionHarv.tmin1"
  #   # list(dSep_test = 4, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #69
  #   # # bs_topdown_inter[[69]] # "coy.tmin1"  "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("coy.latent"), spp = c(".coy"), indices = c(1), lags = c("y-1")),
  #   # Regression #70
  #   # bs_topdown_inter[[70]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #71
  #   # bs_topdown_inter[[71]] # "coy.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent"), spp = c(".elk", ".lion", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #72
  #   # bs_topdown_inter[[72]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv",  "bear.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #73
  #   # bs_topdown_inter[[73]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #74
  #   # bs_topdown_inter[[74]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #75
  #   # # bs_topdown_inter[[75]] # "bearHarv.tmin1" "bear.tmin1"
  #   # list(dSep_test = 3, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #76
  #   # # bs_topdown_inter[[76]] # "bearHarv.tmin1" "wolfHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # # Regression #77
  #   # # bs_topdown_inter[[77]] # "bearHarv.tmin1" "lionHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # # Regression #78
  #   # # bs_topdown_inter[[78]] # "bearHarv.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("bearHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # Regression #79
  #   # bs_topdown_inter[[79]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #80
  #   # bs_topdown_inter[[80]] # "bearHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #81
  #   # bs_topdown_inter[[81]] # "bearHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #82
  #   # bs_topdown_inter[[82]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
  #   # Regression #83
  #   # bs_topdown_inter[[83]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #84         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[84]] # "bear.tmin1"     "wolfHarv.tmin1"
  #   # list(dSep_test = 3, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #85         # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT HARVEST
  #   # # bs_topdown_inter[[85]] # "bear.tmin1"     "lionHarv.tmin1"
  #   # list(dSep_test = 3, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # # Regression #86
  #   # # bs_topdown_inter[[86]] # "bear.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("bear.latent"), spp = c(".bear"), indices = c(1), lags = c("y-1")),
  #   # Regression #87
  #   # bs_topdown_inter[[87]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #88
  #   # bs_topdown_inter[[88]] # "bear.tmin1" "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bear.latent"), spp = c(".elk", ".lion", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #89
  #   # bs_topdown_inter[[89]] # "bear.tmin1" "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bear.latent"), spp = c(".lion", ".coy", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #90
  #   # bs_topdown_inter[[90]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #91
  #   # bs_topdown_inter[[91]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #92
  #   # # bs_topdown_inter[[92]] # "wolfHarv.tmin1" "lionHarv.tmin1"
  #   # list(dSep_test = NA, covariates = c(), spp = c(), indices = c(1), lags = c("y-1")),   # harvest effects harvest...
  #   # # Regression #93
  #   # # bs_topdown_inter[[93]] # "wolfHarv.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("wolfHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # Regression #94
  #   # bs_topdown_inter[[94]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #95
  #   # bs_topdown_inter[[95]] # "wolfHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #96
  #   # bs_topdown_inter[[96]] # "wolfHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #97
  #   # bs_topdown_inter[[97]] # "wolfHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #98
  #   # bs_topdown_inter[[98]] # "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
  #   # # Regression #99
  #   # # bs_topdown_inter[[99]] # "lionHarv.tmin1" "wolf.tmin1"
  #   # list(dSep_test = 2, covariates = c("lionHarv"), spp = c(".harvest"), indices = c(1), lags = c("y-1")),
  #   # Regression #100
  #   # bs_topdown_inter[[100]] # "lionHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
  #   list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
  #   # Regression #101
  #   # bs_topdown_inter[[101]] # "lionHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #102
  #   # bs_topdown_inter[[102]] # "lionHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"
  #   list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #103
  #   # bs_topdown_inter[[103]] # "lionHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
  #   list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
  #   # Regression #104
  #   # bs_topdown_inter[[104]] # "lionHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
  #   list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lionHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
  #   # Regression #105
  #   # bs_topdown_inter[[105]] # "moose.t"     "elk.t"       "moose.tmin1" "wolf.tmin1"  "elk.tmin1"   "lion.tmin1"
  #   list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elk.latent", "lion.latent", "moose.latent"), spp = c(".moose", ".wolf", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #106
  #   # bs_topdown_inter[[106]] # "moose.t"     "coy.t"       "moose.tmin1" "wolf.tmin1"  "lion.tmin1"  "coy.tmin1"
  #   list(dSep_test = 4, covariates = c("moose.latent", "wolf.latent", "lion.latent", "coy.latent", "moose.latent"), spp = c(".moose", ".wolf", ".lion", ".coy", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #107
  #   # bs_topdown_inter[[107]] # "moose.t"        "bear.t"         "moose.tmin1"    "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
  #   list(dSep_test = 3, covariates = c("moose.latent", "wolf.latent", "bearHarv", "bear.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #108
  #   # bs_topdown_inter[[108]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
  #   list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
  #   # Regression #109
  #   # bs_topdown_inter[[109]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1"
  #   list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
  #   # Regression #110
  #   # bs_topdown_inter[[110]] # "elk.t"      "coy.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1" "coy.tmin1"
  #   list(dSep_test = 4, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".coy", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #111
  #   # bs_topdown_inter[[111]] # "elk.t"          "bear.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
  #   list(dSep_test = 3, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv", "bear.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
  #   # Regression #112
  #   # bs_topdown_inter[[112]] # "elk.t"          "wolf.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1"
  #   list(dSep_test = 2, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #113
  #   # bs_topdown_inter[[113]] # "elk.t"          "lion.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
  #   list(dSep_test = 1, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #114
  #   # bs_topdown_inter[[114]] # "coy.t"          "bear.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
  #   list(dSep_test = 3, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
  #   # Regression #115
  #   # bs_topdown_inter[[115]] # "coy.t"          "wolf.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wolfHarv.tmin1"
  #   list(dSep_test = 2, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #116
  #   # bs_topdown_inter[[116]] # "coy.t"          "lion.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "lionHarv.tmin1"
  #   list(dSep_test = 1, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #66
  #   # bs_topdown_inter[[72]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1"
  #   list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #118
  #   # bs_topdown_inter[[118]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
  #   list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y")),
  #   # Regression #119
  #   # bs_topdown_inter[[119]] # "wolf.t"         "lion.t"         "wolfHarv.tmin1" "wolf.tmin1"     "lionHarv.tmin1"
  #   list(dSep_test = 1, covariates = c("wolfHarv", "wolf.latent", "lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf", ".harvest", ".wolf"), indices = c(1,1,2,2), lags = c("y-1","y-1","y-1","y"))
  # )
  

  