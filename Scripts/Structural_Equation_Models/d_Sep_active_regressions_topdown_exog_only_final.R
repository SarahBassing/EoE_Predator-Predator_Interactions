  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations for cor(exog, exog): top-down FINAL
  #'  Sep 2026
  #'  ----------------------------------------------------------------
  #'  List inputs for active regressions to iterate over for each d-separation 
  #'  test when testing correlation between two exogenous variables flagged in 
  #'  the basic set. Because these are predictor variables, we do not have to 
  #'  worry about dealing with the latent process or propagating uncertainty.
  #'  It is also worth noting that there is not causal ordering implied here. 
  #'  This is simply to test for correlations between these variables to be 
  #'  consistent with the dSep process and to incorporate into the Fisher's C 
  #'  statistic if desired.
  #'  
  #'  Input for each iteration are: 
  #'  y_array = observed predictor variable 
  #'  x_array = observed predictor variable 
  #'  
  #'  This list of inputs differs from the input structures of the other two active 
  #'  regression scripts because we are not building an entirely new model with 
  #'  each iteration here. Instead, we are simply testing the correlation between 
  #'  two exogenous variables included in the larger SEM. 
  #'  ----------------------------------------------------------------
  
  dSep_iterations_topdown_exog_only_final <- list(
    # Regression #1
    # bs_topdown_final[[1]] # "deerHarv.tmin1" "elkHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$elkHarv, x_array = data_JAGS_bundle_topdown_final$deerHarv),  
    # Regression #4
    # bs_topdown_final[[4]] # "deerHarv.tmin1" "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$bearHarv, x_array = data_JAGS_bundle_topdown_final$deerHarv),  
    # Regression #7
    # bs_topdown_final[[7]] # "deerHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$wolfHarv, x_array = data_JAGS_bundle_topdown_final$deerHarv),  
    # Regression #16
    # bs_topdown_final[[16]] # "deerHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$lionHarv, x_array = data_JAGS_bundle_topdown_final$deerHarv),  
    # Regression #20
    # bs_topdown_final[[20]] # "elkHarv.tmin1"  "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$bearHarv, x_array = data_JAGS_bundle_topdown_final$elkHarv),  
    # Regression #23
    # bs_topdown_final[[23]] # "elkHarv.tmin1"  "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$wolfHarv, x_array = data_JAGS_bundle_topdown_final$elkHarv),  
    # Regression #32
    # bs_topdown_final[[32]] # "elkHarv.tmin1"  "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$lionHarv, x_array = data_JAGS_bundle_topdown_final$elkHarv),  
    # Regression #63
    # bs_topdown_final[[63]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$wolfHarv, x_array = data_JAGS_bundle_topdown_final$bearHarv),  
    # Regression #72
    # bs_topdown_final[[72]] # "bearHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$lionHarv, x_array = data_JAGS_bundle_topdown_final$bearHarv),
    # Regression #103
    # bs_topdown_final[[103]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_final$lionHarv, x_array = data_JAGS_bundle_topdown_final$wolfHarv)
  )
  
  
