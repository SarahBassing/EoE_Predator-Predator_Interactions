  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations for cor(exog, exog): top-down
  #'  Aug 2026
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
  
  dSep_iterations_topdown_exog_only <- list(
    # Regression #3
    # bs_topdown[[3]] # "deerHarv.tmin1" "elkHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$elkHarv, x_array = data_JAGS_bundle_topdown$deerHarv),   
    # Regression #8
    # bs_topdown[[8]] # "deerHarv.tmin1" "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$bearHarv, x_array = data_JAGS_bundle_topdown$deerHarv),   
    # Regression #11
    # bs_topdown[[11]] # "deerHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$wolfHarv, x_array = data_JAGS_bundle_topdown$deerHarv),   
    # Regression #16
    # bs_topdown[[16]] # "deerHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$lionHarv, x_array = data_JAGS_bundle_topdown$deerHarv),   
    # Regression #54
    # bs_topdown[[54]] # "elkHarv.tmin1"  "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$bearHarv, x_array = data_JAGS_bundle_topdown$elkHarv),   
    # Regression #57
    # bs_topdown[[57]]  # "elkHarv.tmin1"  "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$wolfHarv, x_array = data_JAGS_bundle_topdown$elkHarv),   
    # Regression #61
    # bs_topdown[[61]]  # "elkHarv.tmin1"  "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$lionHarv, x_array = data_JAGS_bundle_topdown$elkHarv),   
    # Regression #120
    # bs_topdown[[120]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$wolfHarv, x_array = data_JAGS_bundle_topdown$bearHarv),   
    # Regression #125
    # bs_topdown[[125]] # "bearHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$lionHarv, x_array = data_JAGS_bundle_topdown$bearHarv),   
    # Regression #144
    # bs_topdown[[144]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown$wolfHarv, x_array = data_JAGS_bundle_topdown$lionHarv)    
  )
  
  
  