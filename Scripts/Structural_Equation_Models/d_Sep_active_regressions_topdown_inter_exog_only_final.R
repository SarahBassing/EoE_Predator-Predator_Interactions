  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations for cor(exog, exog): top-down interference FINAL
  #'  Sept 2026
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
  
  dSep_iterations_topdown_inter_exog_only_final <- list(
    # Regression #3
    # bs_topdown_inter_final[[3]] # "deerHarv.tmin1" "elkHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$elkHarv, x_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #7
    # bs_topdown_inter_final[[7]] # "deerHarv.tmin1" "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$bearHarv, x_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #8
    # bs_topdown_inter_final[[8]] # "deerHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$wolfHarv, x_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #10
    # bs_topdown_inter_final[[10]] # "deerHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$deerHarv),
    # Regression #52
    # bs_topdown_inter_final[[52]] # "elkHarv.tmin1" "bearHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$bearHarv, x_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #53
    # bs_topdown_inter_final[[53]] # "elkHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$wolfHarv, x_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #56
    # bs_topdown_inter_final[[56]] # "elkHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$elkHarv),
    # Regression #96
    # bs_topdown_inter_final[[96]] # "bearHarv.tmin1" "wolfHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$wolfHarv, x_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #99
    # bs_topdown_inter_final[[99]] # "bearHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
    # Regression #108
    # bs_topdown_inter_final[[108]] # "wolfHarv.tmin1" "lionHarv.tmin1"
    list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$wolfHarv)
  )
    
    
    
  #   # Regression #65
  #   # bs_topdown_inter_final[[65]] # "bearHarv.tmin1" "wolfHarv.tmin1"
  #   list(y_array = data_JAGS_bundle_topdown_inter_final$wolfHarv, x_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #66
  #   # bs_topdown_inter_final[[66]] # "bearHarv.tmin1" "lionHarv.tmin1"
  #   list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$bearHarv),
  #   # Regression #90
  #   # bs_topdown_inter_final[[90]] # "wolfHarv.tmin1" "lionHarv.tmin1"
  #   list(y_array = data_JAGS_bundle_topdown_inter_final$lionHarv, x_array = data_JAGS_bundle_topdown_inter_final$wolfHarv)
  # )