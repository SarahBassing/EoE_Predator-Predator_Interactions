  #'  ---------------------------------
  #'  Lefcheck piecewiseSEM practice
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  January 2024
  #'  ---------------------------------
  #'  Familiarizing myself with the piecewiseSEM package following "An Introduction
  #'  to Structural Equation Modeling in R" at https://jslefche.github.io/sem_book/
  #'  ---------------------------------
  
  #install.packages("piecewiseSEM")
  library(piecewiseSEM)
  data(keeley)
  
  #'  Example DAG
  #'  age --> firesev --> cover
  
  #'  Define the list of structural equations, written in linear modeling format
  keeley_psem <- psem(
    lm(cover ~ firesev, data = keeley),
    lm(firesev ~ age, data = keeley),
    data = keeley
  )
  keeley_psem
  
  #'  Derive basis set (identify all claims of conditional independence)
  #'  "a|b(c)" implies a & b are statistically independent, conditioned on c
  basisSet(keeley_psem) # "age | cover ( firesev )"
  
  #'  Conduct test of d-separation
  dSep(keeley_psem, .progressBar = FALSE)
  #'  Output should be the same as if we evaluated the independence claim ourselves
  summary(lm(cover ~ firesev + age, data = keeley))$coefficients[3, ]
  #'  Since p-value is > 0.05, we fail to reject the null, implying the claim of
  #'  conditional independence was reasonable given alpha value = 0.05
  
  #'  Compute Fisher's C statistic using p-values obtained from the d-sep tests
  #'  Fisher's C degrees of freedom = 2k, where k = number of independence claims; df = 2 * 1 = 2
  #'  C = -2 * sum(ln(p[i])) where i is the ith independence claim, 
  #'  p is the p-value from the ith significance test, which is logged using the natural logarithm,
  #'  all of which are summed across all independence claims in the basis set and multiplied by -2
  (C <- -2 * log(summary(lm(cover ~ firesev + age, data = keeley))$coefficients[3,4]))
  df <- 2 * 1
  #'  Calculate p-value by comparing Fisher's C to the Chi-sqrd table with appropriate df
  1-pchisq(C, df)
  #'  0.075 > 0.05 so fail to reject the null. The model and independence claims are supported.
  
  #'  Easier way to compute and compare the Fisher's C test statistic
  fisherC(keeley_psem)
  
  #'  Alternatively, using log-likelihood based Chi-sqrd statistic to evaluate the model
  #'  First create the saturated model, which includes all possible direct paths among variables
  #'  cover ~ firesev + age syntax implies there are direct effects of firesev on cover and age on cover
  keeley_psem2 <- psem(
    lm(cover ~ firesev + age, data = keeley),
    lm(firesev ~ age, data = keeley),
    data = keeley
  )
  
  #'  Compare log-likelihoods of each model by calculating the log-likelihood of 
  #'  each submodel in the original model and subtracting them from the log-likelihood 
  #'  of their corresponding submodel in the saturated model
  LL_1 <- logLik(lm(cover ~ firesev, data = keeley)) - logLik(lm(cover ~ firesev + age, data = keeley))
  LL_2 <- logLik(lm(firesev ~ age, data = keeley)) - logLik(lm(firesev ~ age, data = keeley))
  
  #'  Construct the Chi-sqrd test statistic from the difference of these log-likelihoods
  (ChiSq <- -2 * sum(as.numeric(LL_1), as.numeric(LL_2)))
  
  DF <- 1 # one additional parameter was estimated in the saturated model
  1 - pchisq(ChiSq, DF)
  
  #'  Easier way to conduct the log-likelihood goodness-of-fit test
  LLchisq(keeley_psem)
  #'  P = 0.069 is > 0.05 so again, we fail to reject the model
  #'  Note: this test statistic and p-value are the same as what's obtained in 
  #'  the lavaan package b/c this model assumes multivariate normality given the data
  
  #'  Obtain an AIC score for the model
  #'  Note: default is based on the log-likelihood Chi-sqrd method
  AIC(keeley_psem)
  #'  Get AIC score based on Fisher's C statistic and d-sep tests
  AIC(keeley_psem, AIC.type = "dsep")
  #'  Note: this approach cannot be used if we want to compare the model to the 
  #'  saturated model b/c we cannot run a GoF test on the full saturated (or just 
  #'  identified) model because there are no degrees of freedom
  
  #'  Finally, all of this can be computed at once with the summary() function
  summary(keeley_psem, .progressBar = FALSE)
  #'  NOTE: standardized coefficients are returned with the piecewiseSEM package
  #'  Individual model R^2 values are also returned
  
  
  #'  -------------------------------------------------------
  ####  piecewise SEM with generalized mixed effects models  ####
  #'  -------------------------------------------------------
  data(shipley)
  
  