  #'  --------------------------
  #'  Practicing SEMs
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  January 2024
  #'  --------------------------
  #'  Working through example code in the Shipley 2016 "Cause and Correlation in 
  #'  Biology" book and online tutorials at https://lavaan.ugent.be/
  #'  --------------------------
  
  #install.packages("lavaan", dependencies = TRUE)
  library(lavaan)
  
  #'  Simulate data for causal network:
  #'  A affects B --> B affects C & D --> C, D, & F affect E
  #'  Written as a DAG:
  #'  A --> B --> C 
  #'                --> E <-- F
  #'          --> D 
  #'  Generate data for the two exogenous observed variables
  set.seed(1)
  var.A <- rnorm(100)
  var.F <- rnorm(100)
  
  #'  Generate data for the four endogenous observed variables
  #'  Incorporate variance terms and causal relationships with other variable(s)
  var.B <- 0.5*var.A + rnorm(100, 0, sqrt(1-0.5^2))
  var.C <- -0.5*var.B + rnorm(100, 0, sqrt(1-0.5^2))
  var.D <- 0.5*var.B + rnorm(100, 0, sqrt(1-0.5^2))
  var.E <- 0.5*var.C + 0.5*var.D - 0.5*var.F + rnorm(100, 0, sqrt(1-3*0.5^2))
  sim.data <- data.frame(A = var.A, B = var.B, C = var.C, D = var.D, E = var.E, F = var.F)
  
  #'  Define true model based on data generation process (effect ~ cause)
  true.model1 <- "
  B~A
  C~B
  D~B
  E~C+D+F
  "
  
  #'  Fit true model
  (fit.true <- sem(true.model1, data = sim.data))
  
  #'  Review which free parameters have been included in the model
  #'  ~~ represents variance or covariance
  #'  The user column indicates which parameters were specified by me (1) or 
  #'  implicitly specified as a default (0)
  parTable(fit.true) 
  #'  Note: A ~~ F should not be included b/c did not hypothesize covariance btwn A & F
  #'  If their covariance is estimated, this screws up the degrees of freedom
  
  #'  Refit the true model and explicitly set covariance btwn A & F to 0 (no free covariance)
  true.model2 <- "
  B~A
  C~B
  D~B
  E~C+D+F
  A~~0*F # fix path coefficient between A & F to 0
  "
  (fit.true2 <- sem(true.model2, data = sim.data, fixed.x = FALSE))
  
  #'  Review model outputs
  #'  Parameter estimates (path coefficients)
  summary(fit.true2)
  summary(fit.true2, fit.measures = TRUE) # gives more
  #'  Proportion of the variance of each endogenous variable that is captured by the path model (R^2)
  inspect(fit.true2, "rsquare")
  #'  Parameter estimates with bootstrapped 95% CIs and as standardized or unstandardized coefficients
  #'  Standardized coefficients are calculated using z-transformed variables
  parameterEstimates(fit.true2, ci = TRUE, level = 0.95, boot.ci.type = "perc", 
                     standardized = FALSE)
  #'  Predicted covariance matrix (to be compared to the observed covariance matrix)
  fitted(fit.true2)
  #'  Difference between the observed and predicted covariance matrices
  #'  (works best when estimates are based on standardized values)
  residuals(fit.true2)
  
  #'  Provide better initial values for the model
  true.model3 <- "
  #'  Causal relationships
  B~A
  C~B
  D~B
  E~C+D+F
  #'  Fixed covariance
  A~~0*F
  #'  Exogenous variances
  A~~start(1)*A
  F~~start(1)*F
  #'  Endogenous residual variance
  B~~start(0.7)*B
  C~~start(0.7)*C
  D~~start(0.7)*D
  E~~start(0.3)*E
  "
  (fit.true3 <- sem(true.model3, data = sim.data, fixed.x = FALSE))
  
  #'  Estimate the combined indirect and total effect of A on E through its various pathways
  #'  Adding *xx* in the syntax names the causal relationship so that it can easily be 
  #'  used elsewhere in the analysis, e.g., B~start(0.5)*c1*A names the B~A relationship "c1"
  true.model4 <- "
  #'  Causal relationships
  B~start(0.5)*c1*A          
  C~start(-0.5)*c2*B
  D~start(0.5)*c4*B
  E~start(0.5)*c3*C+start(0.5)*c5*D+start(-0.5)*F
  #'  Fixed covariance
  A~~0*F
  #'  Exogenous variances
  A~~start(1)*A
  F~~start(1)*F
  #'  Endogenous residual variance
  B~~start(0.7)*B
  C~~start(0.7)*C
  D~~start(0.7)*D
  E~~start(0.3)*E
  #'  Definte indirect, total effects
  indirect1:=c1*c2*c3
  indirect2:=c1*c4*c5
  total:=indirect1+indirect2
  "
  fit.true4 <- sem(true.model4, data = sim.data, fixed.x = FALSE)
  summary(fit.true4)
  
  
  