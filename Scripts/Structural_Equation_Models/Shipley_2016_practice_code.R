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
  
  #'  Fitting a measurement model when you have a latent variable
  #'  Pretend to have a data set with 4 imperfect measurements (observed indicators) 
  #'  of an unobserved latent variable
  set.seed(1)
  string.length <- rnorm(100, 0, 10)
  ruler.cent <- 1*string.length + rnorm(100, 0, 0.5)
  hand.lengths <- 0.07*string.length + rnorm(100, 0, 7)
  ruler.inches <- 0.39*string.length + rnorm(100, 0, 3.3)
  eyeball.guess <- 1*string.length + rnorm(100, 0, 10)
  strings <- data.frame(ruler.cent = ruler.cent, hand.lengths = hand.lengths, 
                         ruler.inches = ruler.inches, eyeball.guess = eyeball.guess)
  #'  =~ operator is used instead of ~, meaning the latent variable causes the indicators
  #'  NOTE: string.length is unobserved so not included in any df; this is how lavaan knows it is latent
  #'  Fix scale of the latent variable to unity so it is measured in units of standard deviations
  strings.model <- "
  string.length =~ ruler.cent + hand.lengths + ruler.inches + eyeball.guess"
  #'  Fix variance of latent variable to unity with std.lv = TRUE
  fit <- sem(model = strings.model, data = strings, std.lv = TRUE)
  summary(fit)
  
  #'  If you want to fix latent variable to something other than unity or fix some 
  #'  latent variabces to unity but fix scales of other latent variables using the
  #'  indicator, need to use syntax that says to not fix the path coefficient of first indicator variable & fix the variance of the latent
  strings.model <- "
  string.length =~ NA*ruler.cent + hand.lengths + ruler.inches + eyeball.guess 
  string.length ~~1*string.length" 
  fit <- sem(model = strings.model, data = strings)
  summary(fit)
  #' NA*ruler.cent allows model to estimate path coefficient, otherwise if std.lv 
  #' is not set to TRUE then package default fies the scale of the latent variable 
  #' by fixing the first path coefficient to unity
  #' string.length ~~1*string.length explicitly fixes variance of latent variable to unity
  
  #'  If you want to fix scale fo latent variable to that of first indicator variable 
  #'  you can allow variance of latent variable to be freely estimated but fix
  #'  path coefficient to first indicator variable to 1. This is the default process
  #'  with the most basic syntax for measurement models
  strings.model <- "
  string.length =~ ruler.cent + hand.lengths + ruler.inches + eyeball.guess"
  fit <- sem(model = strings.model, data = strings)
  summary(fit)
  #'  Get total proportion of total variance of each indicator that is explained by
  #'  latent value (highest correlation with ruler.cent, lowest with hand.lengths)
  inspect(fit, "r2") 
  predict(fit)
  