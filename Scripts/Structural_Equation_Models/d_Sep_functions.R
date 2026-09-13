  #'  -----------------------------------
  #'  D-separation functions
  #'  Sarah Bassing & Matt Falcy
  #'  September 2026
  #'  -----------------------------------
  #'  Series of functions to source to test conditional independence claims for 
  #'  each model using the directed-separation process. 
  #'  
  #'  Using two methods to conduct d-Sep tests. 
  #'  
  #'  1. Region of Practical Equivalence (ROPE) 
  #'  This is a bit of a hack to generate pseudo p-values for Fisher's C. Method 
  #'  adapted from Kruschke and Liddell (2018), who used "0.1" of 1 SD to define a ROPE:
  #'  https://link.springer.com/article/10.3758/s13423-016-1221-4
  #'  
  #'  2. Bayesian p-values
  #'  Bayesian p-values also used in Fisher's C goodness-of-fit test.
  #'  
  #'  *Coding assistance and trouble shooting conducted with the help of Claude.ai
  #'  -----------------------------------
  
  #'  -------------------------------
  ####  ROPE method for d-Sep tests  ####
  #'  -------------------------------
  #'  Define function to establish a region of practical equivalence (ROPE) around
  #'  the null value. This expresses a range of parameter values considered equivalent
  #'  to the null value (0). This range is -pct*sd(y) to pct*sd(y)
  #'  If 95% CRI falls entirely within this range, 95% of the posterior is practically
  #'  equivalent to the null value. Taking the mean calculates the proportion
  #'  of the posterior that falls inside the ROPE. In other words, what is the 
  #'  probability that the effect is small enough to be ignored? Larger p.rope 
  #'  values indicate a large percentage of the posterior is "practically" 0; 
  #'  these variables are conditionally independent. Small p.rope values (i.e., < 0.05) 
  #'  indicate a small percentage of the posterior falls inside the ROPE and thus  
  #'  it is not practically equivalent to 0 (i.e., it is "practically" significant).
  #'  These variables are NOT conditionally independent and this relationship is 
  #'  missing in the model. Another way to think about it is, by defining the ROPE 
  #'  as a small percentage of the standard deviation of y, we have defined a threshold 
  #'  of practical significance. If the variable of interest only moves the needle 
  #'  by a 10% of a standard deviation, we can consider its effect to be noise.
  #'  
  #'  ---> p.rope is the probability that the parameter is actually trivial given 
  #'  ---> the data, where small values indicate the prob of being trivial is low. 
  #'  
  #'  y represents the SD of the observed response variable
  #'  post represents the effect of the coefficient of interest - if the parameter 
  #'  effect is "meaningful" it should be larger than the typical "noise" or trivial
  #'  variation associated with the response variable
  p.rope <- function(pct = 0.1, y = NULL, post = NULL){
    #'  Ensure y is a vector
    y_vec <- as.vector(y)
    #'  Calculate SD, after removing any NAs in y
    sd_y <- sd(y_vec, na.rm = TRUE)
    #'  Average number of posterior draws that fall above lower ROPE value AND below upper ROPE value
    mean(-pct * sd_y < post & pct * sd_y > post)
    
  }
  
  #'  --------------------
  ####  Bayesian p-value  ####
  #'  --------------------
  #'  With a Bayesian p-value, the null hypothesis is that the coefficient of the 
  #'  "missing" variable being tested EQUALS zero (not some approximate region around
  #'  0, it IS 0), implying that path really is absent as the DAG asserts. A 
  #'  Bayesian p-value is asking what proportion of the posterior falls on each side
  #'  of zero. Under the formulation below, if the posterior is clustered away from 0, 
  #'  the tail that crosses 0 is small and the resulting p-value is small. If the 
  #'  posterior straddles 0 fairly symetrically, then ~half the posterio should be 
  #'  on either size of 0, giving a p-value close to 1, consistent with the null hypothesis. 
  bayes_pvalue <- function(post) {
    p_pos <- mean(post > 0)
    p_neg <- mean(post < 0)
    2 * min(p_pos, p_neg)
  }
  
  #'  Fisher's C using p-values 
  fishers_C <- function(pval, n_iter = NULL, K = NULL, n = NULL) {
    if(is.null(n_iter)) {
      eps <- 1e-10
    } else {
      eps <- 1 / n_iter # smallest resolvable p-value given draws
    }
    pval[pval <= 0] <- eps
    pval[pval >= 1] <- 1 - eps
    
    k <- length(pval)
    C <- -2 * sum(log(pval))
    df <- 2 * k
    pC <- 1 - pchisq(C, df)
    
    out <- list(claims.p = pval, k = k, C = C, df = df, p.value = pC)
    
    #'  Calculate Fisher's C AIC - must provide K for this part of function to work
    if(!is.null(K)) {
      out$K <- K
      out$AIC <- C + 2 * K
      
      if(!is.null(n)) {
        out$n <- n
        #'  Warn if small-sample correction is advisable (n/K < 40)
        #'  per Shipley (2013) and Burnham & Anderson convention
        if(n / K < 40) {
          message(sprintf("n/K = %.1f (< 40) -- AICc recommended over AIC for this model",
                          n/K))
        }
        out$AICc <- C + 2 * K * (n / (n - K - 1))
      }
    }
    
    return(out)
  }
  
  #'  ----------------------
  ####  Generate Basic set  ####
  #'  ----------------------
  #'  Function to generate and simplify the basic set 
  basic_set <- function(dag) {
    #'  Generate the basic set 
    basicset <- basiSet(dag)
    print(length(basicset))
    View(basicset)
    return(basicset)
  }
  
  #'  -------------------------------------------------------
  #####  Functions for when variable A t-1 --> variable B t  #####
  #'  These functions also work for when variable A t --> variable B t
  #'  -------------------------------------------------------
  #'  Source JAGS template
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_dsep_template.R")
  
  #'  Function to build custom regressions for each iteration
  build_individual_submodels <- function(reg_num, covariates = NULL, spp = NULL, 
                                         indices = NULL, lags = NULL) {   #, beta_prefix = "beta"
    
    
    #'  Create intercept and slope parameter names based on time step
    #'  beta_prefix indicates whether to use main model terms (beta.int) or the
    #'  auxilary beta terms (beta.aux). Only needed for d-Sep claims where t-1
    #'  affects t-1. This helps keep the rest of the main model beta arrays and
    #'  number of regressions unchanged.
    beta0_array <- "beta.int"    # intercept array 
    
    #'  If covariate is null or 0 (i.e., intercept only regressions)
    if(is.null(covariates) || length(covariates) == 0) {
      sprintf("%s[%d]", beta0_array, reg_num)
    } else {
      #'  Build one string per "beta * covariate" by filling placeholders with
      #'  specified character strings or integers (sprintf is vectorized so does
      #'  this in order that strings/integers are provided)
      #'  %s is placeholder for character strings; %d is placeholder for integers
      if(is.null(lags)) lags <- rep("y-1", length(covariates))
      if (length(lags) != length(covariates)) {
        stop("lags must be NULL or the same length as covariates")
      }
      
      terms <- paste0(sprintf("beta%s[%d] * %s[i,%s]", spp, indices, covariates, lags),
                      collapse = " + ")
      sprintf("%s[%d] + %s", beta0_array, reg_num, terms)
    }
  }
  
  #'  Function to assemble full model string for a single iteration
  #'  Requires having sourced the JAGS template model
  build_model_string <- function(iter_config, template, registry) {
    #'  Create empty strings to be filled with each regressions terms
    mu_lines <- character(7)  
    
    #'  Build custom regressions
    for(r in 1:7) {   # for each regression
      #'  If the regression index is the same as the dSep_test value then...
      if(r == iter_config$dSep_test) {
        covs <- iter_config$covariates
        spp <- iter_config$spp
        indices <- iter_config$indices
        lags <- iter_config$lags        # if NULL, defaults to "y-1
        mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags)
      } else {
        #'  Non-focal time t: use original SEM covariates from the registry
        orig <- registry[[r]]
        #'  Grab covariate, species, and index numbers of terms NOT included in d-sep test
        covs <- orig$covs
        spp <- orig$spp
        indices <- orig$indices
        lags <- orig$lags
        mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags)
      }
    }
    
    do.call(sprintf, c(list(template), as.list(mu_lines))) #all_lines
  }
  
  #'  Function to call JAGS and run a single iteration of the model
  run_dSep_iterations <- function(i, iterations, template, registry, data_bundle, listInits, model_name) {  
    
    iter_config <- iterations[[i]]
    
    #'  Call function to build the full model string with custom regressions for d-Sep test
    model_string <- build_model_string(iter_config, template, registry)
    
    #'  Create temporary directory to save all iterations of the template
    temp_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/temp_models", model_name)
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Create name new name for each model where %03d is a placeholder for the 
    #'  iteration number, padded by up to three 0's (e.g., temp_model_001, temp_model_012)
    temp_file <- file.path(temp_dir, sprintf("temp_model_%03d.txt", i))
    writeLines(model_string, temp_file)
    
    #'  Snag bundled data
    data_i <- data_bundle
    #' #'  Indicate which paramters to monitor
    #' monitor_params <- get_monitor_params(iter_config, params)
    
    #'  Fit model in JAGS
    SEM_dSep <- jagsUI::jags(data = data_i, inits = listInits, params, model.file = temp_file, 
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = FALSE, verbose = FALSE)
    # chain-level parallelism is off - parallelizing across
    # iterations instead
    
    #'  Flag convergence issues
    max_rhat <- suppressWarnings(max(unlist(SEM_dSep$Rhat), na.rm = TRUE))
    if(is.finite(max_rhat) && max_rhat > 1.1) {
      warning(sprintf("Iteration %d: max Rhat = %.3f -- check convergence", i, max_rhat))
    }
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = SEM_dSep, config = iter_config, max_rhat = max_rhat), jags_out)
    
    jags_out
  }
  
  #'  ---------------------------------------------------------
  #####  Functions for when variable A t-1 --> variable B t-1  #####
  #'  ---------------------------------------------------------
  #'  Function to extract the latent posterior mean and sd from the main model 
  #'  that was originally fitted. This requires that each spp.latent was monitored 
  #'  when the original model was fit.
  extract_latent_post_summaries <- function(og_fit, spp, nSites, nYear) {
    #'  og_fit: fitted jagsUI output from original SEM (fit in Bayesian_SEMs_relative_density_index_1yLag.R)
    #'  spp: species name - must match naming convention use din spp.latent array name
    
    #'  Grab posteriors
    samples_matrix <- as.matrix(og_fit$samples)
    latent_mean <- matrix(NA_real_, nSites, nYear)
    latent_sd <- matrix(NA_real_, nSites, nYear)
    
    for(i in 1:nSites) {
      for(y in 1:nYear) {
        #'  Create column name for specific spp.latent variable indexed by [nSite,nYear]
        col_name <- sprintf("%s.latent[%d,%d]", spp, i, y)
        #'  If created column name is in the colnames extracted from the og_fit samples
        if(col_name %in% colnames(samples_matrix)) {
          #'  Snag those draws from the posterior and save the mean and sd from each iteration
          draws <- samples_matrix[, col_name]
          latent_mean[i,y] <- mean(draws)
          latent_sd[i,y] <- sd(draws)
        }
      }
    }
    list(mean = latent_mean, sd = latent_sd)
  }
  
  #'  Function to fit a standalone regression using the posterior summaries as
  #'  the noisy "observed" response
  #'  Note: z = NULL is a placeholder for if any of these independence claims had
  #'  a condition set of variables. But in this case, z is never needed.
  fit_one_dSep_claim <- function(y, x, z = NULL, n.chains = n.chains, n.adapt = n.adapt, 
                                 n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin, 
                                 model_name, iter_num) {
    ncondvars <- if(is.null(z) || ncol(as.matrix(z)) == 0) 0 else ncol(as.matrix(z))
    #'  Create and fill in input data for JAGS (in list format)
    jd <- list(y = y, x = x, N = length(y))
    if(ncondvars > 0) {
      jd$z <- as.matrix(z)
      jd$ncondvars <- ncondvars
    }
    
    #'  Create custom regression with added variable for independence claims
    if(ncondvars > 0) {
      #'  Build conditioning variable terms
      z_terms <- paste(sprintf("b_z[%d] * z[i,%d]", 1:ncondvars, 1:ncondvars), collapse = " + ")
      
      #'  Create model string for JAGS that can be updated dynamically for each ind. claim
      #'  IF the independence claim includes a condition set of predictors (the 
      #'  "given blah blah blah" variables) create the first model_string. ELSE 
      #'  create the second model_string with only the y (spp.latent) and x (focal predictor)
      model_string <- sprintf("
                              model {
                              #'  Likelihood to be appended with conditioning claim
                              for(i in 1:N) {
                                y[i] ~ dnorm(mu[i], tau)
                                mu[i] <- b0 + b_x * x[i] + %s }
                              
                              #'  Priors
                              b0 ~ dnorm(0, 1e-4)
                              b_x ~ dnorm(0, 1e-4)
                                for(j in 1:ncondvars) {
                                  b_z[j] ~ dnorm(0, 1e-4) }
                              tau ~ dgamma(0.01, 0.01)
                              }
                              ", z_terms)
    } else {
      model_string <- "
      model {
      #'  Likelihood for all other regressions
      for(i in 1:N) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- b0 + b_x * x[i] }
      #'  Priors
      b0 ~ dnorm(0, 1e-4)
      b_x ~ dnorm(0, 1e-4)
      tau ~ dgamma(0.01, 0.01) 
      }
      "
    }
    
    #'  Create temporary directory to save all iterations of the template
    temp_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/temp_models/tmin1", model_name)
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Create name new name for each model where %03d is a placeholder for the 
    #'  iteration number, padded by up to three 0's (e.g., temp_model_001, temp_model_012)
    temp_file <- file.path(temp_dir, sprintf("temp_model_%03d.txt", iter_num))
    writeLines(model_string, temp_file)
    
    #'  Refit model with added independence claim
    fit <- jagsUI::jags(data = jd, inits = NULL, parameters.to.save = "b_x",
                        model.file = temp_file, n.chains = n.chains, n.adapt = n.adapt, 
                        n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin, 
                        parallel = FALSE, verbose = FALSE)
    
    fit
  }
  
  #'  Function tying this together by extracting posteriors and refitting the model
  fit_aux_claim <- function(i, iterations, og_fit, nSites, nYear, model_name,
                            n.chains, n.adapt, n.burnin, n.iter, n.thin) { #spp, covariate_array, 
    
    #'  Grab details for focal iteration
    iter_deets <- iterations[[i]]
    spp <- iter_deets$spp
    covariate_array <- iter_deets$covariate_array
    
    #'  Grab the posterior samples from the specified variable in the ind. claim
    post <- extract_latent_post_summaries(og_fit, spp, nSites, nYear)
    #'  Grab the spp.latent posterior mean and the "observed" predictor 
    y_vec <- as.vector(post$mean)
    x_vec <- as.vector(covariate_array)
    keep <- !is.na(y_vec) & !is.na(x_vec)
    #print(sum(keep)) # Better not be 0 or close to 0 (means very little data going into model)
    #'  Refit model with added independence claim using spp.latent posterior mean 
    #'  and specified predictor as y and x
    mod_out <- fit_one_dSep_claim(y = y_vec[keep], x = x_vec[keep], z = NULL, 
                                  n.chains = n.chains, n.adapt = n.adapt, n.burnin = n.burnin, 
                                  n.iter = n.iter, n.thin = n.thin, model_name, iter_num = i)
    
    #'  Flag convergence issues
    max_rhat <- suppressWarnings(max(unlist(mod_out$Rhat), na.rm = TRUE))
    if(is.finite(max_rhat) && max_rhat > 1.1) {
      warning(sprintf("Iteration %d: max Rhat = %.3f -- check convergence", i, max_rhat))
    }
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = mod_out, b_x = as.numeric(mod_out$sims.list$b_x), 
                 config = iter_deets, max_rhat = max_rhat), jags_out)
  }
  
  
  #'  ------------------------------------------------------------
  #####  Functions for ind. claims with only exogenous variables  #####
  #'  ------------------------------------------------------------
  #'  Function to grab covariate data and call fit_one_dSep_claim()
  fit_covariate_claim <- function(i, iterations, model_name, #y_array, x_array, iter_num,
                                  n.chains, n.adapt, n.burnin, n.iter, n.thin) {
    iter_deets <- iterations[[i]]
    y_vec <- as.vector(iter_deets$y_array)
    x_vec <- as.vector(iter_deets$x_array)
    keep <- !is.na(y_vec) & !is.na(x_vec)
    
    #'  Use simplified JAGS code and regression in fit_one_dSep_claim() to test claim
    mod_out <- fit_one_dSep_claim(y = y_vec[keep], x = x_vec[keep], z = NULL,
                                  n.chains = n.chains, n.adapt = n.adapt, n.burnin = n.burnin,
                                  n.iter = n.iter, n.thin = n.thin, model_name = model_name,
                                  iter_num = i)
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = mod_out, b_x = as.numeric(mod_out$sims.list$b_x), 
                 config = iter_deets), jags_out)
  }
  
  
  
  