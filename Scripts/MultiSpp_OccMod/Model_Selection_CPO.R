  #'  ---------------------------------
  #'  Calculate CPO for model seleciton
  #'  ---------------------------------
  #'  Attempt at calculating CPO based on code written by Mason Fidino
  #'  https://github.com/mfidino/dcom/blob/master/calculate_cpo.R
  #'  ---------------------------------
  
  library(coda)
  library(stringr)
  
  #'  Load bundled data original provided to JAGS and JAGS output
  load("./Data/MultiSpp_OccMod_Outputs/bundled_predator_data_list.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(setup_habitat_rx)_px(setup_effort).RData")

  #'  Reformat JAGS output to a model matrix
  #'  Columns = different parameters; rows = posterior simulations
  # mm <- as.matrix(coy.bob.hab$sims.list, chains = TRUE)
  mm <- as.matrix(coy.bob.hab$samples, chains = TRUE)
  if(!is.matrix(mm)) {
    stop("The object mm must be a matrix")
  }
  colnames(mm) <- str_replace(colnames(mm), "Spp12", "Both")
  
  #'  Assign all objects in data list so they don't have to be indexed constantly
  bundled_coy.bob <- bundled_pred_list[[6]]
  for(i in 1:length(bundled_coy.bob)) {
    assign(names(bundled_coy.bob)[i], bundled_coy.bob[[i]])
  }
  
  #'  Assign number of covariates on each submodel and number of species
  nspec <- 2
  ncov_psi <- 4
  ncov_rho <- 3
  
  #'  Pull out just covaraites used in the model
  psi_cov2 <- psi_cov[,1:4]
  rho_cov2 <- rho_cov[,,c(1,3,5)]
  
  #'  Create likelihood matrix
  #'  Stores likelihood of each observation based on parameters in model for each
  #'  step of MCMC chain
  lik_mat <- matrix(0, nrow = nrow(mm), ncol = nsites)
  #'  Loop through each step of the MCMC chain
  #'  ncol = number of covariates (including intercept), nrow = number of species
  for(o in 1:nrow(mm)) {
    #'  First order psi
    fo.psi <- matrix(mm[o, grep("betaSpp", colnames(mm))], ncol = ncov_psi, nrow = 2, byrow = T)
    #'  Second order psi interaction
    so.psi <- matrix(mm[o, grep("betaBoth", colnames(mm))], ncol = ncov_psi, nrow = 2, byrow = T)
    if(sum(is.na(so.psi)>0)) so.psi[is.na(so.psi)] <- 0 # make 0 if not in the model
    #'  First order rho
    fo.rho <- matrix(mm[o, grep("alphaSpp", colnames(mm))], ncol = ncov_rho, nrow = 2, byrow = T)
    #' #'  Second order rho interaction
    #' so.rho <- matrix(mm[o, grep("alphaBoth", colnames(mm)], ncol = ncov_rho, nrow = 2, byrow = T)
    #'  Estimated community state at site i
    z <- matrix(mm[o, grep("z", colnames(mm))], ncol = 1, nrow = nsites, byrow = T)
  
    #'  Occupancy linear predictor
    psi <- fo.psi %*% t(psi_cov2)
    #'  Detection linear predictor
    rho <- fo.rho %*% t(rho_cov2)
    
    #'  Set up arrays for the species interactions (psi[species of interest, species conditioning on])
    # psix <- array(0, dim = c(nspec, nspec, nsite))
    psix <- matrix(0, dim = c(nspec, nsite))         # <-- NOT SURE IF THIS IS RIGHT...
    for(k in 1:n_inxs) {  #n_inxs might actually just be 1 since it's only a 2 spp model
      # inxs on occupancy 
      # psix[rows_vec[k], cols_vec[k], ] <-  psi_cov2[j, ] %*% so.psi[k, ]
      psix[rows_vec[k], ] <-  psi_cov2[j, ] %*% so.psi[k, ]   # <--- I'm sure I'm messing up the math here
    }
    
    #'  Numerator of the softmax function for each of the 4 community states a site can be in.
    fsm <- matrix(0, ncol = ncat, nrow = nsites)
    fsm[, 1] <- 1 # -------------------------------------------------------- U
    fsm[, 2] <- exp(psi[1, ]) # -------------------------------------------- Spp1
    fsm[, 3] <- exp(psi[2, ]) # -------------------------------------------- Spp2
    fsm[, 4] <- exp(psix[2, ]) # ------------------------------------------- Both <-- this indexing is wrong
    
    
    
    
    
    }
  
  
  
  
  
  

  
  
  
  