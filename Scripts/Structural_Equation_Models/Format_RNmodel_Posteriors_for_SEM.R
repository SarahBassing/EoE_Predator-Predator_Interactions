  #'  ----------------------------------------------------
  #'  Format Royle-Nichols model results for Bayesian SEM
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2026
  #'  ----------------------------------------------------
  #'  Load posterior distributions from RN models and format for as input for
  #'  Bayesian structural equation models
  #'  ----------------------------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(tidyverse)
  library(ggplot2)
  
  #'  Load JAGS outputs for species and cluster-specific posterior distributions
  file_names <- list.files(path = "./Outputs/Relative_Abundance/RN_model/JAGS_out/", 
                           pattern = "\\.RData$", full.names = T)
  for(i in 1:length(file_names)) {
    load(file_names[i])
  }
  
  #'  List model outputs per year
  mods_yr1 <- list(RN_wolf_20s, RN_lion_20s, RN_bear_20s, RN_coy_20s, RN_elk_20s, RN_moose_20s, RN_wtd_20s)
  mods_yr2 <- list(RN_wolf_21s, RN_lion_21s, RN_bear_21s, RN_coy_21s, RN_elk_21s, RN_moose_21s, RN_wtd_21s)
  mods_yr3 <- list(RN_wolf_22s, RN_lion_22s, RN_bear_22s, RN_coy_22s, RN_elk_22s, RN_moose_22s, RN_wtd_22s)
  
  #'  Function to save estimated posterior distributions of interest for each species,
  #'  cluster, and year
  cluster_mcmc_matrix <- function(mod_post, cluster_out) {
    #'  Empty list to hold model outputs of interest
    mod_out <- NULL
    #'  Loop through posteriors and save outputs of interest
    for(i in 1:length(cluster_out)) {
      mod_out[[i]] <- mod_post$sims.list[names(mod_post$sims.list) %in% cluster_out[i]]
    }
    
    #'  Create an empty matrix with dimensions specific to clusters & iterations in model
    cluster_posteriors <- matrix(data = NA, nrow = length(cluster_out), ncol = length(mod_post$sims.list[[1]]), byrow = TRUE)
    #'  Loop through matrix and fill with estimated posterior distribution for each RDI cluster
    for(i in 1:nrow(cluster_posteriors)) {
      cluster_posteriors[i,] <- mod_out[[i]][[1]]
    }
  
    #' #'  Turn matrix into a data frame with uniqueCluster as first column
    #' uniqueCluster <- as.data.frame(seq(1:nrow(cluster_posteriors)))
    #' names(uniqueCluster) <- "uniqueCluster"
    #' cluster_posteriors <- cbind(uniqueCluster, cluster_posteriors) 
    
    ####  THIS NEEDS TO BE UPDATED WITH BETTER CLUSTER INFO
    ####  ADD IN COVARIATE DATA HERE TOO!!!!
    
    return(cluster_posteriors)
  }
  #'  Names of estimated posteriors of interest per year
  rdi.cl_2020 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9")
  rdi.cl_2021.2022 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10", 
                        "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20", 
                        "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  #'  Run function for each year and species
  posteriors_20s <- lapply(mods_yr1, cluster_mcmc_matrix, cluster_out = rdi.cl_2020)
  posteriors_21s <- lapply(mods_yr2, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  posteriors_22s <- lapply(mods_yr3, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  
  #'  Name lists by species
  spp <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
  names(posteriors_20s) <- spp
  names(posteriors_21s) <- spp
  names(posteriors_22s) <- spp
  
  #'  Double check these look right
  #'  First 10 iterations for each cluster in matrix vs model output
  print(posteriors_20s[[1]][, 1:10]) # wolf (2020)
  RN_wolf_20s$sims.list$rdi.cl1[1:10]; RN_wolf_20s$sims.list$rdi.cl9[1:10]
  print(posteriors_21s[[3]][, 1:10]) # bear (2021)
  RN_bear_21s$sims.list$rdi.cl1[1:10]; RN_bear_21s$sims.list$rdi.cl9[1:10]
  print(posteriors_22s[[7]][, 1:10]) # white-tailed deer (2022)
  RN_wtd_22s$sims.list$rdi.cl1[1:10]; RN_wtd_22s$sims.list$rdi.cl9[1:10]
  
  #'  Stack species-specific estimated posteriors across iterations 
  #'  (e.g., time t = wolf 2021 and 2022 vs. time t-1 = wolf 2020 and 2021)
  stacked_posteriors <- function(yr1, yr2, yr3) {
    #'  Posteriors from time t 
    time_t <- rbind(yr2, yr3)
    #'  Posteriors from time t-1, given an annual lag from time t
    time_tmin1 <- rbind(yr1, yr2)
    
    #'  List posteriors from time t and t-1
    timelag_list <- list(time_t, time_tmin1)
    return(timelag_list)
  }
  #'  Run function for each species
  wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]])
  lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]])
  bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]])
  coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]])
  elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]])
  moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]])
  wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]])
  
  #'  Double check everything looks right
  print(wolf_timelag[[1]][,1:10])
  print(wolf_timelag[[2]][,1:10])
  print(elk_timelag[[1]][,1:10])
  print(elk_timelag[[2]][,1:10])
  
  
  
  
  
  