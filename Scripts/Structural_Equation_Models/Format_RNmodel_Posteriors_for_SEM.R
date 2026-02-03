  #'  ----------------------------------------------------
  #'  Format Royle-Nichols model results for Bayesian SEM
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2026
  #'  ----------------------------------------------------
  #'  Load posterior distributions from RN models and format for as input for
  #'  Bayesian structural equation models
  #'  ----------------------------------------------------
  
  #' #'  Clean workspace
  #' rm(list = ls())
  
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
  
  
  #'  Function to save posterior means and SD for each site per species, cluster, and year
  cluster_posterior_summary <- function(mod_post, cluster_out) {
    #'  Empty lists to hold model outputs of interest
    posterior_mu <- posterior_sd <- c()
    #'  Loop through posteriors and save outputs of interest that match cluster names
    for(i in 1:length(cluster_out)) {
      posterior_mu[i] <- unlist(mod_post$mean[names(mod_post$mean) %in% cluster_out[i]])
      posterior_sd[i] <- unlist(mod_post$sd[names(mod_post$sd) %in% cluster_out[i]])
    }
    
    #'  Bind mu and sd outputs into a matrix
    cluster_posteriors <- cbind(posterior_mu, posterior_sd)
    return(cluster_posteriors)
  }
  #'  Names of estimated posteriors of interest per year
  rdi.cl_2020 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9",
                   "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  rdi.cl_2021.2022 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10", 
                        "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20", 
                        "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  
  #'  Run function for each year and species
  posteriors_20s <- lapply(mods_yr1, cluster_posterior_summary, cluster_out = rdi.cl_2020)
  posteriors_21s <- lapply(mods_yr2, cluster_posterior_summary, cluster_out = rdi.cl_2021.2022)
  posteriors_22s <- lapply(mods_yr3, cluster_posterior_summary, cluster_out = rdi.cl_2021.2022)
  
  #'  Name lists by species
  spp <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
  names(posteriors_20s) <- spp
  names(posteriors_21s) <- spp
  names(posteriors_22s) <- spp
  
  #'  Double check these look right
  #'  Mean and SD matrix vs cluster specific summary outputs
  print(posteriors_20s[[1]]) # wolf (2020)
  RN_wolf_20s$mean[9:10]; RN_wolf_20s$sd[9:10]  # 1st two clusters
  print(posteriors_21s[[3]]) # bear (2021)
  RN_bear_21s$mean[9:10]; RN_bear_21s$sd[9:10]  # 1st two clusters
  print(posteriors_22s[[7]]) # white-tailed deer (2022)
  RN_wtd_22s$mean[31:32]; RN_wtd_22s$sd[31:32]  # last two clusters
  
  #'  Create new posteriors for 2021 that exclude GMU 1 RDI posteriors 
  #'  Needed below when stacking data for time t vs t-1
  #'  posteriors_21s[[i]][10:20,] needs to go (corresponds to GMU 1 clusters)
  drop_gmu1 <- function(post_rdi) {
    #'  Snag posteriors from GMU 10A and GMU 6
    gmu10a_rdi <- post_rdi[1:9,]
    gmu6_rdi <- post_rdi[21:24,]
    
    #'  Bind into a single matrix that no long includes posteriors from GMU 1
    rdi_2021_no_gmu1 <- rbind(gmu10a_rdi, gmu6_rdi)
    return(rdi_2021_no_gmu1)
  }
  posteriors_21s_noGMU1 <- lapply(posteriors_21s, drop_gmu1)
  
  #'  Stack species-specific estimated posterior summaries across iterations 
  #'  (e.g., time t = wolf 2021 and 2022 vs. time t-1 = wolf 2020 and 2021)
  stacked_posteriors <- function(yr1, yr2, yr3, yr2_noGMU1) {
    #'  Posteriors from time t 
    #'  NOTE: yr2 for time_t excludes GMU 1 posteriors because there were no
    #'  posteriors for GMU 1 in the corresponding time_tmin1 (yr1, 2020)
    time_t <- rbind(yr2_noGMU1, yr3)
    
    #'  Posteriors from time t-1, given an annual lag from time 
    #'  NOTE: yr2 for time_tmin1 includes GMU 1 posteriors because the corresponding 
    #'  yr3 posteriors in time_t also include GMU 1
    time_tmin1 <- rbind(yr1, yr2)
    
    #'  List posteriors from time t and t-1
    timelag_list <- list(time_t, time_tmin1)
    return(timelag_list)
  }
  #'  Run function for each species
  wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]], posteriors_21s_noGMU1[[1]])
  lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]], posteriors_21s_noGMU1[[2]])
  bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]], posteriors_21s_noGMU1[[3]])
  coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]], posteriors_21s_noGMU1[[4]])
  elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]], posteriors_21s_noGMU1[[5]])
  moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]], posteriors_21s_noGMU1[[6]])
  wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]], posteriors_21s_noGMU1[[7]])
  
  #'  Double check everything looks right (both time steps should have 37 rows)
  #'  [[1]] = time t (yr2 and yr3); [[2]] = time t-1 (yr1 and yr2)
  print(wolf_timelag[[1]]) # observations [1:13,] in [[1]] should be observations [c(14:22, 34:37)] in [[2]]
  print(wolf_timelag[[2]]) 
  print(elk_timelag[[1]])
  print(elk_timelag[[2]])
  
  
  
  
  #' ####  Save full posteriors (all iterations)  ####
  #' ####  Not necessary for how uncertainty from RN models are being propagated in SEM  ####
  #' #'  Function to save estimated posterior distributions of interest for each species,
  #' #'  cluster, and year
  #' cluster_mcmc_matrix <- function(mod_post, cluster_out) {
  #'   #'  Empty list to hold model outputs of interest
  #'   mod_out <- NULL
  #'   #'  Loop through posteriors and save outputs of interest
  #'   for(i in 1:length(cluster_out)) {
  #'     mod_out[[i]] <- mod_post$sims.list[names(mod_post$sims.list) %in% cluster_out[i]]
  #'   }
  #'   
  #'   #'  Create an empty matrix with dimensions specific to clusters & iterations in model
  #'   cluster_posteriors <- matrix(data = NA, nrow = length(cluster_out), ncol = length(mod_post$sims.list[[1]]), byrow = TRUE)
  #'   #'  Loop through matrix and fill with estimated posterior distribution for each RDI cluster
  #'   for(i in 1:nrow(cluster_posteriors)) {
  #'     cluster_posteriors[i,] <- mod_out[[i]][[1]]
  #'   }
  #' 
  #'   #' #'  Turn matrix into a data frame with uniqueCluster as first column
  #'   #' uniqueCluster <- as.data.frame(seq(1:nrow(cluster_posteriors)))
  #'   #' names(uniqueCluster) <- "uniqueCluster"
  #'   #' cluster_posteriors <- cbind(uniqueCluster, cluster_posteriors) 
  #'   
  #'   ####  ADD IN COVARIATE DATA HERE!!!!
  #'   
  #'   return(cluster_posteriors)
  #' }
  #' #'  Names of estimated posteriors of interest per year
  #' rdi.cl_2020 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9",
  #'                  "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  #' rdi.cl_2021.2022 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10", 
  #'                       "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20", 
  #'                       "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  #' #'  Run function for each year and species
  #' posteriors_20s <- lapply(mods_yr1, cluster_mcmc_matrix, cluster_out = rdi.cl_2020)
  #' posteriors_21s <- lapply(mods_yr2, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  #' posteriors_22s <- lapply(mods_yr3, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  #' 
  #' #'  Name lists by species
  #' spp <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
  #' names(posteriors_20s) <- spp
  #' names(posteriors_21s) <- spp
  #' names(posteriors_22s) <- spp
  #' 
  #' #'  Double check these look right
  #' #'  First 10 iterations for each cluster in matrix vs model output
  #' print(posteriors_20s[[1]][, 1:10]) # wolf (2020)
  #' RN_wolf_20s$sims.list$rdi.cl1[1:10]; RN_wolf_20s$sims.list$rdi.cl24[1:10] # note that rdi.cl10 - rdi.cl20 are not in 2020 posteriors
  #' print(posteriors_21s[[3]][, 1:10]) # bear (2021)
  #' RN_bear_21s$sims.list$rdi.cl1[1:10]; RN_bear_21s$sims.list$rdi.cl24[1:10]
  #' print(posteriors_22s[[7]][, 1:10]) # white-tailed deer (2022)
  #' RN_wtd_22s$sims.list$rdi.cl1[1:10]; RN_wtd_22s$sims.list$rdi.cl24[1:10]
  #' 
  #' #'  Create new posteriors for 2021 that exclude GMU 1 RDI posteriors 
  #' #'  Needed below when stacking data for time t vs t-1
  #' #'  posteriors_21s[[i]][10:20,] needs to go (corresponds to GMU 1 clusters)
  #' drop_gmu1 <- function(post_rdi) {
  #'   #'  Snag posteriors from GMU 10A and GMU 6
  #'   gmu10a_rdi <- post_rdi[1:9,]
  #'   gmu6_rdi <- post_rdi[21:24,]
  #'   
  #'   #'  Bind into a single matrix that no long includes posteriors from GMU 1
  #'   rdi_2021_no_gmu1 <- rbind(gmu10a_rdi, gmu6_rdi)
  #'   return(rdi_2021_no_gmu1)
  #' }
  #' posteriors_21s_noGMU1 <- lapply(posteriors_21s, drop_gmu1)
  #' 
  #' #'  Stack species-specific estimated posteriors across iterations 
  #' #'  (e.g., time t = wolf 2021 and 2022 vs. time t-1 = wolf 2020 and 2021)
  #' stacked_posteriors <- function(yr1, yr2, yr3, yr2_noGMU1) {
  #'   #'  Posteriors from time t 
  #'   #'  NOTE: yr2 for time_t excludes GMU 1 posteriors because there were no
  #'   #'  posteriors for GMU 1 in the corresponding time_tmin1 (yr1, 2020)
  #'   time_t <- rbind(yr2_noGMU1, yr3)
  #'   
  #'   #'  Posteriors from time t-1, given an annual lag from time 
  #'   #'  NOTE: yr2 for time_tmin1 includes GMU 1 posteriors because the corresponding 
  #'   #'  yr3 posteriors in time_t also include GMU 1
  #'   time_tmin1 <- rbind(yr1, yr2)
  #'   
  #'   #'  List posteriors from time t and t-1
  #'   timelag_list <- list(time_t, time_tmin1)
  #'   return(timelag_list)
  #' }
  #' #'  Run function for each species
  #' wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]], posteriors_21s_noGMU1[[1]])
  #' lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]], posteriors_21s_noGMU1[[2]])
  #' bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]], posteriors_21s_noGMU1[[3]])
  #' coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]], posteriors_21s_noGMU1[[4]])
  #' elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]], posteriors_21s_noGMU1[[5]])
  #' moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]], posteriors_21s_noGMU1[[6]])
  #' wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]], posteriors_21s_noGMU1[[7]])
  #' 
  #' #'  Double check everything looks right (both time steps should have 37 rows)
  #' #'  [[1]] = time t (yr2 and yr3); [[2]] = time t-1 (yr1 and yr2)
  #' print(wolf_timelag[[1]][,1:10])
  #' print(wolf_timelag[[2]][,1:10])
  #' print(elk_timelag[[1]][,1:10])
  #' print(elk_timelag[[2]][,1:10])
  
  

  