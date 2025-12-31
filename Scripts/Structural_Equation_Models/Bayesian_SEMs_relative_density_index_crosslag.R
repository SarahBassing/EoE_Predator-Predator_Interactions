  #'  ----------------------------------------------------
  #'  Bayesian Structural Equation Models with cross lags
  #'  Sarah Bassing
  #'  December 2025
  #'  ----------------------------------------------------
  #'  Source data formatting script and run structural equation models (SEM) in
  #'  a Bayesian framework to test hypotheses about how predator-prey and predator- 
  #'  predator interactions influence wildlife populations in northern Idaho. This
  #'  formulation includes a 1-year time lag where 2020 affects 2021, and 2021 
  #'  affects 2022. Species index order: 
  #'    1 = wolf; 2 = cougar; 3 = black bear; 4 = coyote; 5 = elk; 6 = moose; 7 = wtd
  #'  ----------------------------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Run script that formats covariate data
  source("./Scripts/Structural_Equation_Models/Format_covariate_data_for_SEMs.R") 
  
  #'  Run script that formats density data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  
  #'  Take a quick look
  head(density_long_annual_combo)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  #'  Generate unique ID for each cluster across all GMUs and rearrange data frame
  dat_final <- density_long_annual_combo %>%
    mutate(GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
    group_by(GMU_cluster) %>%
    mutate(uniqueCluster = cur_group_id()) %>%
    ungroup() %>%
    relocate(uniqueCluster, .after = ClusterID) %>%
    dplyr::select(-GMU_cluster) %>%
    arrange(uniqueCluster, timestep) %>%
    mutate_all(~replace(., is.na(.), 0)) # filling in NAs just to see if this helps with compilation error
  
  c <- length(unique(dat_final$uniqueCluster))
  t <- length(unique(dat_final$timestep))
  
  wolf.matrix <- matrix(dat_final$wolf, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  lion.matrix <- matrix(dat_final$mountain_lion, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  bear.matrix <- matrix(dat_final$bear_black, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  coy.matrix <- matrix(dat_final$coyote, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  elk.matrix <- matrix(dat_final$elk, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  moose.matrix <- matrix(dat_final$moose, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  deer.matrix <- matrix(dat_final$whitetailed_deer, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))
  harv.matrix <- matrix(dat_final$annual_harvest, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))   ####  LOOK INTO GETTNG GMU1 YR1 harvest data
  for.matrix <- matrix(dat_final$DisturbedForest_last20Yrs, nrow = c, ncol = t, byrow = TRUE, 
                        dimnames = list(NULL, c("t1", "t2", "t3")))   ####  LOOK INTO GETTING GMU1 YR1 forest data
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle data for JAGS
  bundle_dat <- function(dat, nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor) {
    #'  Bundle data for JAGS
    bundled <- list(nWolf = nwolf,
                    nLion = nlion, 
                    nBear = nbear, 
                    nCoy = ncoy, 
                    nElk = nelk,
                    nMoose = nmoose, 
                    nDeer = nwtd, 
                    nharvest = nharv,
                    nforest = nfor,
                    nCluster = as.numeric(length(unique(dat$uniqueCluster))),
                    nTimestep = as.numeric(length(unique(dat$timestep))),
                    nSpp = 7,
                    #'  Standardize and set to numeric for each variable
                    wolf = as.matrix(wolf.matrix), #as.numeric(scale(dat$wolf)),  # THINK ABOUT SCALING!!!
                    lion = as.matrix(lion.matrix), #as.numeric(scale(dat$mountain_lion)),
                    bear = as.matrix(bear.matrix), #as.numeric(scale(dat$bear_black)), 
                    coy = as.matrix(coy.matrix), #as.numeric(scale(dat$coyote)), 
                    elk = as.matrix(elk.matrix), #as.numeric(scale(dat$elk)), 
                    moose = as.matrix(moose.matrix), #as.numeric(scale(dat$moose)), 
                    wtd = as.matrix(deer.matrix), #as.numeric(scale(dat$whitetailed_deer)), 
                    harv = as.matrix(harv.matrix), #as.numeric(scale(dat$annual_harvest)), 
                    forest = as.matrix(for.matrix))#as.numeric(scale(dat$DisturbedForest_last20Yrs)))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle <- bundle_dat(dat_final, nwolf = 3, nlion = 1, nbear = 2, ncoy = 0, 
                                 nelk = 0, nmoose = 0, nwtd = 0, nharv = 1, nfor = 0)
                                 #nwolf = 7, nlion = 7, nbear = 4, ncoy = 2, 
                                 #nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  
  
  # save(data_JAGS_bundle, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_20s.RData")
  
  #'  Initial values
  initial_n <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv) { #spp #nfor
    ninit <- list(
      # sigma.spp <- as.numeric(runif(spp, -1, 1))
      # beta.wolf = runif(nwolf, -0.5, 0.5),
      # beta.lion = runif(nlion, -0.5, 0.5),
      # beta.bear = runif(nbear, -0.5, 0.5),
      # beta.coy = runif(ncoy, -0.5, 0.5),
      # beta.elk = runif(nelk, -0.5, 0.5),
      # beta.moose = runif(nmoose, -0.5, 0.5),
      # beta.wtd = runif(nwtd, -0.5, 0.5),
      # beta.harvest = runif(nharv, -0.5, 0.5)
      # beta.forest = runif(nfor, -0.5, 0.5)
    )
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit <- initial_n()#spp = 7, nwolf = 7, nlion = 7, nbear = 4, ncoy = 2, nelk = 1, 
  #nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  
  #'  Parameters monitored
  params <- c("beta.wolf", "beta.lion", "beta.bear") #"beta.coy", "beta.elk", "beta.moose", 
              #"beta.wtd", "beta.harvest", "beta.forest", "sigma.spp", "sigma.cluster") 
  
  
  #'  MCMC settings
  nc <- 3
  ni <- 500
  nb <- 100
  nt <- 1
  na <- 500
  
  
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter_crosslag.R")
  start.time = Sys.time()
  inits_sem_tst <- function(){list(sppinits = ninit)}
  SEM_tst <- jags(data_JAGS_bundle, inits = NULL, params,
                  "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_crosslag.txt",
                  n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                  n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_tst$summary)
  which(SEM_tst$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_tst$samples)
  save(SEM_tst, file = paste0("./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_crosslag_", Sys.Date(), ".RData"))
  
  