  #'  -------------------------------------
  #'  Bayesian Structural Equation Models
  #'  Sarah Bassing
  #'  December 2025
  #'  -------------------------------------
  #'  Source data formatting script and run structural equation models (SEM) in
  #'  a Bayesian framework to test hypotheses about how predator-prey and predator- 
  #'  predator interactions influence wildlife populations in northern Idaho. This
  #'  formulation relies on original data structure (stacking 2020/2021 and 
  #'  2021/2022 for Yr1 --> Yr2 effect).
  #'  Species index order: 
  #'    1 = wolf; 2 = cougar; 3 = black bear; 4 = coyote; 5 = elk; 6 = moose; 7 = wtd
  #'  -------------------------------------
  
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
  head(density_wide_1YrLag_20s_22s)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  dat_final <- density_wide_1YrLag_20s_22s %>%
    mutate(obs = seq(1:nrow(.)),
           GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
    group_by(GMU_cluster) %>%
    mutate(uniqueCluster = cur_group_id()) %>%
    ungroup() #%>%
    # arrange(obs)
    
  
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
                    nWtd = nwtd, 
                    nharvest = nharv,
                    nforest = nfor,
                    nCluster = as.numeric(length(dat$obs)), #as.numeric(length(unique(dat$uniqueCluster))),
                    nSpp = 7,
                    #'  Standardize and set to numeric for each variable
                    wolf.t = as.numeric(scale(dat$wolf.T)),
                    wolf.tmin1 = as.numeric(scale(dat$wolf.Tminus1)),
                    lion.t = as.numeric(scale(dat$mountain_lion.T)),
                    lion.tmin1 = as.numeric(scale(dat$mountain_lion.Tminus1)),
                    bear.t = as.numeric(scale(dat$bear_black.T)), 
                    bear.tmin1 = as.numeric(scale(dat$bear_black.Tminus1)), 
                    coy.t = as.numeric(scale(dat$coyote.T)), 
                    coy.tmin1 = as.numeric(scale(dat$coyote.Tminus1)), 
                    elk.t = as.numeric(scale(dat$elk.T)), 
                    elk.tmin1 = as.numeric(scale(dat$elk.Tminus1)), 
                    moose.t = as.numeric(scale(dat$moose.T)), 
                    moose.tmin1 = as.numeric(scale(dat$moose.Tminus1)), 
                    wtd.t = as.numeric(scale(dat$whitetailed_deer.T)), 
                    wtd.tmin1 = as.numeric(scale(dat$whitetailed_deer.Tminus1)), 
                    harvest.t = as.numeric(scale(dat$annual_harvest.T)), 
                    harvest.tmin1 = as.numeric(scale(dat$annual_harvest.Tminus1)), 
                    forest.t = as.numeric(scale(dat$DisturbedForest_last20Yrs.T)), 
                    forest.tmin1 = as.numeric(scale(dat$DisturbedForest_last20Yrs.Tminus1)))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle <- bundle_dat(dat_final, nwolf = 7, nlion = 7, nbear = 4, ncoy = 2, 
                                 nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
                                 
  
  # save(data_JAGS_bundle, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_20s.RData")
  
  #'  Initial values
  initial_n <- function(spp) {
    ninit <- list(
      # sigma.spp <- as.numeric(runif(spp, -1, 1))
      wolf.t = runif(1, -0.5, 0.5),
      lion.t = runif(1, -0.5, 0.5),
      bear.t = runif(1, -0.5, 0.5),
      coy.t = runif(1, -0.5, 0.5),
      elk.t = runif(1, -0.5, 0.5),
      moose.t = runif(1, -0.5, 0.5),
      wtd.t = runif(1, -0.5, 0.5),
      harvest.t = runif(1, -0.5, 0.5)
      # forest.t = runif(1, -0.5, 0.5)
    )
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit <- initial_n()#spp = 7, nwolf = 7, nlion = 7, nbear = 4, ncoy = 2, nelk = 1, 
                     #nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  
  #'  Parameters monitored
  params <- c("wolf.t_1", "lion.t_1", "bear.t_1", "coy.t_1", "elk.t_1", "moose.t_1", 
              "wtd.t_1", "harvest.t_1", "forest.t_1", "sigma.spp") 
  #"wolf.t", "lion.t", "bear.t", "coy.t", "elk.t", "moose.t", "wtd.t", "harvest.t", "forest.t",
   
  
  #'  MCMC settings
  nc <- 3
  ni <- 500
  nb <- 100
  nt <- 1
  na <- 500
  
  
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter.R")
  start.time = Sys.time()
  inits_sem_tst <- function(){list(sppinits = ninit)}
  SEM_tst <- jags(data_JAGS_bundle, inits = NULL, params,
                      "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_tst$summary)
  which(SEM_tst$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_tst$samples)
  save(SEM_tst, file = paste0("./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_", Sys.Date(), ".RData"))
  
  