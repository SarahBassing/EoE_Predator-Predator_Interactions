  #'  --------------------------------------------
  #'  JAGS data setup functions for Bayesian SEM
  #'  July 2026
  #'  Sarah B. Bassing
  #'  --------------------------------------------
  #'  Set of functions to format and bundle input for JAGS to fit SEMs. Source
  #'  this script when fitting original SEMs and conducting d-Sep tests.
  #'  --------------------------------------------
  
  #'  ------------------------
  ####  Bundle data for JAGS  ####
  #'  ------------------------
  bundle_dat <- function(dat, covs, nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
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
                    nWSI = nwsi,
                    nCluster = as.numeric(length(unique(dat$wolf$time_t$cluster))),
                    nSpp = 7,
                    #'  Standardized posterior means and SD for each species and time step
                    wolf.t_hat = dat$wolf$time_t$posterior_mu_z,
                    wolf.t.sigma_hat = dat$wolf$time_t$posterior_sd_z,
                    wolf.tmin1_hat = dat$wolf$time_tmin1$posterior_mu_z,
                    wolf.tmin1.sigma_hat = dat$wolf$time_tmin1$posterior_sd_z,
                    lion.t_hat = dat$lion$time_t$posterior_mu_z,
                    lion.t.sigma_hat = dat$lion$time_t$posterior_sd_z,
                    lion.tmin1_hat = dat$lion$time_tmin1$posterior_mu_z,
                    lion.tmin1.sigma_hat = dat$lion$time_tmin1$posterior_sd_z,
                    bear.t_hat = dat$bear$time_t$posterior_mu_z,
                    bear.t.sigma_hat = dat$bear$time_t$posterior_sd_z,
                    bear.tmin1_hat = dat$bear$time_tmin1$posterior_mu_z,
                    bear.tmin1.sigma_hat = dat$bear$time_tmin1$posterior_sd_z,
                    coy.t_hat = dat$coy$time_t$posterior_mu_z,
                    coy.t.sigma_hat = dat$coy$time_t$posterior_sd_z, 
                    coy.tmin1_hat = dat$coy$time_tmin1$posterior_mu_z,
                    coy.tmin1.sigma_hat = dat$coy$time_tmin1$posterior_sd_z,
                    elk.t_hat = dat$elk$time_t$posterior_mu_z,
                    elk.t.sigma_hat = dat$elk$time_t$posterior_sd_z,
                    elk.tmin1_hat = dat$elk$time_tmin1$posterior_mu_z,
                    elk.tmin1.sigma_hat = dat$elk$time_tmin1$posterior_sd_z,
                    moose.t_hat = dat$moose$time_t$posterior_mu_z,
                    moose.t.sigma_hat = dat$moose$time_t$posterior_sd_z,
                    moose.tmin1_hat = dat$moose$time_tmin1$posterior_mu_z,
                    moose.tmin1.sigma_hat = dat$moose$time_tmin1$posterior_sd_z,
                    wtd.t_hat = dat$wtd$time_t$posterior_mu_z,
                    wtd.t.sigma_hat = dat$wtd$time_t$posterior_sd_z,
                    wtd.tmin1_hat = dat$wtd$time_tmin1$posterior_mu_z,
                    wtd.tmin1.sigma_hat = dat$wtd$time_tmin1$posterior_sd_z,
                    #'  Standardized harvest and habitat variables for 1 year time lag
                    wolfHarv.tmin1 = covs$wolf_harv$time_tmin1$wolf_harv_z,
                    lionHarv.tmin1 = covs$lion_harv$time_tmin1$lion_harv_z,
                    bearHarv.tmin1 = covs$bear_harv$time_tmin1$bear_harv_z,
                    elkHarv.tmin1 = covs$elk_harv$time_tmin1$elk_harv_z,
                    mooseHarv.tmin1 = covs$moose_harv$time_tmin1$moose_harv_z,
                    deerHarv.tmin1 = covs$deer_harv$time_tmin1$deer_harv_z,
                    wsi.tmin1 = covs$wsi$time_tmin1$wsi_z,
                    forest.tmin1 = covs$prop_disturbed$time_tmin1$prop_disturb_z,
                    road.tmin1 = covs$road_density$time_tmin1$road_density_z,
                    public.tmin1 = covs$public_land$time_tmin1$public_land_z)
    str(bundled)
    return(bundled)
  }
  
  #'  ------------------------------------------------------------
  ####  Generate initial values for each parameter (random node)  ####
  #'  ------------------------------------------------------------
  generate_inits <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
    
    #'  Generate random values for each species-specific beta (nwolf, nlion, etc.
    #'  based on number of species-specific betas to be estimated)
    list(
      wolf.t_1 = runif(nwolf, -0.5, 0.5),  # consider -1, 1
      lion.t_1 = runif(nlion, -0.5, 0.5),
      bear.t_1 = runif(nbear, -0.5, 0.5),
      coy.t_1 = runif(ncoy, -0.5, 0.5),
      elk.t_1 = runif(nelk, -0.5, 0.5),
      moose.t_1 = runif(nmoose, -0.5, 0.5),
      wtd.t_1 = runif(nwtd, -0.5, 0.5),
      wolfharvest.t_1 = runif(nharv, -0.5, 0.5),
      lionharvest.t_1 = runif(nharv, -0.5, 0.5),
      bearharvest.t_1 = runif(nharv, -0.5, 0.5),
      elkharvest.t_1 = runif(nharv, -0.5, 0.5),
      mooseharvest.t_1 = runif(nharv, -0.5, 0.5),
      deerharvest.t_1 = runif(nharv, -0.5, 0.5),
      wsi.t_1 = runif(nwsi, -0.5, 0.5),
      forest.t_1 = runif(nfor, -0.5, 0.5),
      road.t_1 = runif(nwsi, -0.5, 0.5),
      public.t_1 = runif(nwsi, -0.5, 0.5) # doesn't matter that using nwsi or other b/c all = 1
    )
  }
  