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
  # 
  # library(piecewiseSEM)
  # library(semEff)
  # library(labelled)
  # library(DiagrammeR)
  # library(lme4)
  library(jagsUI)
  library(tidyverse)
  
  #'  Run script that formats covariate data
  source("./Scripts/Structural_Equation_Models/Format_covariate_data_for_SEMs.R") 
  
  #'  Run script that formats density data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  
  #'  Take a quick look
  head(density_wide_1YrLag_20s_22s)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle data for JAGS
  bundle_dat <- function(dh, nsite, nsurvey, cov, effort) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #'  Count number of sites per GMU
    ncams_perGMU <- cov %>%
      group_by(GMU) %>%
      summarise(nsites = n()) %>%
      ungroup()
    #'  Split up covariates by GMU
    covs_GMU10A <- filter(cov, GMU == 1)
    covs_GMU6 <- filter(cov, GMU == 2)
    covs_GMU1 <- filter(cov, GMU == 3)
    #'  Bundle data for JAGS
    bundled <- list(y = dh, 
                    nsites = dim(dh)[1], 
                    nsurveys = dim(dh)[2], 
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    ncams1 = as.numeric(ncams_perGMU[1,2]), # GMU10A
                    ncams2 = as.numeric(ncams_perGMU[2,2]), # GMU6
                    ncams3 = as.numeric(ifelse(is.na(ncams_perGMU[3,2]), 0, ncams_perGMU[3,2])), #GMU1
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup),
                    forest = as.numeric(cov$PercFor), 
                    elev = as.numeric(cov$Elev), 
                    ndays = as.numeric(cov$nDays), 
                    seffort = as.numeric(effort),
                    #'  GMU-specific covariates for predicting N
                    forest1 = as.numeric(covs_GMU10A$PercFor),
                    forest2 = as.numeric(covs_GMU6$PercFor),
                    forest3 = as.numeric(covs_GMU1$PercFor),
                    elev1 = as.numeric(covs_GMU10A$Elev),
                    elev2 = as.numeric(covs_GMU6$Elev),
                    elev3 = as.numeric(covs_GMU1$Elev),
                    #'  Area of each (km2)
                    area1 = as.numeric(8527.31),
                    area2 = as.numeric(5905.44),
                    area3 = as.numeric(14648.92))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_20s <- lapply(DH_npp20s_RNmod, bundle_dat, cov = stations_npp20s, effort = zeffort_RNmod[[1]])
  
  save(data_JAGS_bundle_20s, file = "./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_20s.RData")
  
  #'  Initial values
  #'  Using naive occupancy as a starting point for local abundance
  initial_n <- function(dh) {
    #'  Max value per row
    ninit <- apply(dh, 1, max, na.rm = TRUE)
    ninit <- as.vector(ninit)
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit <- lapply(DH_npp20s_RNmod, initial_n)
  
  #'  Parameters monitored
  params <- c("harvest.t_1", "wolf.t_1", "lion.t_1", "bear.t_1", "coy.t_1", "elk.t_1", "moose.t_1", "wtd.t_1", "sigma.spp")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept, i.e., mean lambda for GMU10A 
  #'  mean.r = the intercept, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 500
  nb <- 100
  nt <- 1
  na <- 500
  
  
  
  #'  JAGS model (top-down, interference model)
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      harvest.t_1 ~ dnorm(0, 0.01)
      
      #'  Priors for species lag effects
      for(w in 1:nWolf) {
        wolf.t_1[w] ~ dnorm(0, 0.01)  # look more into precision value...
      }
      for(l in 1:nLion) {
        lion.t_1[l] ~ dnorm(0, 0.01)
      }
      for(b in 1:nBear) {
        bear.t_1[b] ~ dnorm(0, 0.01)
      }
      for(c in 1:nCoy) {
        coy.t_1[c] ~ dnorm(0, 0.01)
      }
      for(e in 1:nElk) {
        elk.t_1[e] ~ dnorm(0, 0.01)
      }
      for(m in 1:nMoose) {
        moose.t_1[m] ~ dnorm(0, 0.01)
      }
      for(wtd in 1:nWtd) {
        wtd.t_1[wtd] ~ dnorm(0, 0.01)
      }
      
      #'  SD prior for each regression
      for(k in 1:nSpp) {
        sigma.spp[7] ~ dunif(0, 10)
        tau.spp[7] <- 1 / pow(sigma.spp[7], 2)
      }
      
      
      #'  Likelihood
      #'  ----------
      for(i in 1:nCluster) {
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        mu.wolf.t[i] <- wolf.t_1[1, i] + harvest.t_1[i]  #  May need to fix indexing because of wolf.t_1
        
        lion.t ~ dnorm(mu.lion.t[i], tau.spp[2])
        mu.lion.t <- lion.t_1[i] + wolf.t_1[i] + bear.t_1[i]
        
        bear.t ~ dnorm(mu.bear.t[i], tau.spp[3])
        mu.bear.t <- bear.t_1[i] + wolf.t_1[i]
        
        coy.t ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t <- coy.t_1[i] + wolf.t_1[i] + lion.t_1[i]
        
        elk.t ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t <- elk.t_1[i] + wolf.t_1[i] + lion.t_1[i] + bear.t_1[i]
        
        moose.t ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t <- moose.t_1[i] + wolf.t_1[i]
        
        wtd.t ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t <- wtd.t_1[i] + wolf.t_1[i] + lion.t_1[i] + bear.t_1[i] + coy.t_1[i]
      
      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")
  
  
  top_down_inter <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter)
  