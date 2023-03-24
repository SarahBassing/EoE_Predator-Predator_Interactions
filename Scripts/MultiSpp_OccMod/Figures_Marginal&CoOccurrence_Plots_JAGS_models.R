  #'  ------------------------------
  #'  Plot co-occurrence results
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------
  #'  Script to create figures of predicted co-occurrence and covariate effects
  #'  from multispecies occupancy models run in JAGS.
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
    
  #'  Load libraries
  library(ggplot2)
  library(tidyverse)
  
  #'  Load covariate data
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s.RData")
  
  #'  Identify top models
  load("./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.RData")
  print(topmodels)
  
  #'  -------------------------------------------------
  ####  Predict Pr(occupancy) across covariate values  ####
  #'  -------------------------------------------------
  #'  MASSIVE function the predict marginal and conditional occupancy for each species
  predict_occupancy <- function(mod, ncat, npoints, focal_cov, psi_cov, psi_inxs_cov, psi_cov_index, psi_inxs_cov_index) {
    #'  Rename model and review
    out1 <- mod
    # print(out1$summary[1:50,])
    #'  Number of states
    ncat <- ncat
    #'  Number of draws in the MCMC chain
    ndraws <- out1$mcmc.info$n.samples
    #'  Number of observations to predict
    npoints <- npoints
    
    #'  Create range of covariate values to predict across
    print(r <- range(focal_cov))
    print(mean.cov <- mean(focal_cov))
    print(sd.cov <- sd(focal_cov))
    cov.pred.orig <- seq(r[1], r[2], length.out = npoints)
    #'  Scale new data to be consistent with data used in model
    cov.pred <- (cov.pred.orig - mean.cov) / sd.cov
    
    #'  Create model matrices for prediction covariates on psi and psix
    #'  Fill matrix with categorical or mean value of each covariate
    psi_cov <- psi_cov
    psi_covs <- matrix(psi_cov, nrow = npoints, ncol = length(psi_cov), byrow = TRUE)
    #'  Same thing but for covariates on interaction
    psi_inxs_cov <- psi_inxs_cov
    psi_inxs_covs <- matrix(psi_inxs_cov, nrow = npoints, ncol = length(psi_inxs_cov), byrow = TRUE)
    #'  Replace focal covariate column with new scaled data
    #'  New data must be in correct column based on order of covariates in original model
    psi_covs[,psi_cov_index] <- cov.pred
    print(head(psi_covs))
    psi_inxs_covs[,psi_inxs_cov_index] <- cov.pred
    print(head(psi_inxs_covs))
    
    #'  Assemble linear predictors across all iterations of the MCMC chains
    psiSpp1 <- out1$sims.list$betaSpp1 %*% t(psi_covs)
    psiSpp2 <- out1$sims.list$betaSpp2 %*% t(psi_covs)
    psiSpp12 <- psiSpp1 + psiSpp2 + out1$sims.list$betaSpp12 %*% t(psi_inxs_covs) 
    
    #'  Compute latent state vector (lsv) based on new data and means of other covariates
    lsv <- array(NA, dim = c(ndraws, npoints, ncat))
    lsv[,,1] <- 1
    lsv[,,2] <- exp(psiSpp1)
    lsv[,,3] <- exp(psiSpp2)
    lsv[,,4] <- exp(psiSpp12)
    
    #'  Compute latent state probability of each state from MCMC chains (lsp)
    #'  denominator  = sum of lsv for each draw x point
    denom <- apply(lsv, 1:2, sum)
    lspMC <- sweep(lsv, 1:2, denom, "/")
    dimnames(lspMC)[[3]] <- c("U", "Spp1", "Spp2", "Both")
    str(lspMC)
    
    #'  Fill in summary values for 'lsp'
    lsp <- array(NA, dim = c(npoints, ncat, 4))
    dimnames(lsp) <- list(NULL, c("U", "Spp1", "Spp2", "Both"),
                          c("mean", "SD", "lower", "upper"))
    lsp[,,1] <- apply(lspMC, 2:3, mean)
    lsp[,,2] <- apply(lspMC, 2:3, sd)
    lsp[,,3] <- apply(lspMC, 2:3, quantile, probs = 0.025)
    lsp[,,4] <- apply(lspMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute and plot marginal probabilities
    #'  --------------------------------------- 
    #'  Occupancy probability for each species, regardless of presence or
    #'  absence of other species
    #'  Get MCMC chains
    marginalMC <- array(NA, dim = c(ndraws, npoints, 2))
    dimnames(marginalMC)[[3]] <- c("Spp1", "Spp2")
    marginalMC[,,"Spp1"] <- apply(lspMC[,,c("Spp1", "Both")], 1:2, sum)
    marginalMC[,,"Spp2"] <- apply(lspMC[,,c("Spp2", "Both")], 1:2, sum)
    
    #'  Fill in summary values for marginal probabilities
    marginal <- array(NA, dim = c(npoints, 2, 4))
    dimnames(marginal) <- list(NULL, c("Spp1", "Spp2"),
                               c("mean", "sd", "lower", "upper"))
    marginal[,,1] <- apply(marginalMC, 2:3, mean)
    marginal[,,2] <- apply(marginalMC, 2:3, sd)
    marginal[,,3] <- apply(marginalMC, 2:3, quantile, probs = 0.025)
    marginal[,,4] <- apply(marginalMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute conditional occupancy for all species
    #'  --------------------------------------------- 
    #'  Occupancy probability for each species, conditional on the presence or
    #'  absence of the other species
    #'  Get MCMC chains
    conditionalMC <- array(NA, dim = c(ndraws, npoints, 4))
    dimnames(conditionalMC)[[3]] <- c("Spp1.alone", "Spp1.given.Spp2",
                                      "Spp2.alone", "Spp2.given.Spp1")
    
    conditionalMC[,,"Spp1.alone"] <- lspMC[,,"Spp1"]/apply(lspMC[,,c("U", "Spp1")], 1:2, sum)
    conditionalMC[,,"Spp1.given.Spp2"] <- lspMC[,,"Both"]/marginalMC[,,"Spp2"]
    conditionalMC[,,"Spp2.alone"] <- lspMC[,,"Spp2"]/apply(lspMC[,,c("U", "Spp2")], 1:2, sum)
    conditionalMC[,,"Spp2.given.Spp1"] <- lspMC[,,"Both"]/marginalMC[,,"Spp1"]
    
    #'  Fill in summary values
    conditional <- array(NA, dim = c(npoints, 4, 4))
    dimnames(conditional) <- list(NULL, c("Spp1.alone", "Spp1.given.Spp2",
                                          "Spp2.alone", "Spp2.given.Spp1"),
                                  c("mean", "sd", "lower", "upper"))
    conditional[,,1] <- apply(conditionalMC, 2:3, mean)
    conditional[,,2] <- apply(conditionalMC, 2:3, sd)
    conditional[,,3] <- apply(conditionalMC, 2:3, quantile, probs = 0.025)
    conditional[,,4] <- apply(conditionalMC, 2:3, quantile, probs = 0.975)
    
    #'  List predicted probabilities and new covariate data to plot
    predicted_probabilities <- list(marginal, conditional, cov.pred.orig, cov.pred)
    
    return(predicted_probabilities)
  }
  #####  Wolf-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(setup_preygroup_rx)_px(setup_effort).RData") 
  wolf.coy.smdeer.pred1 <- predict_occupancy(mod = wolf.coy.preygroup, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nsmall_deer,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0), cov_index = 5)
  wolf.coy.bgdeer.pred1 <- predict_occupancy(mod = wolf.coy.preygroup, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nbig_deer,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0), cov_index = 6)
  wolf.coy.smdeer.pred0 <- predict_occupancy(mod = wolf.coy.preygroup, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$Nsmall_deer,
                                           psi_cov = c(1, 0, 0, 0, 0, 0, 0), cov_index = 5)
  wolf.coy.bgdeer.pred0 <- predict_occupancy(mod = wolf.coy.preygroup, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nbig_deer,
                                             psi_cov = c(1, 0, 0, 0, 0, 0, 0), cov_index = 6)
  
  #####  Wolf-Lion predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(setup_anthro_rx)_px(setup_effort).RData")
  wolf.lion.elev.pred1 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$Elev,
                                           psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), cov_index = 3)
  wolf.lion.distrd.pred1 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$logNearestRd,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), cov_index = 6)
  wolf.lion.cattle.pred1 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nlivestock,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), cov_index = 8)
  wolf.lion.elev.pred0 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$Elev,
                                           psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), cov_index = 3)
  wolf.lion.distrd.pred0 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$logNearestRd,
                                             psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), cov_index = 6)
  wolf.lion.cattle.pred0 <- predict_occupancy(mod = wolf.lion.anthro, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nlivestock,
                                             psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), cov_index = 8)
  
  #####  Coyote-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort).RData")
  coy.bob.preyabundx$summary[1:30,]
  coy.bob.for.pred0 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$PercForest,
                                         psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 4,
                                         psi_inxs_cov = c(1, 0, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.wtd.pred0 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Nwtd,
                                         psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 5,
                                         psi_inxs_cov = c(1, 0, 0, 0), psi_inxs_cov_index = 3)
  coy.bob.lago.pred0 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Nlagomorph,
                                          psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 6,
                                          psi_inxs_cov = c(1, 0, 0, 0), psi_inxs_cov_index = 4)
  
  coy.bob.for.pred1 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$PercForest,
                                         psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 4,
                                         psi_inxs_cov = c(1, 1, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.wtd.pred1 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Nwtd,
                                         psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 5,
                                         psi_inxs_cov = c(1, 1, 0, 0), psi_inxs_cov_index = 3)
  coy.bob.lago.pred1 <- predict_occupancy(mod = coy.bob.preyabundx, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Nlagomorph,
                                          psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 6,
                                          psi_inxs_cov = c(1, 1, 0, 0), psi_inxs_cov_index = 4)
  
  
  #'  --------------------------------
  ####  Plot marginal Pr(occupancy)  ####
  #'  --------------------------------
  #'  Function to plot marginal probabilities
  marginal_prob_plots <- function(focal_cov, marg, covariate_name, spp1, spp2, plotit = T) { 
    #'  Plot marginal probability for each species
    matplot(focal_cov, marg[,,"mean"],
            type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, xlab = covariate_name, #col = 1:2
            ylab = "Marginal occupancy probability", ylim = c(0, 1), #ylim = c(0, max(marg[,,"upper"]))
            main = paste("Marginal Occupancy for", spp1, "and", spp2))
    #'  Add CRIs
    polygon(c(focal_cov, rev(focal_cov)), c(marg[,1,"lower"], rev(marg[,1,"upper"])),
            border = NA, col = adjustcolor("black", alpha = 0.2))
    polygon(c(focal_cov, rev(focal_cov)), c(marg[,2,"lower"], rev(marg[,2,"upper"])),
            border = NA, col = adjustcolor("darkblue", alpha = 0.2))
    # for(i in 1:2) {
    #   polygon(c(focal_cov, rev(focal_cov)), c(marg[,i,"lower"], rev(marg[,i,"upper"])),
    #           border = NA, col = adjustcolor(i, alpha = 0.2))
    # }
    legend("topright", c(spp1, spp2), lwd = 3, lty = 1:2, col = c("black", "darkblue"), horiz = TRUE, bty = "n") #col = 1:2
  }
  #####  Wolf-Lion marginal occupancy  ####
  marginal_prob_plots(focal_cov = wolf.lion.elev.pred1[[3]], 
                      marg = wolf.lion.elev.pred1[[1]], covariate_name = "Elevation (m)", 
                      spp1 = "Wolf", spp2 = "Lion")
  marginal_prob_plots(focal_cov = wolf.lion.distrd.pred1[[3]], 
                      marg = wolf.lion.distrd.pred1[[1]], covariate_name = "log distance to neareast road", 
                      spp1 = "Wolf", spp2 = "Lion")
  marginal_prob_plots(focal_cov = wolf.lion.cattle.pred1[[3]], 
                      marg = wolf.lion.cattle.pred1[[1]], covariate_name = "Cattle activity", 
                      spp1 = "Wolf", spp2 = "Lion")
  marginal_prob_plots(focal_cov = wolf.lion.elev.pred0[[3]], 
                      marg = wolf.lion.elev.pred0[[1]], covariate_name = "Elevation (m)", 
                      spp1 = "Wolf", spp2 = "Lion")
  marginal_prob_plots(focal_cov = wolf.lion.distrd.pred0[[3]], 
                      marg = wolf.lion.distrd.pred0[[1]], covariate_name = "log distance to neareast road", 
                      spp1 = "Wolf", spp2 = "Lion")
  marginal_prob_plots(focal_cov = wolf.lion.cattle.pred0[[3]], 
                      marg = wolf.lion.cattle.pred0[[1]], covariate_name = "Cattle activity", 
                      spp1 = "Wolf", spp2 = "Lion")
  
  #####  Wolf-Coyote marginal occupancy  ####
  marginal_prob_plots(focal_cov = wolf.coy.smdeer.pred1[[3]], 
                      marg = wolf.coy.smdeer.pred1[[1]], covariate_name = "Deer activity", 
                      spp1 = "Wolf", spp2 = "Coyote")
  marginal_prob_plots(focal_cov = wolf.coy.bgdeer.pred1[[3]], 
                      marg = wolf.coy.bgdeer.pred1[[1]], covariate_name = "Elk/moose activity", 
                      spp1 = "Wolf", spp2 = "Coyote")
  marginal_prob_plots(focal_cov = wolf.coy.smdeer.pred0[[3]], 
                      marg = wolf.coy.smdeer.pred0[[1]], covariate_name = "Deer activity", 
                      spp1 = "Wolf", spp2 = "Coyote")
  marginal_prob_plots(focal_cov = wolf.coy.bgdeer.pred0[[3]], 
                      marg = wolf.coy.bgdeer.pred0[[1]], covariate_name = "Elk/moose activity", 
                      spp1 = "Wolf", spp2 = "Coyote")
  
  
  #####  Coyote-Bobcat marginal occupancy  ####
  coy.bob.for0_margocc <- marginal_prob_plots(focal_cov = coy.bob.for.pred0[[3]], 
                                              marg = coy.bob.for.pred0[[1]], 
                                              covariate_name = "Percent forest cover", 
                                              spp1 = "Coyote", spp2 = "Bobcat")
  coy.bob.wtd0_margocc <- marginal_prob_plots(focal_cov = coy.bob.wtd.pred0[[3]], 
                                               marg = coy.bob.wtd.pred0[[1]], 
                                               covariate_name = "White-tailed deer RAI", 
                                               spp1 = "Coyote", spp2 = "Bobcat")
  coy.bob.lago0_margocc <- marginal_prob_plots(focal_cov = coy.bob.lago.pred0[[3]],
                                              marg = coy.bob.lago.pred0[[1]], 
                                              covariate_name = "Lagomorph RAI", 
                                              spp1 = "Coyote", spp2 = "Bobcat")
  coy.bob.for1_margocc <- marginal_prob_plots(focal_cov = coy.bob.for.pred1[[3]],
                                              marg = coy.bob.for.pred1[[1]],
                                              covariate_name = "Percent forest cover",
                                              spp1 = "Coyote", spp2 = "Bobcat")
  # coy.bob.wtd1_margocc <- marginal_prob_plots(focal_cov = coy.bob.wtd.pred1[[3]], 
  #                                              marg = coy.bob.wtd.pred1[[1]], 
  #                                              covariate_name = "White-tailed deer RAI", 
  #                                              spp1 = "Coyote", spp2 = "Bobcat")
  # coy.bob.lago1_margocc <- marginal_prob_plots(focal_cov = coy.bob.lago.pred1[[3]], 
  #                                              marg = coy.bob.lago.pred1[[1]], 
  #                                              covariate_name = "Lagomorph RAI", 
  #                                              spp1 = "Coyote", spp2 = "Bobcat")
  
  
  #'  -----------------------------------
  ####  Plot conditional Pr(occupancy)  ####
  #'  -----------------------------------
  #'  Function to plot conditional probabilities
  conditional_prob_plot <- function(focal_cov, cond, covariate_name, spp1, spp2, plotit = T) { 
    #'  Plot conditional probabilities of occupancy
    op <- par(mfrow = c(1,2))
    matplot(focal_cov, cond[,1:2,"mean"],
            type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, #col = 1:2
            xlab = covariate_name, ylab = "Conditional occupancy probability",
            ylim = c(0, 1), main = spp1)
    polygon(c(focal_cov, rev(focal_cov)), c(cond[,1,"lower"], rev(cond[,1,"upper"])),
            border = NA, col = adjustcolor("black", alpha = 0.2))
    polygon(c(focal_cov, rev(focal_cov)), c(cond[,2,"lower"], rev(cond[,2,"upper"])),
            border = NA, col = adjustcolor("darkblue", alpha = 0.2))
    # for(i in 1:2){
    #   polygon(c(focal_cov, rev(focal_cov)), c(cond[,i,"lower"], rev(cond[,i,"upper"])),
    #           border = NA, col = adjustcolor(i, alpha = 0.2))
    # }
    legend("topright", c("alone", paste("given", spp2)), lwd = 3, lty = 1:3, col = c("black", "darkblue"), bty = 'n', cex = 1.5) #col = 1:2
    
    matplot(focal_cov, cond[,3:4,"mean"],
            type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, #col = 1:2
            xlab = covariate_name, ylab = "Conditional occupancy probability",
            ylim = c(0, 1), main = spp2)
    polygon(c(focal_cov, rev(focal_cov)), c(cond[,3,"lower"], rev(cond[,3,"upper"])),
            border = NA, col = adjustcolor("black", alpha = 0.2))
    polygon(c(focal_cov, rev(focal_cov)), c(cond[,4,"lower"], rev(cond[,4,"upper"])),
            border = NA, col = adjustcolor("darkblue", alpha = 0.2))
    # for(i in 1:2){
    #   polygon(c(focal_cov, rev(focal_cov)), c(cond[,i+2,"lower"], rev(cond[,i+2,"upper"])),
    #           border = NA, col = adjustcolor(i, alpha = 0.2))
    # }
    legend("topright", c("alone", paste("given", spp1)), lwd = 3, lty = 1:3, col = c("black", "darkblue"), bty = 'n', cex = 1.5) #col = 1:2
    par(op)
  }
  #####  Wolf-Lion conditional occupancy  ####
  #'  Anthropogenic model
  conditional_prob_plot(focal_cov = wolf.lion.elev.pred1[[3]], 
                        covariate_name = "Elevation", cond = wolf.lion.pred1[[2]],
                        spp1 = "wolf", spp2 = "lion")
  conditional_prob_plot(focal_cov = wolf.lion.distrd.pred1[[3]], 
                        covariate_name = "log distance to nearest road", cond = wolf.lion.distrd.pred1[[2]],
                        spp1 = "wolf", spp2 = "lion")
  conditional_prob_plot(focal_cov = wolf.lion.cattle.pred1[[3]], 
                        covariate_name = "Cattle activity", cond = wolf.lion.cattle.pred1[[2]],
                        spp1 = "wolf", spp2 = "lion")
  conditional_prob_plot(focal_cov = wolf.lion.elev.pred0[[3]], 
                        covariate_name = "Elevation", cond = wolf.lion.elev.pred0[[2]],
                        spp1 = "wolf", spp2 = "lion")
  conditional_prob_plot(focal_cov = wolf.lion.distrd.pred0[[3]], 
                        covariate_name = "log distance to nearest road", cond = wolf.lion.distrd.pred0[[2]],
                        spp1 = "wolf", spp2 = "lion")
  conditional_prob_plot(focal_cov = wolf.lion.cattle.pred0[[3]], 
                        covariate_name = "Cattle activity", cond = wolf.lion.cattle.pred0[[2]],
                        spp1 = "wolf", spp2 = "lion")
  
  
  #####  Wolf-Coyote conditional occupancy  ####
  #'  Prey group model
  conditional_prob_plot(focal_cov = wolf.coy.smdeer.pred1[[3]], 
                        covariate_name = "Deer activity", cond = wolf.coy.smdeer.pred1[[2]],
                        spp1 = "wolf", spp2 = "coyote")
  conditional_prob_plot(focal_cov = wolf.coy.bgdeer.pred1[[3]], 
                        covariate_name = "Elk/moose activity", cond = wolf.coy.bgdeer.pred1[[2]],
                        spp1 = "wolf", spp2 = "coyote")
  conditional_prob_plot(focal_cov = wolf.coy.smdeer.pred0[[3]], 
                        covariate_name = "Deer activity", cond = wolf.coy.smdeer.pred0[[2]],
                        spp1 = "wolf", spp2 = "coyote")
  conditional_prob_plot(focal_cov = wolf.coy.bgdeer.pred0[[3]], 
                        covariate_name = "Elk/moose activity", cond = wolf.coy.bgdeer.pred0[[2]],
                        spp1 = "wolf", spp2 = "coyote")
  
  
  #####  Coyote-Bobcat conditional occupancy  ####
  coy.bob.wtd0_coocc <- conditional_prob_plot(focal_cov = coy.bob.wtd.pred0[[3]],
                        covariate_name = "White-tailed deer RAI", cond = coy.bob.wtd.pred0[[2]],
                        spp1 = "coyote", spp2 = "bobcat")
  coy.bob.lago0_coocc <- conditional_prob_plot(focal_cov = coy.bob.lago.pred0[[3]],
                        covariate_name = "Lagomorph RAI", cond = coy.bob.lago.pred0[[2]],
                        spp1 = "coyote", spp2 = "bobcat")
  # coy.bob.wtd1_coocc <- conditional_prob_plot(focal_cov = coy.bob.wtd.pred1[[3]],
  #                       covariate_name = "White-tailed deer RAI", cond = coy.bob.wtd.pred1[[2]],
  #                       spp1 = "coyote", spp2 = "bobcat")
  # coy.bob.lago1_coocc <- conditional_prob_plot(focal_cov = coy.bob.lago.pred1[[3]], 
  #                       covariate_name = "Lagomorph RAI", cond = coy.bob.lago.pred1[[2]],
  #                       spp1 = "coyote", spp2 = "bobcat")
  
  
  
  
  
  
  
  