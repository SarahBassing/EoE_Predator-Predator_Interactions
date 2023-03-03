  #'  ------------------------------
  #'  Plot co-occurrence results
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------
  #'  Script to create figures of predicted co-occurrence and covariate effects
  #'  from multispecies occupancy models run in JAGS.
  #'  -------------------------------
  
  #'  Load libraries
  library(ggplot2)
  
  #'  Load saved JAGS output and covariate data
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preydiversity)_p(setup_effort).RData") #coybob_psi(habitat)_p(effort)
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s.RData")
  
  out1 <- coy.bob.preydiversity
  
  #'  Define number of states
  ncat <- 4
  #'  Number of draws in the MCMC chain
  ndraws <- out1$mcmc.info$n.samples
  #'  Number of expected observations to predict
  npoints <- 1000
  
  #'  Grab posterior means of parameters
  # (tmp <- out1$mean[1:5])
  print(out1$summary[1:44,])
  
  #'  Create prediction covariates for occupancy as a function of focal covariate
  (r <- range(stations_skinny_eoe20s21s$Nwtd))
  (mean.cov <- mean(stations_skinny_eoe20s21s$Nwtd))
  (sd.cov<- sd(stations_skinny_eoe20s21s$Nwtd))
  #'  Create new data based on range of sampled covariate data
  cov.pred.orig <- seq(r[1], r[2], length.out = npoints) 
  #'  Scale new data to be consistent with model input
  cov.pred <- (cov.pred.orig - mean.cov) / sd.cov
  
  #'  Create model matrices for the prediction covariates in order of model coefficients
  #'  Vary one covariate while holding others constant
  head(psi_cov <- cbind(1, 1, cov.pred, 0, 0, 0, 0, 0, 0)) #' intercept, setup, elevation, forest, elk, moose, md, wtd, lagomorph
  
  #'  Assemble linear predictors across all iterations of the MCMC chains
  psiSpp1 <- out1$sims.list$betaSpp1 %*% t(psi_cov)
  psiSpp2 <- out1$sims.list$betaSpp2 %*% t(psi_cov)
  psiSpp12 <- psiSpp1 + psiSpp2 + out1$sims.list$betaSpp12[,1:9] %*% t(psi_cov) #drop the indixing once I rerun without the unused priors
  
  ####  Assemble the state probability array  ####
  #'  Compute latent state vector (lsv)
  lsv <- array(NA, dim = c(ndraws, npoints, ncat))
  lsv[,,1] <- 1
  lsv[,,2] <- exp(psiSpp1)
  lsv[,,3] <- exp(psiSpp2)
  lsv[,,4] <- exp(psiSpp12)

  #'  Probability of each state as MCMC chains (lsp)
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
  
  #'  Plot the latent state probabilities as a function of covariate
  matplot(cov.pred.orig, lsp[,,"mean"], type = "l", lty = 1, lwd = 5, col = 1:4,
          frame = FALSE, xlab = "Covariate", ylab = "Probability of State",
          ylim = c(0, max(lsp[,,-2])), las = 1)
  #'  Add CRIs
  for(i in 1:4){
    polygon(c(cov.pred.orig, rev(cov.pred.orig)), c(lsp[,i,"lower"], rev(lsp[,i,"upper"])),
            border = NA, col = adjustcolor(i, alpha = 0.2))
  }
  legend("topright", lwd = 3, lty = 1, col = 1:4, colnames(lsp), bty = "n", cex = 1.2, ncol = 2)
  
  #'  Or plot them separately
  op <- par(mfrow=c(2,2))
  for(i in 1:4) {
    plot(cov.pred.orig, lsp[,i,'mean'], type = 'n',
         frame = FALSE, xlab = 'Covariate', ylab = 'Probability of state',
         ylim = c(0, max(lsp[,,-2])), las = 1)
    polygon(c(cov.pred.orig, rev(cov.pred.orig)), c(lsp[,i,'lower'], rev(lsp[,i,'upper'])),
            border=NA, col=adjustcolor(i, alpha=0.2))
    lines(cov.pred.orig, lsp[,i,'mean'], lwd = 2, col = i)
    legend('topright', legend=colnames(lsp)[i], lwd = 2, col = i, bty = 'n')
  }
  par(op)
  
  ####  Compute and plot marginal probabilities  ####
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
  
  #'  Plot marginal occupancy for each species
  matplot(cov.pred.orig, marginal[,,"mean"],
          type = "l", lty = 1:2, lwd = 3, col = 1:2, frame = FALSE, xlab = "Covariate",
          ylab = "Marginal occupancy probability", ylim = c(0, max(marginal[,,"upper"])))
  #'  Add CRIs
  for(i in 1:2) {
    polygon(c(cov.pred.orig, rev(cov.pred.orig)), c(marginal[,i,"lower"], rev(marginal[,i,"upper"])),
            border = NA, col = adjustcolor(i, alpha = 0.2))
  }
  legend("bottom", c("Spp1", "Spp2"), lwd = 3, lty = 1:2, col = 1:2, horiz = TRUE, bty = "n")
  
  ####  Compute conditional occupancy for all species  ####
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
  
  #'  Plot conditional probabilities of occupancy
  op <- par(mfrow = c(1,2))
  matplot(cov.pred.orig, conditional[,1:2,"mean"],
          type = "l", lty = 1:2, lwd = 3, col = 1:2, frame = FALSE,
          xlab = "Covariate", ylab = "Occupancy probability",
          ylim = c(0, 1), main = "Spp1")
  for(i in 1:2){
    polygon(c(cov.pred.orig, rev(cov.pred.orig)), c(conditional[,i,"lower"], rev(conditional[,i,"upper"])),
            border = NA, col = adjustcolor(i, alpha = 0.2))
  }
  legend("bottomright", c("alone", "given Spp2"), lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
  
  matplot(cov.pred.orig, conditional[,3:4,"mean"],
          type = "l", lty = 1:2, lwd = 3, col = 1:2, frame = FALSE,
          xlab = "Covariate", ylab = "Occupancy probability",
          ylim = c(0, 1), main = "Spp2")
  for(i in 1:2){
    polygon(c(cov.pred.orig, rev(cov.pred.orig)), c(conditional[,i+2,"lower"], rev(conditional[,i+2,"upper"])),
            border = NA, col = adjustcolor(i, alpha = 0.2))
  }
  legend("topright", c("alone", "given Spp1"), lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
  par(op)
  
  
  
  
  
  