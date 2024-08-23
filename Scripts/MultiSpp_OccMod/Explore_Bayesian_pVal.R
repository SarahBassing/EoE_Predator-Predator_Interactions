  #'  Review GoF outputs ---- do these models really suck that bad???
  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_altGoF_2024-08-21.RData") 
  mod <- bear.coy.habx
  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_altGoF_2024-08-21.RData")
  mod <- wolf.coy.hab
  
  # load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_GoF_2024-08-07.RData")
  # mod <- wolf.bear.null
  # 
  # load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_GoF_2024-08-07.RData")
  # mod <- wolf.bear.hab
  # 
  # load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_GoF_2024-08-08.RData")
  # mod <- wolf.lion.null
  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_altGoF_2024-08-23.RData") # need to run without z.sim & fixed FT
  mod <- coy.bob.habx
  
  mod$summary
  
  #'  Observed and simulated data
  y.hat <- mod$sims.list$y.hat; y.hat[15000,1,,]; hist(y.hat[1,1:1138,,])
  y.hat.index <- mod$sims.list$y.hat.index; y.hat.index[15000,1,,]; hist(y.hat[1,1:1138,,])
  y.hat.maxindex <- mod$sims.list$y.hat.maxindex; y.hat.maxindex[100,1:100,]; hist(y.hat.maxindex[1,1:1138,])
  chi2.obs <- mod$sims.list$chi2.obs; chi2.obs[1:10]; hist(chi2.obs); summary(chi2.obs)
  chi2.sim <- mod$sims.list$chi2.sim; chi2.sim[1:10]; hist(chi2.sim); summary(chi2.sim)
  (mean(chi2.sim > chi2.obs))
  
  pl <- range(c(chi2.obs, chi2.sim))
  plot(chi2.obs, chi2.sim, xlab = "Chi2 observed data", ylab = "Chi2 expected data",
       main = "Chi2 discrepency ratio", xlim = pl, ylim = pl, frame.plot = FALSE)
  abline(0,1, lwd = 2) #1:1 ratio would be perfect fit
  text(quantile(pl, 0.1), max(pl), paste('Bpv = ', round(mean(chi2.sim > chi2.obs), 2)), cex = 1.5)
  
  #' y <- mod$sims.list$y; y[1, 1:10, ]; hist(y[1,,])
  #' y.sim <- mod$sims.list$y.sim; y.sim[1, 1:10, ]; hist(y.sim[1,,])
  #' 
  #' #'  Observed data broken down by species
  #' y_A <- mod$sims.list$y_A; y_A[1, 1:10, ]; hist(y_A[1,,])
  #' y_B <- mod$sims.list$y_B; y_B[1, 1:10, ]; hist(y_B[1,,])
  #' 
  #' #'  Simulated data broken down by species
  #' yrep_A <- mod$sims.list$yrep_A; yrep_A[1, 1:10, ]; hist(yrep_A[1,,])
  #' yrep_B <- mod$sims.list$yrep_B; yrep_B[1, 1:10, ]; hist(yrep_B[1,,])
  #' 
  #' #'  Observed detection frequencies per species (sum number of dets per site and iteration)
  #' detfreq_A <- mod$sims.list$detfreq_A; detfreq_A[1:10, 1:10]; hist(detfreq_A)
  #' detfreq_B <- mod$sims.list$detfreq_B; detfreq_B[1:10, 1:10]; hist(detfreq_B)
  #' 
  #' #'  Simulated detection frequencies per species based on model
  #' detfreqrep_A <- mod$sims.list$detfreqrep_A; detfreqrep_A[1:10, 1:10]; hist(detfreqrep_A)
  #' detfreqrep_B <- mod$sims.list$detfreqrep_B; detfreqrep_B[1:10, 1:10]; hist(detfreqrep_B)
  #' 
  #' #'  Expected detection matrix (z[i]*rdm[i,j]) per species
  #' tmp_A <- mod$sims.list$tmp_A; tmp_A[1, 1:10, ]; hist(tmp_A[1:1000,,])   # how are there values >1 is rdm should be (0-1)???
  #' tmp_B <- mod$sims.list$tmp_B; tmp_B[1, 1:10, ]; hist(tmp_B[1:1000,,])
  #' 
  #' #'  Expected detection frequency per species (sum expected number of dets per site and iteration)
  #' E_A <- mod$sims.list$E_A; E_A[1:10, 1:10]; hist(E_A)
  #' E_B <- mod$sims.list$E_B; E_B[1:10, 1:10]; hist(E_B)
  
  #'  Chi2 discrepancy measure between observed and expected 
  x2_A <- mod$sims.list$x2_A; x2_A[1:10, 1:10]; hist(x2_A); max(x2_A)
  x2_B <- mod$sims.list$x2_B; x2_B[1:10, 1:10]; hist(x2_B); max(x2_B)
  
  #'  Chi2 discrepency measure between simulated and expected
  x2rep_A <- mod$sims.list$x2rep_A; x2rep_A[1:10, 1:10]; hist(x2rep_A); max(x2rep_A)
  x2rep_B <- mod$sims.list$x2rep_B; x2rep_B[1:10, 1:10]; hist(x2rep_B); max(x2rep_B)
  
  #'  Overall Chi2 test statistic for observed data (sum x2 measures across sites per interaction)
  chi2.obs_A <- mod$sims.list$chi2.obs_A; chi2.obs_A[1:10]; hist(chi2.obs_A); summary(chi2.obs_A)
  chi2.obs_B <- mod$sims.list$chi2.obs_B; chi2.obs_B[1:10]; hist(chi2.obs_B); summary(chi2.obs_B)
  
  #'  Overall Chi2 test statistic for simulated data
  chi2.sim_A <- mod$sims.list$chi2.sim_A; chi2.sim_A[1:10]; hist(chi2.sim_A); summary(chi2.sim_A)
  chi2.sim_B <- mod$sims.list$chi2.sim_B; chi2.sim_B[1:10]; hist(chi2.sim_B); summary(chi2.sim_B)
  
  #'  Chi2 Bayesian p-value for each species
  (mean(chi2.sim_A > chi2.obs_A))
  (mean(chi2.sim_B > chi2.obs_B))
  
  #'  Chi2 fit statistic ratio
  (mean(chi2.obs_A/chi2.sim_A))
  (mean(chi2.obs_B/chi2.sim_B))
  
  #'  Overall Freeman-Tukey test statistic for observed data
  ft.obs_A <- mod$sims.list$ft.obs_A; ft.obs_A[1:10]; hist(ft.obs_A); summary(ft.obs_A)
  ft.obs_B <- mod$sims.list$ft.obs_B; ft.obs_B[1:10]; hist(ft.obs_B); summary(ft.obs_B)
  
  #'  Overall Freeman-Tukey test statistic for simulated data
  ft.sim_A <- mod$sims.list$ft.sim_A; ft.sim_A[1:10]; hist(ft.sim_A); summary(ft.sim_A)
  ft.sim_B <- mod$sims.list$ft.sim_B; ft.sim_B[1:10]; hist(ft.sim_B); summary(ft.sim_B)
  
  #'  FT Bayesian p-value
  (mean(ft.sim_A > ft.obs_A))
  (mean(ft.sim_B > ft.obs_B))
  
  #'  FT fit statistic ratio
  (mean(ft.obs_A/ft.sim_A))
  (mean(ft.obs_B/ft.sim_B))
  
  
  #'  Plot Chi2 discrepancy measures for observed and simulated data for SppA
  pl <- range(c(chi2.obs_A, chi2.sim_A))
  plot(chi2.obs_A, chi2.sim_A, xlab = "Chi2 observed data", ylab = "Chi2 expected data",
       main = "Chi2 discrepency ratio for Species A", xlim = pl, ylim = pl, frame.plot = FALSE)
  abline(0,1, lwd = 2) #1:1 ratio would be perfect fit
  text(quantile(pl, 0.1), max(pl), paste('Bpv = ', round(mean(chi2.sim_A > chi2.obs_A), 2)), cex = 1.5)
  
  #'  And SppB
  pl <- range(c(chi2.obs_B, chi2.sim_B))
  plot(chi2.obs_B, chi2.sim_B, xlab = "Chi2 observed data", ylab = "Chi2 expected data",
       main = "Chi2 discrepancy ratio for Species B", xlim = pl, ylim = pl, frame.plot = FALSE)
  abline(0,1, lwd = 2)
  text(quantile(pl, 0.1), max(pl), paste('Bpv = ', round(mean(chi2.sim_B > chi2.obs_B), 2)), cex = 1.5)
  
  #'  Plot Freeman-Tukey discrepancy measures for observed and simulated data for SppA
  pl <- range(c(ft.obs_A, ft.sim_A))
  plot(ft.obs_A, ft.sim_A, xlab = "FT observed data", ylab = "FT expected data",
       main = "FT discrepency ratio for Species A", xlim = pl, ylim = pl, frame.plot = FALSE)
  abline(0,1, lwd = 2)
  text(quantile(pl, 0.1), max(pl), paste('Bpv = ', round(mean(ft.sim_A > ft.obs_A), 2)), cex = 1.5)

  #'  and SppB
  pl <- range(c(ft.obs_B, ft.sim_B))
  plot(ft.obs_B, ft.sim_B, xlab = "FT observed data", ylab = "FT expected data",
       main = "FT discrepancy ratio for Species B", xlim = pl, ylim = pl, frame.plot = FALSE)
  abline(0,1, lwd = 2)
  text(quantile(pl, 0.1), max(pl), paste('Bpv = ', round(mean(ft.sim_B > ft.obs_B), 2)), cex = 1.5)

  
  