  #'  -----------------------------------
  #'  D-separation tests with ROPE method
  #'  Sarah Bassing & Matt Falcy
  #'  July 2026
  #'  -----------------------------------
  #'  Using Region of Practical Equivalence (ROPE) to test d-sep claims. This is 
  #'  a bit of a hack to generate pseudo p-values for Fisher's C. Method adapted
  #'  from Kruschke and Liddell (2018), who used "0.1" of 1 SD to define a ROPE:
  #'  https://link.springer.com/article/10.3758/s13423-016-1221-4
  #'  Make sure data have been properly formatted using Bayesian_SEMs_relative_density_index_1yLag.R
  #'  -----------------------------------
  
  #'  Load libries
  install.packages("BiocManager") # Needed for 'graph' which is needed by ggm
  BiocManager::install("graph")
  library(ggm) # for DAG() and basiSet(), which retrieves the basis set
  library(jagsUI)
  #library(R2jags)
  
  #'  Define function to establish a region of practical equivalence (ROPE) around
  #'  the null value. This expresses a range of parameter values considered equivalent
  #'  to the null value. This range is -pct*sd(y) to pct*sd(y)
  #'  If 95% CRI falls entirely within this range, 95% of the posterior is pratically
  #'  equivalent to the null value. Taking the mean here calculates the average 
  #'  portion of the posterior that falls outside the ROPE. Larger p.rope values
  #'  (essentially p-values) indicate.....
  p.rope <- function(pct = 0.1, y = NULL, post = NULL){
    mean(-pct*sd(y) < post & pct*sd(y) > post)
  }
  
  #'  Generate DAG for topdown model (exclude time t intercepts b/c these are 
  #'  statistically necessary but not causal)
  dag_topdown <- DAG(mu.wolf.t ~ beta.wolf[1] + beta.harvest[1], #beta.int[2] + 
                     # mu.wolf.tmin1 ~ beta.int.tmin1[2],
                     mu.lion.t ~ beta.lion[1] + beta.harvest[2], #beta.int[1] + 
                     # mu.lion.tmin1 ~ beta.int.tmin1[1],
                     mu.bear.t ~ beta.bear[1] * bear.tmin1[i] + beta.harvest[3], #beta.int[3] + 
                     # mu.bear.tmin1 ~ beta.int.tmin1[3], 
                     mu.coy.t ~ beta.coy[1], #beta.int[4] + 
                     # mu.coy.tmin1 ~ beta.int.tmin1[4], 
                     mu.elk.t ~ beta.elk[1] + beta.wolf[2] + beta.lion[2] + beta.bear[2] + beta.harvest[4], #beta.int[5] + 
                     # mu.elk.tmin1 ~ beta.int.tmin1[5],
                     mu.moose.t ~ beta.moose[1] + beta.wolf[3], #beta.int[6] + 
                     # mu.moose.tmin1 ~ beta.int.tmin1[6],
                     mu.wtd.t ~ beta.wtd[1] + beta.wolf[4] + beta.lion[3] + beta.bear[3] + beta.coy[2] + beta.harvest[6])#, #beta.int[7] + 
                     # mu.wtd.tmin1 ~ beta.int.tmin1[7])
  
  #'  Generate the basic set for this model
  bs_topdown <- basiSet(dag_topdown); length(bs_topdown)
  #'  Return basic set with conditional independence claims involving 3+ variables
  bs_topdown_3plus <- sapply(bs_topdown, length) >= 3
  bs_topdown_skinny <- bs_topdown[bs_topdown_3plus]; length(bs_topdown_skinny)
  View(bs_topdown_skinny)
  #'  Only testing conditional independence claims that are consistent with a priori hypotheses or biologically plausible
  
  #'  Call JAGS and fit model (MUST RUN LINES 1 - 210 IN Bayesian_SEMs_relative_density_index_1yLag.R)
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_dsep.R")
  SEM_topdown_dsep <- jags(data_JAGS_bundle_top, inits = initsList_topdown, params, 
                           "./Outputs/SEM/JAGS_out/d_Sep/JAGS_SEM_topdown_dsep.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)

  #'  bs_topdown_skinny[[7]]: beta.coy[2] and mu.moose.t are independent given beta.wolf[3] and beta.moose[1]
  #'  Apply ROPE method to standardized original data and estimated parameter posterior of interest
  bs_topdown_skinny[[7]]
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = SEM_topdown_dsep$sims.list$beta.coy)  # 0.1275667
  
  #'  bs_topdown_skinny[[8]]: beta.coy[2] and mu.elk.t are independent given beta.harvest[4], beta.bear[2], beta.lion[2], beta.wolf[2], and beta.elk[1]
  bs_topdown_skinny[[8]]
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = SEM_topdown_dsep$sims.list$beta.coy)  # 0.06733333
  
  
  
  
  #bs_topdown: "mu.wtd.tmin1" and "beta.harvest[6]" are independent given "beta.int.tmin1[7]".  Let's test:
  
  #'  Call JAGS and fit model (MUST HAVE ALL JAGS PARAMETERS RUN IN Bayesian_SEMs_relative_density_index_1yLag.R)
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_dsep.R")
  SEM_topdown_dsep <- jags(data_JAGS_bundle_top, inits = initsList_topdown, params, 
                           "./Outputs/SEM/JAGS_out/d_Sep/JAGS_SEM_topdown_dsep.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  #'  Apply ROPE method to standardized original data and estimated parameter posterior of interest
  p.T_bs_topdown <- p.rope(y = data_JAGS_bundle_top$wtd.tmin1_hat, post = SEM_topdown_dsep$sims.list$beta.harvest)
  print(p.T_bs_topdown)
  
  
  
  