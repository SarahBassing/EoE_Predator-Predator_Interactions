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
  # install.packages("BiocManager") # Needed for 'graph' which is needed by ggm
  # BiocManager::install("graph")
  library(ggm) # for DAG() and basiSet(), which retrieves the basis set
  library(jagsUI)
  #library(R2jags)
  
  #'  Define function to establish a region of practical equivalence (ROPE) around
  #'  the null value. This expresses a range of parameter values considered equivalent
  #'  to the null value (0). This range is -pct*sd(y) to pct*sd(y)
  #'  If 95% CRI falls entirely within this range, 95% of the posterior is practically
  #'  equivalent to the null value. Taking the mean calculates the proportion
  #'  of the posterior that falls inside the ROPE. In other words, what is the 
  #'  probability that the effect is small enough to be ignored? Larger p.rope 
  #'  values indicate a large percentage of the posterior is "practically" 0; 
  #'  these variables are conditionally independent. Small p.rope values (i.e., < 0.05) 
  #'  indicate a small percentage of the posterior falls inside the ROPE and thus  
  #'  it is not practically equivalent to 0 (i.e., it is "practically" significant).
  #'  These variables are NOT conditionally independent and this relationship is 
  #'  missing in the model. p.rope can be treated like a p-value in this context.
  p.rope <- function(pct = 0.1, y = NULL, post = NULL){
    #'  Average number of posterior draws that fall above lower ROPE value AND below upper ROPE value
    mean(-pct*sd(y) < post & pct*sd(y) > post)
  }
  
  
  #'  Generate DAG for top-down model 
  dag_topdown <- DAG(wolf.t ~ wolf.tmin1 + wolfHarv.tmin1,  
                     lion.t ~ lion.tmin1 + lionHarv.tmin1,  
                     bear.t ~ bear.tmin1 + bearHarv.tmin1,  
                     coy.t ~ coy.tmin1,  
                     elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1 + bear.tmin1 + elkHarv.tmin1,  
                     moose.t ~ moose.tmin1 + wolf.tmin1,  
                     wtd.t ~ wtd.tmin1 + wolf.tmin1 + lion.tmin1 + bear.tmin1 + coy.tmin1 + deerHarv.tmin1) 
  
  #'  Generate the basic set for this model
  bs_topdown <- basiSet(dag_topdown); length(bs_topdown)
  #'  Return basic set with conditional independence claims involving 3+ variables
  bs_topdown_3plus <- sapply(bs_topdown, length) >= 3
  bs_topdown_skinny <- bs_topdown[bs_topdown_3plus]; length(bs_topdown_skinny)
  View(bs_topdown_skinny)
  #'  Only testing conditional independence claims that are consistent with a priori hypotheses or biologically plausible
  
  #'  NOTE: MUST RUN LINES 1 - 210 IN Bayesian_SEMs_relative_density_index_1yLag.R BEFORE THIS
  #'  Update bundled data for JAGS
  #'  d-Sep test requires a few extra priors for some parameters that weren't needed
  #'  in the original model so rerunning code to update the bundled data
  data_JAGS_bundle_top <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 4, nlion = 3, nbear = 3, 
                                     ncoy = 2, nelk = 2, nmoose = 2, nwtd = 2, nharv = 6, nfor = 0, nwsi = 0)
  for(i in 1:num.chains){
    initsList_topdown[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 3, ncoy = 2, nelk = 2, 
                                             nmoose = 2, nwtd = 2, nharv = 6, nfor = 0, nwsi = 0)
  }
  
  #'  Function to source JAGS script and fit model 
  #'  NOTE: The JAGS script must be updated for each conditional independence claim 
  #'  below. Comment out all unneeded independence claims for each test below. This
  #'  is tedious but better than having 100s of the same script with 1 line of code difference.
  run_jags <- function() {
    source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_dsep.R")
    SEM_topdown_dsep <- jags(data_JAGS_bundle_top, inits = initsList_topdown, params, 
                             "./Outputs/SEM/JAGS_out/d_Sep/JAGS_SEM_topdown_dsep.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
    return(SEM_topdown_dsep)
  }

  #'  Review specific conditional ind. claim and update JAGS_SEM_topdown_dsep.R script
  bs_topdown_skinny[[1]] # deerHarv.tmin1 and coy.t are independent given coy.tmin1
  #'  Call run_jags function to source JAGS script updated for this conditional ind. claim and fit model
  d_Sep1 <- run_jags()
  #'  Apply ROPE method to standardized original data and posterior of interest
  #'  NOTE THE INDEXING ON POSTERIOR! Make sure to grab correct posterior for each test 
  #'  (JAGS script written flexibly to estimate posteriors for >1 parameter with 
  #'  same name that go unused in this d-Sep testing process)
  p.rope(y = data_JAGS_bundle_top$coy.t_hat, post = d_Sep1$sims.list$beta.harvest[,1])  # 0.3752
  
  bs_topdown_skinny[[2]] # deerHarv.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
  d_Sep2 <- run_jags()
  print(SEM_topdown_dsep$summary)
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep2$sims.list$beta.harvest[,2])  # 0.2181333
  
  bs_topdown_skinny[[3]] # deerHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep3 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep3$sims.list$beta.harvest[,2])  # 0.1299333
  
  bs_topdown_skinny[[4]] # deerHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep4 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep4$sims.list$beta.harvest[,1])  # 0.2194
  
  bs_topdown_skinny[[5]] # deerHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep5 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep5$sims.list$beta.harvest[,2])  # 0.2586
  
  bs_topdown_skinny[[6]] # deerHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep6 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep6$sims.list$beta.harvest[,2])  # 0.03406667
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[7]] # wtd.tmin1 and coy.t are independent given coy.tmin1
  d_Sep7 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$coy.t_hat, post = d_Sep7$sims.list$beta.wtd[,1])  # 0.1310667
  
  bs_topdown_skinny[[8]] # wtd.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
  d_Sep8 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep8$sims.list$beta.wtd[,1])  # 0.3638667
  
  bs_topdown_skinny[[9]] # wtd.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep9 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep9$sims.list$beta.wtd[,1])  # 0.3979333
  
  bs_topdown_skinny[[10]] # wtd.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep10 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep10$sims.list$beta.wtd[,1])  # 0.3182
  
  bs_topdown_skinny[[11]] # wtd.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep11 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep11$sims.list$beta.wtd[,1])  # 0.2496667
  
  bs_topdown_skinny[[12]] # wtd.tmin1 and wolf.t are independent given wolfHarv.tmin and wolf.tmin1
  d_Sep12 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep12$sims.list$beta.wtd[,1])  # 0.3222
  
  bs_topdown_skinny[[13]] # moose.tmin1 and coy.t are independent given coy.tmin1
  d_Sep13 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$coy.t_hat, post = d_Sep13$sims.list$beta.moose[,1])  # 0.1149333
  
  bs_topdown_skinny[[14]] # moose.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
  d_Sep14 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep14$sims.list$beta.moose[,1])  # 0.1387333
  
  bs_topdown_skinny[[15]] # moose.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep15 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep15$sims.list$beta.moose[,1])  # 0.3180667
  
  bs_topdown_skinny[[16]] # moose.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep16 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep16$sims.list$beta.moose[,1])  # 0.01773333
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[17]] # moose.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep17 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep17$sims.list$beta.moose[,1])  # 0.08826667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[18]] # moose.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep18 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep18$sims.list$beta.moose[,1])  # 0.1264667
  
  bs_topdown_skinny[[19]] # elkHarv.tmin1 and coy.t are independent given coy.tmin1
  d_Sep19 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$coy.t_hat, post = d_Sep19$sims.list$beta.harvest[,1])  # 0.09426667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[20]] # elkHarv.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
  d_Sep20 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep20$sims.list$beta.harvest[,2])  # 0.08666667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[21]] # elkHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep21 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep21$sims.list$beta.harvest[,2])  # 0.1417333
  
  bs_topdown_skinny[[22]] # elkHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1 
  d_Sep22 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep22$sims.list$beta.harvest[,2])  # 0.06093333
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[23]] # elkHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1 
  d_Sep23 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep23$sims.list$beta.harvest[,1])  # 0.0906
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[24]] # elkHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1 
  d_Sep24 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep24$sims.list$beta.harvest[,2])  # 0.07266667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[25]] # elk.tmin1 and coy.t are independent given coy.tmin1 
  d_Sep25 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$coy.t_hat, post = d_Sep25$sims.list$beta.elk[,1])  # 0.1714667
  
  bs_topdown_skinny[[26]] # elk.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1 
  d_Sep26 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep26$sims.list$beta.elk[,1])  # 0.2322
  
  bs_topdown_skinny[[27]] # elk.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1 
  d_Sep27 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep27$sims.list$beta.elk[,1])  # 0.361
  
  bs_topdown_skinny[[28]] # elk.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep28 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep28$sims.list$beta.elk[,1])  # 0.05626667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[29]] # elk.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep29 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep29$sims.list$beta.elk[,1])  # 0.2372667
  
  bs_topdown_skinny[[30]] # elk.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep30 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep30$sims.list$beta.elk[,1])  # 0.168
  
  bs_topdown_skinny[[31]] # coy.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
  d_Sep31 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep31$sims.list$beta.coy[,1])  # 0.2836667
  
  bs_topdown_skinny[[32]] # coy.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep32 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep32$sims.list$beta.coy[,1])  # 0.3178667
  
  bs_topdown_skinny[[33]] # coy.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep33 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep33$sims.list$beta.coy[,1])  # 0.2546
  
  bs_topdown_skinny[[34]] # coy.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep34 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep34$sims.list$beta.coy[,1])  # 0.1318667
  
  bs_topdown_skinny[[35]] # coy.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep35 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep35$sims.list$beta.coy[,1])  # 0.1266
  
  #' #'  Excluding this test b/c bearHarv.tmin1 is a purely exogenous variable 
  #' bs_topdown_skinny[[36]] # coy.tmin1 and bearHarv.tmin1 are independent given coy.tmin1
  #' d_Sep36 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$bearHarv.tmin1, post = d_Sep36$sims.list$beta.coy[,1])   
  
  #' #'  Excluding this test b/c not biologically possible (time t cannot affect time t-1)
  #' bs_topdown_skinny[[37]] # coy.t and bear.tmin1 are independent given coy.tmin1
  #' d_Sep37 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$bear.tmin1_hat, post = d_Sep37$sims.list$beta.coy[,1])
  
  bs_topdown_skinny[[38]] # coy.t and bear.t are independent given coy.tmin1, bearHarv.tmin1 and bear.tmin1
  d_Sep38 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$bear.t_hat, post = d_Sep38$sims.list$beta.coy[,2])  # 0.001933333
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  #' #'  Excluding this test b/c not biologically possible (time t cannot affect time t-1)
  #' bs_topdown_skinny[[39]] # coy.t and lionHarv.tmin1 are independent given coy.tmin1
  #' d_Sep39 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$lionHarv.tmin1, post = d_Sep39$sims.list$beta.coy[,2])
  #' 
  #' bs_topdown_skinny[[40]] # coy.t and lion.tmin1 are independent given coy.tmin1
  #' d_Sep40 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$lion.tmin1_hat, post = d_Sep40$sims.list$beta.coy[,2])  
  
  bs_topdown_skinny[[41]] # coy.t and lion.t are independent given lionHarv.tmin1, lion.tmin1 and coy.tmin1
  d_Sep41 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep41$sims.list$beta.coy[,2])  # 0
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  #' bs_topdown_skinny[[42]] # coy.t and wolfHarv.tmin1 are independent given coy.tmin1
  #' d_Sep42 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$wolfHarv.tmin1, post = d_Sep42$sims.list$beta.coy[,2])  
  #' 
  #' bs_topdown_skinny[[43]] # coy.t and wolf.tmin1 are independent given coy.tmin1
  #' d_Sep43 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$wolf.tmin1_hat, post = d_Sep43$sims.list$beta.coy[,2])
  
  bs_topdown_skinny[[44]] # coy.t and wtd.t are independent given coy.tmin1, deerHarv.tmin1, wtd.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep44 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep44$sims.list$beta.coy[,2]) # 0.001066667
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  bs_topdown_skinny[[45]] # coy.t and moose.t are independent given coy.tmin1, moose.tmin1 and wolf.tmin1
  d_Sep45 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep45$sims.list$beta.coy[,2]) # 0
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  bs_topdown_skinny[[46]] # coy.t and elk.t are independent given coy.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep46 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep46$sims.list$beta.coy[,2]) # 6e-04
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  bs_topdown_skinny[[47]] # coy.t and wolf.t are independent given coy.tmin1, wolfHarv.tmin1 and wolf.tmin1
  d_Sep47 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep47$sims.list$beta.coy[,2]) # 0.005933333
  ############  ADD TO SEM.... BUT PROBABLY NOT - Not biologically reasonable  #############
  
  bs_topdown_skinny[[48]] # bearHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep48 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep48$sims.list$beta.harvest[,2]) # 0.1323333
  
  bs_topdown_skinny[[49]] # bearHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep49 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep49$sims.list$beta.harvest[,2]) # 0.05426667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[50]] # bearHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep50 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep50$sims.list$beta.harvest[,1]) # 0.1934667
  
  bs_topdown_skinny[[51]] # bearHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep51 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep51$sims.list$beta.harvest[,2]) # 0.08106667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[52]] # bearHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep52 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep52$sims.list$beta.harvest[,2]) # 0.1676
  
  bs_topdown_skinny[[53]] # bear.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
  d_Sep53 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep53$sims.list$beta.bear[,1]) # 0.1090667
  
  bs_topdown_skinny[[54]] # bear.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep54 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep54$sims.list$beta.bear[,1]) # 0.07173333
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[55]] # bear.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep55 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep55$sims.list$beta.bear[,1]) # 0.009
  ############  ADD TO SEM  #############
  
  #' bs_topdown_skinny[[56]] # bear.t and lionHarv.tmin1 are independent given bearHarv.tmin1 and bear.tmin1
  #' d_Sep56 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$lionHarv.tmin1, post = d_Sep56$sims.list$beta.bear[,2]) 
  #'
  #' bs_topdown_skinny[[57]] # bear.t and lion.tmin1 are independent given bearHarv.tmin1 and bear.tmin1
  #' d_Sep57 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$lion.tmin1, post = d_Sep57$sims.list$beta.bear[,2]) 
  
  bs_topdown_skinny[[58]] # bear.t and lion.t are independent given bearHarv.tmin1, bear.tmin1, lionHarv.tmin1 and lion.tmin1
  d_Sep58 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$lion.t_hat, post = d_Sep58$sims.list$beta.bear[,2]) # 0
  ############  ADD TO SEM  #############
  
  # bs_topdown_skinny[[59]] # bear.t and wolfHarv.tmin1 are independent given bearHarv.tmin1, bear.tmin1
  # d_Sep59 <- run_jags()
  # p.rope(y = data_JAGS_bundle_top$wolfHarv.tmin1, post = d_Sep59$sims.list$beta.bear[,2]) 
  # 
  # bs_topdown_skinny[[60]] # bear.t and wolf.tmin1 are independent given bearHarv.tmin1, bear.tmin1
  # d_Sep60 <- run_jags()
  # p.rope(y = data_JAGS_bundle_top$wolf.tmin1_hat, post = d_Sep60$sims.list$beta.bear[,2]) 
  
  bs_topdown_skinny[[61]] # bear.t and wtd.t are independent given bearHarv.tmin1, bear.tmin1, deerHarv.tmin1, wtd.tmin1, coy.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep61 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep61$sims.list$beta.bear[,2]) # 6.666667e-05
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[62]] # bear.t and moose.t are independent given bearHarv.tmin1, bear.tmin1, moose.tmin1 and wolf.tmin1
  d_Sep62 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep62$sims.list$beta.bear[,2]) # 0.01746667
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[63]] # bear.t and elk.t are independent given bearHarv.tmin1, bear.tmin1, elkHarv.tmin1, elk.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep63 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep63$sims.list$beta.bear[,2]) # 0.0002666667
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[64]] # bear.t and wolf.t are independent given bearHarv.tmin1, bear.tmin1, wolfHarv.tmin1 and wolf.tmin1
  d_Sep64 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep64$sims.list$beta.bear[,2]) # 0.001266667
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[65]] # lionHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep65 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep65$sims.list$beta.harvest[,2]) # 0.02313333
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[66]] # lionHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep66 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep66$sims.list$beta.harvest[,1]) # 0.09626667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[67]] # lionHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep67 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep67$sims.list$beta.harvest[,2]) # 0.09446667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  bs_topdown_skinny[[68]] # lionHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep68 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep68$sims.list$beta.harvest[,2]) # 0.01553333
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[69]] # lion.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep69 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep69$sims.list$beta.lion[,1]) # 0.1137333
  
  bs_topdown_skinny[[70]] # lion.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
  d_Sep70 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep70$sims.list$beta.lion[,1]) # 0.08226667
  ############  CONSIDER ADDING TO SEM - TRENDING  #############
  
  #' bs_topdown_skinny[[71]] # lion.t and wolfHarv.tmin1 are independent given lionHarv.tmin1 and lion.tmin1
  #' d_Sep71 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$wolfHarv.tmin1, post = d_Sep71$sims.list$beta.lion[,1]) 
  #' 
  #' bs_topdown_skinny[[72]] # lion.t and wolf.tmin1 are independent given lionHarv.tmin1 and lion.tmin1
  #' d_Sep72 <- run_jags()
  #' p.rope(y = data_JAGS_bundle_top$wolf.tmin1_hat, post = d_Sep72$sims.list$beta.lion[,1]) 
  
  bs_topdown_skinny[[73]] # lion.t and wtd.t are independent given lionHarv.tmin1, lion.tmin1, deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1 and wolf.tmin1
  d_Sep73 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep73$sims.list$beta.lion[,2]) # 0.002
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[74]] # lion.t and moose.t are independent given lionHarv.tmin1, lion.tmin1, moose.tmin1 and wolf.tmin1
  d_Sep74 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep74$sims.list$beta.lion[,2]) # 4e-04
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[75]] # lion.t and elk.t are independent given lionHarv.tmin1, lion.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1 and wolf.tmin1
  d_Sep75 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep75$sims.list$beta.lion[,2]) # 4e-04
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[76]] # lion.t and wolf.t are independent given lionHarv.tmin1, lion.tmin1, wolfHarv.tmin1 and wolf.tmin1
  d_Sep76 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep76$sims.list$beta.lion[,2]) # 0.0001333333
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[77]] # wolfHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep77 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wtd.t_hat, post = d_Sep77$sims.list$beta.harvest[,2]) # 0.3290667
  
  bs_topdown_skinny[[78]] # wolfHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  d_Sep78 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep78$sims.list$beta.harvest[,1]) # 0.4398667
  
  bs_topdown_skinny[[79]] # wolfHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
  d_Sep79 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep79$sims.list$beta.harvest[,2]) # 0.3474667
  
  bs_topdown_skinny[[80]] # wtd.t and moose.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, moose.tmin1 and wolf.tmin1
  d_Sep80 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$moose.t_hat, post = d_Sep80$sims.list$beta.wtd[,2]) # 0.0476
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[81]] # wtd.t and elk.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, elk.tmin1, elkHarv.tmin1, and wolf.tmin1
  d_Sep81 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep81$sims.list$beta.wtd[,2]) # 0.0066
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[82]] # wtd.t and wolf.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, wolf.tmin1 and wolfHarv.tmin1
  d_Sep82 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep82$sims.list$beta.wtd[,2]) # 0.02166667
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[83]] # moose.t and elk.t are independent given moose.tmin1, wolf.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1 and lion.tmin1
  d_Sep83 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$elk.t_hat, post = d_Sep83$sims.list$beta.moose[,2]) # 0.0003333333
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[84]] # moose.t and wolf.t are independent given moose.tmin1, wolf.tmin1 and wolfHarv.tmin1
  d_Sep84 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep84$sims.list$beta.moose[,2]) # 0
  ############  ADD TO SEM  #############
  
  bs_topdown_skinny[[85]] # elk.t and wolf.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1, wolf.tmin1 and wolfHarv.tmin1
  d_Sep85 <- run_jags()
  p.rope(y = data_JAGS_bundle_top$wolf.t_hat, post = d_Sep85$sims.list$beta.elk[,2]) # 0.0146
  ############  ADD TO SEM  #############
  
  