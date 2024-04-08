  #'  --------------------------------
  #'  Structural Equation Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2024
  #'  --------------------------------
  #'  Source data formatting script and run structural equation models (SEM) to
  #'  test hypotheses about how predator-prey and predator-predator interactions
  #'  influence wildlife populations in northern Idaho.
  #'  --------------------------------
  
  library(piecewiseSEM)
  library(labelled)
  library(lme4)
  library(tidyverse)

  #'  Run script to prepare data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_data_for_SEMs.R")

  #'  Take a quick look
  head(localN_z); head(localN_z_1YrLag)
  
  #'  ---------------------------
  ####  SEM direct effects only  ####
  #'  ---------------------------
  #'  For simplicity, writing models to only include interference competition 
  #'  (predators directly interact) in addition to the bottom-up vs top-down relationships.
  #'  But if d-sep tests indicate no direct interaction btwn competitors this would
  #'  support exploitation hypothesis that all competitive interactions occur through
  #'  use of shared prey. Effects of another species at T-1 hypothesized to refelct
  #'  demographic effect of species interactions, effects of another species at t
  #'  hypothesized to reflect behavioral response of interacting species.
  
  #'  Top down, predator affects primary prey
  h.td1 <- psem(
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1 + wolf.T, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T + bear_black.Tminus1 + bear_black.T, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + mountain_lion.Tminus1 + mountain_lion.T + bear_black.Tminus1 + bear_black.T, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bobcat.Tminus1 + bobcat.T + coyote.Tminus1 + coyote.T, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1 + wolf.T, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + coyote.Tminus1 + coyote.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + wolf.T + bear_black.Tminus1 + bear_black.T, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(h.td1)
  
  #'  Updated h.td1 model based on results from d-sep tests
  h.td1.a <- psem(
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + wolf.T + bear_black.T + bear_black.Tminus1, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.T + bear_black.T + coyote.Tminus1 + bobcat.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bobcat.Tminus1 + bobcat.T, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1 + whitetailed_deer.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + mountain_lion.Tminus1 + coyote.Tminus1 + coyote.T + whitetailed_deer.Tminus1 + lagomorphs.Tminus1, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + moose.T, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + bobcat.Tminus1 + wolf.T + lagomorphs.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + elk.Tminus1, weights = precision_coyote.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  h.td1.a <- update(h.td1.a, moose.T %~~% bobcat.Tminus1)
  h.td1.a <- update(h.td1.a, moose.T %~~% coyote.Tminus1)
  h.td1.a <- update(h.td1.a, elk.T %~~% moose.Tminus1)
  h.td1.a <- update(h.td1.a, elk.T %~~% moose.T)
  h.td1.a <- update(h.td1.a, elk.T %~~% bobcat.Tminus1)
  h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% moose.Tminus1)
  h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% moose.T)
  h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% elk.Tminus1)
  h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% elk.T)
  h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% lagomorphs.Tminus1)
  # h.td1.a <- update(h.td1.a, whitetailed_deer.T %~~% lagomorphs.T)
  h.td1.a <- update(h.td1.a, lagomorphs.T %~~% bear_black.Tminus1)
  h.td1.a <- update(h.td1.a, lagomorphs.T %~~% bear_black.T)
  h.td1.a <- update(h.td1.a, lagomorphs.T %~~% mountain_lion.T)
  # h.td1.a <- update(h.td1.a, lagomorphs.T %~~% whitetailed_deer.T)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% moose.Tminus1)
  # h.td1.a <- update(h.td1.a, bobcat.T %~~% elk.Tminus1)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% elk.T)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% bear_black.Tminus1)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% bear_black.T)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% wolf.Tminus1)
  h.td1.a <- update(h.td1.a, bobcat.T %~~% wolf.T)
  h.td1.a <- update(h.td1.a, coyote.T %~~% moose.T)
  h.td1.a <- update(h.td1.a, mountain_lion.T %~~% bobcat.Tminus1)
  h.td1.a <- update(h.td1.a, bear_black.T %~~% lagomorphs.Tminus1)
  h.td1.a <- update(h.td1.a, bear_black.T %~~% coyote.Tminus1)
  summary(h.td1.a)
  
  
  #'  Top down most abundant prey
  h.td2 <- psem(
    lm(moose.T ~ moose.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T + bear_black.Tminus1 + bear_black.T, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1 + wolf.T, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + coyote.Tminus1 + coyote.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + wolf.T + bear_black.Tminus1 + bear_black.T, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(h.td2)
  
  #'  Bottom up prey affects primary predator
  h.bu1 <- psem(
    lm(moose.T ~ moose.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + moose.T + elk.Tminus1 + elk.T + mountain_lion.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + elk.T + wolf.Tminus1 + wolf.T + PercDisturbedForest.Tminus1 + PercDisturbedForest.T, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + lagomorphs.Tminus1 + lagomorphs.T + coyote.Tminus1 + coyote.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + elk.T + whitetailed_deer.Tminus1 + whitetailed_deer.T + wolf.Tminus1 + wolf.T + bear_black.Tminus1 + bear_black.T, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(h.bu1)
  
  #'  Bottom up most abundant prey
  h.bu2 <- psem(
    lm(moose.T ~ moose.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + PercDisturbedForest.Tminus1 + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + mountain_lion.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + wolf.Tminus1 + wolf.T + PercDisturbedForest.Tminus1 + PercDisturbedForest.T, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + coyote.Tminus1 + coyote.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + wolf.Tminus1 + wolf.T + bear_black.Tminus1 + bear_black.T, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + wolf.Tminus1 + wolf.T + mountain_lion.Tminus1 + mountain_lion.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(h.bu2)
  
  
  
  
  
  
  
  #'  ---------------------------------------------
  ####  SEM based on hypothesized causal networks  ####
  #'  ---------------------------------------------
  #'  Models build around several broad ecological hypotheses
  #'  1) Top-down vs bottom-up system
  #'  - food web is structured through top-down factors (predator driven)
  #'  - food web is structured through bottom-up factors (habitat driven)
  #'  - food web is structured through both top-down (predator) and bottom-up (habitat) factors
  #'  2) Primary prey vs prey availability shape predator-prey relationships
  #'  - predator-prey interactions center on predators and their primary prey species
  #'  - predator-prey interactions center on predators and the most abundant prey species
  #'  3) Interference vs exploitative competition shape predator-predator relationships
  #'  - predators affect competitors through interference competition
  #'  - predators affect competitors through exploitative competition via their
  #'      primary prey or the most abundant prey species
  #'  ----------------------------------------------
   
  #####  Top-down system  #####
  #'  Models center around hypothesis that predators structure the food web. Models 
  #'  include a combination of possible species interactions within the top-down framework.
  #'  ----------------------------------------------
  ######  H.td1  ######
  #'  Predators affect their primary prey species; predators affect their competitors
  #'  through interference competition
  H.td1_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z, weights = precision_wolf.yr2),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z, weights = precision_wolf.yr3),
    lm(mountain_lion.yr2 ~ wolf.yr1 + bear_black.yr1 + mountain_lion.yr1, data = localN_z, weights = precision_mountain_lion.yr2),
    lm(mountain_lion.yr3 ~ wolf.yr2 + bear_black.yr2 + mountain_lion.yr2, data = localN_z, weights = precision_mountain_lion.yr3),
    lm(bear_black.yr2 ~ wolf.yr1 + bear_black.yr1, data = localN_z, weights = precision_bear_black.yr2),
    lm(bear_black.yr3 ~ wolf.yr2 + bear_black.yr2, data = localN_z, weights = precision_bear_black.yr3),
    lm(coyote.yr2 ~ wolf.yr1 + mountain_lion.yr1 + coyote.yr1, data = localN_z, weights = precision_coyote.yr2),
    lm(coyote.yr3 ~ wolf.yr2 + mountain_lion.yr2 + coyote.yr2, data = localN_z, weights = precision_coyote.yr3),
    lm(bobcat.yr2 ~ coyote.yr1 + bobcat.yr1, data = localN_z, weights = precision_bobcat.yr2),
    lm(bobcat.yr3 ~ coyote.yr2 + bobcat.yr2, data = localN_z, weights = precision_bobcat.yr3),
    lm(moose.yr2 ~ wolf.yr1 + moose.yr1, data = localN_z, weights = precision_moose.yr2),
    lm(moose.yr3 ~ wolf.yr2 + moose.yr2, data = localN_z, weights = precision_moose.yr3),
    lm(elk.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + elk.yr1, data = localN_z, weights = precision_elk.yr2),
    lm(elk.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + elk.yr2, data = localN_z, weights = precision_elk.yr3),
    lm(whitetailed_deer.yr2 ~ mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z, weights = precision_whitetailed_deer.yr2),
    lm(whitetailed_deer.yr3 ~ mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z, weights = precision_whitetailed_deer.yr3),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z, weights = precision_lagomorphs.yr2),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z, weights = precision_lagomorphs.yr3),
    data = localN_z)
  # basisSet(H.td1_psem)
  # dSep(H.td1_psem)
  summary(H.td1_psem)
  
  ######  H.td2  ######
  #'  Predators affect their primary prey species; predators affect their competitors
  #'  through exploitative competition via their primary prey
  H.td2_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ elk.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ elk.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ elk.yr1 + bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ elk.yr2 + bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(moose.yr2 ~ wolf.yr1 + moose.yr1, data = localN_z),
    lm(moose.yr3 ~ wolf.yr2 + moose.yr2, data = localN_z),
    lm(elk.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + elk.yr1, data = localN_z),
    lm(elk.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + elk.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td2_psem)
  
  ######  H.td3  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through interference competition
  H.td3_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ wolf.yr1 + bear_black.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ wolf.yr2 + bear_black.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ wolf.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ wolf.yr2 + coyote.yr2, data = localN_z),
    lm(bobcat.yr2 ~ coyote.yr1 + bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ coyote.yr2 + bobcat.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    data = localN_z)
  summary(H.td3_psem)
  
  ######  H.td4  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through exploitative competition via the most abundant prey (wolves dominant)
  H.td4_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ bear_black.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ bear_black.yr2 + mountain_lion.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td4_psem)
  
  ######  H.td5  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through exploitative competition via the most abundant prey (lions dominant)
  H.td5_psem <- psem(
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(wolf.yr2 ~ whitetailed_deer.yr1 + wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ whitetailed_deer.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td5_psem)
  
  #'  ----------------------------------------------
  #####  Bottom-up system  #####
  #'  Models center around hypothesis that prey structure the food web. Models include
  #'  a combination of possible species interactions within the bottom-up framework.
  #'  ----------------------------------------------
  ######  H.bu1  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by their primary prey species; predators affect their 
  #'  competitors through interference competition
  H.bu1_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1 + bear_black.yr1 + wolf.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2 + bear_black.yr2 + wolf.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1 + wolf.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + wolf.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu1_psem)
  
  ######  H.bu2  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by their primary prey species; predators affect their 
  #'  competitors through exploitative competition
  H.bu2_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu2_psem)
  
  ######  H.bu3  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through interference competition
  H.bu3_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu3_psem)
  
  ######  H.bu4  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (wolves dominant)
  H.bu4_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu4_psem)
  
  ######  H.bu5  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (lions dominant)
  H.bu5_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + mountain_lion.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + mountain_lion.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu5_psem)
  
  #'  ----------------------------------------------
  #####  Top-down and bottom-up system  #####
  #'  Models center around hypothesis that prey structure the food web. Models include
  #'  a combination of possible species interactions within the bottom-up framework.
  #'  ----------------------------------------------
  ######  H.tdbu1  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect their primary prey species; predators affect their 
  #'  competitors through interference competition
  H.tdbu1_psem <- psem(        # + DecFeb_WSI.yr2
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + whitetailed_deer.yr1, data = localN_z), # + bear_black.yr1 + mountain_lion.yr1
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + mountain_lion.yr2, data = localN_z), # + bear_black.yr2 + whitetailed_deer.yr2
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr2 + moose.yr1 + whitetailed_deer.yr1, data = localN_z), # + bear_black.yr1 + wolf.yr1 + moose.yr2 + whitetailed_deer.yr2
    lm(mountain_lion.yr3 ~ wolf.yr3, data = localN_z), # + bear_black.yr2 + wolf.yr2 + mountain_lion.yr2 + 
    lm(bear_black.yr2 ~ bear_black.yr1 + whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z), # + wolf.yr1
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z), # + whitetailed_deer.yr2
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + elk.yr1 + bear_black.yr2 + bear_black.yr1, data = localN_z), # + wolf.yr1 + mountain_lion.yr1 + wolf.yr2
    lm(coyote.yr3 ~ coyote.yr2 + bear_black.yr3 + wolf.yr3, data = localN_z), # + wolf.yr2 + mountain_lion.yr2 + whitetailed_deer.yr2 + elk.yr2 + bear_black.yr2
    lm(bobcat.yr2 ~ bobcat.yr1 + coyote.yr2 + mountain_lion.yr2 + whitetailed_deer.yr1, data = localN_z), # + coyote.yr1 + mountain_lion.yr1 + lagomorphs.yr1
    lm(bobcat.yr3 ~ bobcat.yr2 + coyote.yr3 + mountain_lion.yr2 + mountain_lion.yr3 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z), # + coyote.yr2
    lm(moose.yr2 ~ moose.yr1 + wolf.yr2, data = localN_z), # + wolf.yr1 + habitat_class.yr1
    lm(moose.yr3 ~ moose.yr2 + wolf.yr3, data = localN_z), # + wolf.yr2 + habitat_class.yr2
    lm(elk.yr2 ~ elk.yr1 + bear_black.yr2 + moose.yr1 + moose.yr2 + habitat_class.yr1, data = localN_z), # + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1
    lm(elk.yr3 ~ elk.yr2 + bear_black.yr3 + moose.yr2 + moose.yr3 + habitat_class.yr2, data = localN_z), # + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + mountain_lion.yr1+ coyote.yr2, data = localN_z), # + habitat_class.yr1 + bear_black.yr1 + wolf.yr1 + moose.yr2 
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2 + moose.yr3 + wolf.yr2 + coyote.yr3, data = localN_z), # + mountain_lion.yr2 + bear_black.yr2
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + habitat_class.yr1, data = localN_z), # + coyote.yr1 + bobcat.yr1
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + habitat_class.yr2, data = localN_z), #  + coyote.yr2 + bobcat.yr2
    data = localN_z)
  
  #'  Add correlated errors 
  #'  Accounting for autoregression at different time lags than t to t+1
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% mountain_lion.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr3 %~~% bear_black.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% bobcat.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, moose.yr3 %~~% moose.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, whitetailed_deer.yr3 %~~% whitetailed_deer.yr1)
  
  #'  Among species already hypothesized to interact but at a different time lag
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% mountain_lion.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, coyote.yr3 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, moose.yr3 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr2 %~~% habitat_class.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr3 %~~% habitat_class.yr3)
  # H.tdbu1_psem <- update(H.tdbu1_psem, elk.yr3 %~~% habitat_class.yr1)
  
  #'  Among species/habitats not previously hypothesize to interact but interactions are biologically plausible
  #'  Predator abundances are likely correlated with certain habitat types related
  #'  to their own preferences or prey distributions
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% habitat_class.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% elk.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr2 %~~% habitat_class.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr2 %~~% habitat_class.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% habitat_class.yr3)
  
  #'  Among species not previously hypothesized to interact and not sure what's driving this correlation structure
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% bobcat.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% bobcat.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% wolf.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, lagomorphs.yr2 %~~% moose.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, lagomorphs.yr3 %~~% moose.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, whitetailed_deer.yr3 %~~% elk.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, elk.yr2 %~~% whitetailed_deer.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, coyote.yr3 %~~% mountain_lion.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% bear_black.yr1)
  
  #'  Add causal relationships identified by Dsep that were missing in original model 
  #'  Among species already hypothesized to interact but at a different time lag
  #'  Bear t ~ habitat t --> nixed due to convergence issues (probably too many categories) so added correlated error
  #'  Moose t ~ wolf t
  #'  Bobcat t ~ coyote t + mountain_lion t
  #'  Elk t ~ bear_black t + moose t
  #'  Coyote t ~ bear_balck t + wolf t
  #'  Mountain lion t ~ wolf t
  #'  Whitetailed_deer t ~ moose t + coyote t
  
  #'  Among species not previously hypothesized to interact but interactions are biologically plausible
  #'  Bobcat t ~ Mountain lion (t-1) + whitetailed_deer (t-1) + lagomorphs (t-1)
  #'  Bear t ~ whitetailed_deer (t-1)
  #'  Coyote t ~ whitetailed_deer (t-1) + elk (t-1) + bear (t-1)
  #'  Mountain lion t ~ moose (t-1) + whitetailed_deer (t-1)
  #'  wolf t ~ moose (t-1) + whitetailed_deer (t-1) + mountain_lion (t-1) + bear (t-1) 
  #'  Elk t ~ moose (t-1)
  #'  Whitetailed_deer t ~ wolf (t-1)
  
  out <- summary(H.tdbu1_psem, getOption("max.print")); print(out)
  out_coeffs <- out$coefficients
  out_ds <- dSep(H.tdbu1_psem)
  
  ######  H.tdbu2  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect their primary prey species; predators affect their 
  #'  competitors through exploitative competition
  H.tdbu2_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu2_psem)
  
  ######  H.tdbu3  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect most abundant prey species; predators affect their 
  #'  competitors through interference competition
  H.tdbu3_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + coyote.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu3_psem)
  
  ######  H.tdbu4  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (wolves dominant)
  H.tdbu4_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu4_psem)
  
  ######  H.tdbu5  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (lions dominant)
  H.tdbu5_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu5_psem)