  #'  ---------------------------------
  #'  Model selection & result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ---------------------------------
  #'  Script to identify most supported model using DIC and format result tables.
  #'  ---------------------------------
  
  #'  Clean workspace & load libraries
  
  library(jagsUI)
  library(tidyverse)
  
  #'  Load model outputs
  #'  Wolf-Bear
  #'  Wolf-Coyote
  #'  Wolf-Lion
  #'  Lion-Bear
  #'  Lion-Bobcat
  #'  Coyote-Bobcat