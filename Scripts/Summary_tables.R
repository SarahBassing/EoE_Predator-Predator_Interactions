  #'  ------------------------------
  #'  Summary tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  ------------------------------
  #'  Summarize deployment and detection data
  #'  ------------------------------
  
  #'  Load libraries
  library(tidyverse)
  
  #'  Load unique detection event data
  load("./Data/Detection_Histories/eoe20s_det_events.RData") 
  eoe20s_det_events <- mutate(eoe20s_det_events, Season = "Smr20")
  load("./Data/Detection_Histories/eoe21s_det_events.RData")
  eoe21s_det_events <- mutate(eoe21s_det_events, Season = "Smr21")
  
  #'  Merge data
  eoe_det_events <- rbind(eoe20s_det_events, eoe21s_det_events) %>%
    mutate(GMU = sub("_.*", "", NewLocationID)) %>%
    filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" | 
             Species == "mountain_lion" | Species == "wolf")

  #'  Detection events per species, GMU, and season
  eoe_det_table1 <- eoe_det_events %>%
    group_by(Season, GMU, Species) %>%
    summarise(n_detections = n()) %>%
    ungroup()
  
  #'  Detection events per species and season, across GMUs
  eoe_det_table2 <- eoe_det_events %>%
    group_by(Season, Species) %>%
    summarise(n_detections = n()) %>%
    ungroup()
  
  #'  Detection events per species, across seasons and GMUs
  eoe_det_table3 <- eoe_det_events %>%
    group_by(Species) %>%
    summarise(n_detections = n()) %>%
    ungroup()
  
  #'  Detection events per species and GMUs, across seasons 
  eoe_det_table4 <- eoe_det_events %>%
    group_by(Species, GMU) %>%
    summarise(n_detections = n()) %>%
    ungroup()
  
  write.csv(eoe_det_table1, "./Outputs/SummaryTable_EoESmr20-21a.csv")
  write.csv(eoe_det_table2, "./Outputs/SummaryTable_EoESmr20-21b.csv")
  write.csv(eoe_det_table3, "./Outputs/SummaryTable_EoESmr20-21c.csv")
  write.csv(eoe_det_table4, "./Outputs/SummaryTable_EoESmr20-21d.csv")
  
  