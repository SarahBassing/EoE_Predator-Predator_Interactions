  #'  --------------------------------
  #'  Reformat 50 years of PRISM data
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2024
  #'  --------------------------------
  #'  Load 50 years of monthly PRISM weather data, reformat and generate 
  #'  standardized monthly averages per winter from Dec. 1972 - Feb. 2023.
  #'  Standardizing across 50 years of data means each monthly average is 
  #'  represented relative to the 50-year average and 50-year variability.
  #'  Winter months = December, January, February
  #'  --------------------------------
  
  #'  Load libraries
  library(tidyverse)
  
  #' #'  PRISM total monthly precipitation (averaged per GMU) per camera site
  #' precip <- read_csv("./Data/GEE outputs/PRISM_camBuff_monthly_total_precip_1972_2023.csv") %>%  
  #'   dplyr::select(c(NewLocationID, date, meanMonthlyValue))
  #' #'  PRISM minimum monthly temperature (averaged per GMU) per camera site
  #' temp <- read_csv("./Data/GEE outputs/PRISM_camBuff_monthly_min_temp_1972_2023.csv") %>%  
  #'   dplyr::select(c(NewLocationID, date, meanMonthlyValue))
  
  #'  PRISM total monthly precipitation (averaged per GMU) per cluster
  precip <- read_csv("./Data/GEE outputs/PRISM_ClusterAvg_monthly_total_precip_1972_2023.csv") %>%  
    dplyr::select(c(featureID, Clusters, date, meanMonthlyValue)) %>%
    mutate(GMU = ifelse(str_detect(featureID, "1_1_"), "GMU1", featureID),
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU10A", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU6", GMU))
  #'  PRISM minimum monthly temperature (averaged per GMU) per cluster
  temp <- read_csv("./Data/GEE outputs/PRISM_ClusterAvg_monthly_min_temp_1972_2023.csv") %>%  
    dplyr::select(c(featureID, Clusters, date, meanMonthlyValue)) %>%
    mutate(GMU = ifelse(str_detect(featureID, "1_1_"), "GMU1", featureID),
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU10A", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU6", GMU))
  
  #'  Function to add season column to weather data
  add_season <- function(prism) {
    annual_winter_season <- prism %>%
      #'  Reformat
      mutate(Year = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%Y")),
             Month = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%m"))) %>%
      #'  Filter to winter months (Dec, Jan, Feb) of each year
      filter(Month <= 2 | Month == 12) %>%
      #'  Label monthly data with a yearly seasonal classification (e.g., Dec. '72 - Feb. '73 = Wtr7273)
      #'  Is there a more efficient way to do this? Probably.
      #'  Am I going to spend the time to figure it out? Nope.
      mutate(Season = ifelse(date <= "1973-02-01", "Wtr7273", "Wtr2223"),
             Season = ifelse(date >= "1973-12-01" & date <= "1974-02-01", "Wtr7374", Season),
             Season = ifelse(date >= "1974-12-01" & date <= "1975-02-01", "Wtr7475", Season),
             Season = ifelse(date >= "1975-12-01" & date <= "1976-02-01", "Wtr7576", Season),
             Season = ifelse(date >= "1976-12-01" & date <= "1977-02-01", "Wtr7677", Season),
             Season = ifelse(date >= "1977-12-01" & date <= "1978-02-01", "Wtr7778", Season),
             Season = ifelse(date >= "1978-12-01" & date <= "1979-02-01", "Wtr7879", Season),
             Season = ifelse(date >= "1979-12-01" & date <= "1980-02-01", "Wtr7980", Season),
             Season = ifelse(date >= "1980-12-01" & date <= "1981-02-01", "Wtr8081", Season),
             Season = ifelse(date >= "1981-12-01" & date <= "1982-02-01", "Wtr8182", Season),
             Season = ifelse(date >= "1982-12-01" & date <= "1983-02-01", "Wtr8283", Season),
             Season = ifelse(date >= "1983-12-01" & date <= "1984-02-01", "Wtr8384", Season),
             Season = ifelse(date >= "1984-12-01" & date <= "1985-02-01", "Wtr8485", Season),
             Season = ifelse(date >= "1985-12-01" & date <= "1986-02-01", "Wtr8586", Season),
             Season = ifelse(date >= "1986-12-01" & date <= "1987-02-01", "Wtr8687", Season),
             Season = ifelse(date >= "1987-12-01" & date <= "1988-02-01", "Wtr8788", Season),
             Season = ifelse(date >= "1988-12-01" & date <= "1989-02-01", "Wtr8889", Season),
             Season = ifelse(date >= "1989-12-01" & date <= "1990-02-01", "Wtr8990", Season),
             Season = ifelse(date >= "1990-12-01" & date <= "1991-02-01", "Wtr9091", Season),
             Season = ifelse(date >= "1991-12-01" & date <= "1992-02-01", "Wtr9192", Season),
             Season = ifelse(date >= "1992-12-01" & date <= "1993-02-01", "Wtr9293", Season),
             Season = ifelse(date >= "1993-12-01" & date <= "1994-02-01", "Wtr9394", Season),
             Season = ifelse(date >= "1994-12-01" & date <= "1995-02-01", "Wtr9495", Season),
             Season = ifelse(date >= "1995-12-01" & date <= "1996-02-01", "Wtr9596", Season),
             Season = ifelse(date >= "1996-12-01" & date <= "1997-02-01", "Wtr9697", Season),
             Season = ifelse(date >= "1997-12-01" & date <= "1998-02-01", "Wtr9798", Season),
             Season = ifelse(date >= "1998-12-01" & date <= "1999-02-01", "Wtr9899", Season),
             Season = ifelse(date >= "1999-12-01" & date <= "2000-02-01", "Wtr9900", Season),
             Season = ifelse(date >= "2000-12-01" & date <= "2001-02-01", "Wtr0001", Season),
             Season = ifelse(date >= "2001-12-01" & date <= "2002-02-01", "Wtr0102", Season),
             Season = ifelse(date >= "2002-12-01" & date <= "2003-02-01", "Wtr0203", Season),
             Season = ifelse(date >= "2003-12-01" & date <= "2004-02-01", "Wtr0304", Season),
             Season = ifelse(date >= "2004-12-01" & date <= "2005-02-01", "Wtr0405", Season),
             Season = ifelse(date >= "2005-12-01" & date <= "2006-02-01", "Wtr0506", Season),
             Season = ifelse(date >= "2006-12-01" & date <= "2007-02-01", "Wtr0607", Season),
             Season = ifelse(date >= "2007-12-01" & date <= "2008-02-01", "Wtr0708", Season),
             Season = ifelse(date >= "2008-12-01" & date <= "2009-02-01", "Wtr0809", Season),
             Season = ifelse(date >= "2009-12-01" & date <= "2010-02-01", "Wtr0910", Season),
             Season = ifelse(date >= "2010-12-01" & date <= "2011-02-01", "Wtr1011", Season),
             Season = ifelse(date >= "2011-12-01" & date <= "2012-02-01", "Wtr1112", Season),
             Season = ifelse(date >= "2012-12-01" & date <= "2013-02-01", "Wtr1213", Season),
             Season = ifelse(date >= "2013-12-01" & date <= "2014-02-01", "Wtr1314", Season),
             Season = ifelse(date >= "2014-12-01" & date <= "2015-02-01", "Wtr1415", Season),
             Season = ifelse(date >= "2015-12-01" & date <= "2016-02-01", "Wtr1516", Season),
             Season = ifelse(date >= "2016-12-01" & date <= "2017-02-01", "Wtr1617", Season),
             Season = ifelse(date >= "2017-12-01" & date <= "2018-02-01", "Wtr1718", Season),
             Season = ifelse(date >= "2018-12-01" & date <= "2019-02-01", "Wtr1819", Season),
             Season = ifelse(date >= "2019-12-01" & date <= "2020-02-01", "Wtr1920", Season),
             Season = ifelse(date >= "2020-12-01" & date <= "2021-02-01", "Wtr2021", Season),
             Season = ifelse(date >= "2021-12-01" & date <= "2022-02-01", "Wtr2122", Season),
             Season = ifelse(date >= "2022-12-01" & date <= "2023-02-01", "Wtr2223", Season))
        
    return(annual_winter_season)
  }
  precip_season <- add_season(precip)
  temp_season <- add_season(temp)
  
  #'  Function to reformat and summarize weather data
  format_weather <- function(prism) {
    avg_winter_weather <- prism %>%
      #'  Reformat columns
      transmute(Clusters = Clusters, #NewLocationID = NewLocationID,
                GMU = GMU,
                Season = Season,
                Date = date,
                meanWeather = meanMonthlyValue) %>%
      #'  Average monthly weather by year and site
      # group_by(NewLocationID, Season) %>%
      group_by(GMU, Clusters, Season) %>%
      summarise(DecFeb_meanWeather = mean(meanWeather),
                DecFeb_meanWeather_se = sd(meanWeather)/sqrt(nrow(.))) %>%
      ungroup() %>%
      #'  Standardize average monthly weather (mean = 0, SD = 1)
      #'  This represents the monthly mean relative to the 50 year average and SD
      mutate(DecFeb_meanWeather_z = as.numeric(scale(DecFeb_meanWeather)))
    
    #'  Double check standardized data
    print(mean(avg_winter_weather$DecFeb_meanWeather_z))
    print(sd(avg_winter_weather$DecFeb_meanWeather_z))
    
    return(avg_winter_weather)
  }
  wtr_totalPrecip <- format_weather(precip_season) %>%
    rename("DecFeb_meanPPT_mm" = "DecFeb_meanWeather") %>%
    rename("DecFeb_meanPPT_se" = "DecFeb_meanWeather_se") %>%
    rename("DecFeb_meanPPT_z" = "DecFeb_meanWeather_z")
  wtr_minTemp <- format_weather(temp_season) %>%
    rename("DecFeb_meanMinTemp_C" = "DecFeb_meanWeather") %>%
    rename("DecFeb_meanMinTemp_se" = "DecFeb_meanWeather_se") %>%
    rename("DecFeb_meanMinTemp_z" = "DecFeb_meanWeather_z") 
  