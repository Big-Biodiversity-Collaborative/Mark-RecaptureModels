# Summary of general results for thesis

# Gaby Samaniego
# gabysamaniego@arizona.edu
# 2025-06-11

# Load packages
library(tidyverse)

# Clear environment
rm(list = ls()) 

# Load banding data
dat.raw <- read.csv('output/capture-data/cleanded-capture-data-RMNP-full.csv')

# Identify birds banded in September
sept.birds <- dat.raw %>% 
  filter(band_status == 1 & month == 9) %>% 
  select(band)

# Prepare data set for survival analysis 
dat <- dat.raw %>%
  filter(!band %in% sept.birds$band,# exclude September birds
         !band_site %in% c('WB2','WB1', 'WPK1', 'NFPC', 'POLC', 'SHIP'),
         month != 9) %>% # exclude September recaptures
  select(band, band_status, year, sex, obssite, band_age, band_site) %>% 
  rename(age = band_age) %>% 
  distinct()

# Unique individuals in the data
length(unique(dat$band))
# 9507 individuals 

# How many unique individuals per banding site
sites <- dat %>%
  group_by(band_site) %>%
  summarize(unique_bands = n_distinct(band))

# How many birds were banded as juveniles? 
juveniles <- dat %>% 
  filter(band_status == 1,
         age == 'HY') %>% 
  select(band) %>% 
  distinct()
length(unique(juveniles$band))
# 1338 juvenile individuals

# How many individuals were banded as adults?
adults <- dat %>% 
  filter(band_status == 1,
         age == 'AHY') %>% 
  select(band) %>% 
  distinct()
length(unique(adults$band))
# 8169 adult individuals

# How many of the banded individuals are females?
females <- dat %>% 
  filter(band_status == 1,
         sex == 'F') %>% 
  select(band) %>% 
  distinct()
length(unique(females$band))  
# 6884 are females

# How many of the banded individuals are males?
males <- dat %>% 
  filter(band_status == 1,
         sex == 'M') %>% 
  select(band) %>% 
  distinct()
length(unique(males$band))  
# 2623 are males

# Load covariate data 
# Load data
winter.mx <- read.csv('output/weather-data/covariates-output/winter-covar-mexico.csv')
summer.co <- read.csv('output/weather-data/covariates-output/summer-covar-colorado.csv')
summer.co.swe <- read.csv('output/weather-data/covariates-output/winter-swe-colorado.csv')

# Look at basic stats
summary(winter.mx)
summary(summer.co)
summary(summer.co.swe)

mean(summer.co.swe$total_swe_winter_co)
