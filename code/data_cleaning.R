########## 
##########
# This code contains some data cleaning and data preparation for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2020-03-28
##########
##########

# set-up =======================================================================

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)

dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data' #private data repo

offspring = 
  read_csv(
    'C:/Users/brookson/Documents/Github/Schur-etal-Data/offspring_data.csv',
    guess_max = 20000)
offspring = offspring %>% 
  select(-`X1`)

# Create data groupings ========================================================
str(offspring)
offspring$generation = as.factor(offspring$generation)
levels(offspring$generation)[levels(offspring$generation)== 'F3'] = 3

for(i in unique(offspring$generation)) {
  
  temp = offspring[which(offspring$generation == i),]
  assign(paste('generation_f', i, sep = ''), temp)
  rm(temp)
} 

# Put data in proper format ====================================================
nec_fit_data_f0 = generation_f0 %>% 
  filter(treatment %in% c('HFC', 'PS_400', 'PS_2000', 'PS_10000'),
         day == 21) 

nec_fit_data_f0$treatment = as.factor(nec_fit_data_f0$treatment)
levels(nec_fit_data_f0$treatment)
levels(nec_fit_data_f0$treatment) = c('0', '10000', '2000', '400')
nec_fit_data_f0$treatment = 
  as.numeric(levels(nec_fit_data_f0$treatment))[nec_fit_data_f0$treatment]

nec_fit_data_f0 = nec_fit_data_f0 %>% 
  select(treatment, survival) %>% 
  rename(x = treatment, y = survival)

