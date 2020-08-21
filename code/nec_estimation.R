########## 
##########
# This code contains the analysis for the fitting of the NEC parameter as a
# part of the larger analysis of the effect of microplastics on daphnia
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2020-07-04
##########
##########

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)

#read in data
dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data' #private data repo

offspring = 
  read_csv(
    'C:/Users/brookson/Documents/Github/Schur-etal-Data/offspring_data.csv',
    guess_max = 20000)
offspring = offspring %>% 
  select(-`X1`)

#create data groupings
str(offspring)
offspring$generation = as.factor(offspring$generation)
levels(offspring$generation)[levels(offspring$generation)== 'F3'] = 3

for(i in unique(offspring$generation)) {
  
  temp = offspring[which(offspring$generation == i),]
  assign(paste('generation_f', i, sep = ''), temp)
  rm(temp)
} 

nec_model = stan_model(file = here('./code/stan_files/nec_estimation.stan'))

nec_fit = sampling(nec_model, 
                   data = Data, 
                   seed = 12); beep(3)
