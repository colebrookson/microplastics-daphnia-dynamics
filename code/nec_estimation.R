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

### read in data
dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data' #private data repo

offspring = 
  read_csv(
    'C:/Users/brookson/Documents/Github/Schur-etal-Data/offspring_data.csv',
    guess_max = 20000)
offspring = offspring %>% 
  select(-`X1`)

### create data groupings
str(offspring)
offspring$generation = as.factor(offspring$generation)
levels(offspring$generation)[levels(offspring$generation)== 'F3'] = 3

for(i in unique(offspring$generation)) {
  
  temp = offspring[which(offspring$generation == i),]
  assign(paste('generation_f', i, sep = ''), temp)
  rm(temp)
} 

### pull the data and get into proper format
#note the response is the proportion of surviving/original
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

#put into a list
n = nrow(nec_fit_data_f0)
x = nec_fit_data_f0$x
y = nec_fit_data_f0$y
nec_fit_data_f0_sample = list(n=n,
                              x=x,
                              y=y)

nec_model = stan_model(file = here('./code/stan_files/nec_estimation.stan'))

nec_fit = sampling(nec_model, 
                   data = nec_fit_data_f0_sample, 
                   seed = 12); beep(3)
