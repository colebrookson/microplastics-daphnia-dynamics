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


