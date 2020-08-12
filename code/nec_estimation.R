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

write('
data {
  int <lower = 0> N; // sample size    
  vector[N] x; // concentration                 
  vector[N] y; // response (death rate)            
  real<lower = 0> y[N, 2];    
}

parameters {
  real <lower = 0> a; // basal response    
  real <lower = 0> b; // rate of decay of the response  
  real <lower = 0> g; // threshold of NEC   
}

model {
  a ~ gamma(0.0001, 0.0001) // arguments are alpha, beta
  b ~ gamma(0.0001, 0.0001)
  g ~ gamma(0.0001, 0.0001)
  
  y ~ binomial 
}
generated quantities {
  }
}', 
here('./code/stan_files/nec_estimation.stan'))

