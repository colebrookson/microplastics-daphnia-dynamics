//////////////////// 
////////////////////
// This code contains the analysis for the fitting of the NEC parameter as a
// part of the larger analysis of the effect of microplastics on daphnia
////////////////////
////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2020-07-04
////////////////////
////////////////////

data {
  int <lower = 0> n; // sample size    
  vector[n] x; // concentration                 
  int y[n]; // response (death rate)            
 
}

parameters {
  real <lower = 0> a; // basal response    
  real <lower = 0> b; // rate of decay of the response  
  real <lower = 0> g; // threshold of NEC  
  real <lower = 0, upper = 1> theta;
  
}

model {
  //declare the response, proportion of surviving/original
  vector[n] yhat;
  
  
  //priors
  a ~ gamma(0.0001, 0.0001); // arguments are alpha, beta
  b ~ gamma(0.0001, 0.0001);
  g ~ gamma(0.0001, 0.0001);
  theta ~ beta(1,1);
  
  //define the mean
  for(i in 1:n) {
  //define likilihood/probability model
  yhat[i] = a * exp(-b * (x[i] - g) * int_step(x[i] - g));
  
  }
  
  //define likilihood/probability model
  y ~ bernoulli_logit(yhat);
  
}

