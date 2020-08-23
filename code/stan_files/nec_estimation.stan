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
  int <lower = 0> N; // sample size    
  vector[N] x; // concentration                 
  int <lower = 0> y[N]; // response (death rate) 
  real n[N]; // Nber of individuals
 
}

parameters {
  real <lower = 0> a; // basal response    
  real <lower = 0> b; // rate of decay of the response  
  real <lower = 0> g; // threshold of NEC  
  real <lower = 0, upper = 1> theta;
  
}

transformed parameters {
    //declare the response, proportion of surviving/original
  vector[N] yhat;
  
  //define the mean
  for(i in 1:N) {
  //define model in loop to look at yhat, int_step works like step in jags
  yhat[i] = a * exp(-b * (x[i] - g) * int_step(x[i] - g)); //
  }
}

model {

  
  
  //priors
  a ~ gamma(0.0001, 0.0001); // arguments are alpha, beta
  b ~ gamma(0.0001, 0.0001);
  g ~ gamma(0.0001, 0.0001);
  theta ~ beta(1,1);
  

  //define model outside loop to get estimate of y
  y ~ binomial_logit(n, yhat);
  
  
}

