//////////////////// 
////////////////////
// This code contains the stan code for the individual component as a
// part of the larger analysis of the effect of microplastics on daphnia
////////////////////
////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2020-11-18↕
////////////////////
////////////////////

// define the functions first 
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real b = theta[3];
    real c = theta[4];
    real u = theta[5];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*0.075*P))-u);
    return { dP_dt, dH_dt };
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

