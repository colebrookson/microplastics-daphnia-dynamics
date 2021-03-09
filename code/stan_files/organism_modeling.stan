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
    real l = z[1];
    real c = z[2];

    real cstar = theta[1];  
    real cq = theta[2];
    real NEC = theta[3];
    real ke = theta[4];


    real LFC_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    real LFC_cq = ke*(c-cq);
    real PS400_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    real PS400_cq = ke*(c-cq);
    real PS2000_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    real PS2000_cq = ke*(c-cq);
    real PS10000_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    real PS10000_cq = ke*(c-cq);
    return {LFC_ll, LFC_cq, PS400_ll, PS400_cq,
            PS2000_ll, PS2000_cq, PS10000_ll, PS10000_cq};
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  real ts[N]; // time points
  vector[N] y; // this is the length data
  vector[N] x; // this is the concentration data
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

