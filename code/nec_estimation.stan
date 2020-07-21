/////////////////////////
// This code contains the Stan file for fitting NEC values for an analysis
// investigating the effects of microplastic pollution on Daphnia magna
// presented in Brookson et al. ()
/////////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2020-07-20
/////////////////////////
/////////////////////////

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

