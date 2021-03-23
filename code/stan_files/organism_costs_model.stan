//////////////////// 
////////////////////
// This code contains the stan code for the individual component as a
// part of the larger analysis of the effect of microplastics on daphnia
////////////////////
////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2021-03-18
////////////////////
////////////////////

// define the functions first 
functions { // dz_dt holds all state variables (in our case 6)
  real[] dz_dt(real t, 
               real[] z, // specifying the output   
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real l = z[1];
    real c = z[2];

    real cstar = theta[1];  
    real cq = theta[2];
    real NEC = theta[3];
    real ke = theta[4];

    real d_cq = ke*(c-cq);
    real d_LL = 1*((1+1)/(1+1*(1+(cstar^-1)*fmax(0, (cq-NEC)))))*(1-l);

    // real d_R = (R_m/(1-(lp^3)))*(1*(l^2)(1*(1+(cstar*(max(0, (cq-nec))))))/())
    // real d_PS400_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    // real d_PS400_cq = ke*(c-cq);
    // real d_PS2000_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    // real d_PS2000_cq = ke*(c-cq);
    // real d_PS10000_ll = 1*((1+1)/(1+1*(1+(cstar*(cq-NEC)))))*(1-l);
    // real d_PS10000_cq = ke*(c-cq);
    return {d_cq, d_LL};
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // length of dataset - 1
  real ts[N]; // time points
  real y_init[2];              
  real<lower = 0> y[N, 2];
  
  ///////// BEGIN NOTE /////////////////////////////////////////////////////////
  // In a complex system you have to understand both the parameter values 
  // and abundances while still getting the dynamics wrong because you start in
  // a different place. So passing it the initial value is helping the model 
  // find the correct place to start, from which you'll have a much easier time
  // getting the model to converge. You can then test both ways
  ///////// END NOTE ///////////////////////////////////////////////////////////
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower = 0> theta[4];   // theta = {cstar, cq, nec, ke}
  real<lower = 0> z_init[2];  // initial values
  real<lower = 0> sigma[2];   // error scale
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                          rep_array(0.0, 0), rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  theta[{1}] ~ uniform(0, 8000); // cstar
  theta[{2}] ~ normal(0.5,0.5); // cq
  theta[{3}] ~ uniform(0, 500); // NEC
  theta[{4}] ~ uniform(0.5, 5); // ke
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(140), 1);
  for (k in 1:2) {
    y_init[k] ~ normal(z_init[k], sigma[k]);
    y[ , k] ~ normal(z[, k], sigma[k]);
  }
}
// Uncertainty due to parameter estimation is rolled into the values of z_init
// z, and sigma. The uncertainty due to unexplained variation and measurement 
// error is captured through the use of the normal pseudorandom number generator
// *normal_rng*. The additional noise in the measurements y over that of the 
// underlying population predictions is visualized in plots
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = normal_rng(z_init[k], sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = normal_rng(z[n, k], sigma[k]);
  }
}


