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
functions { // dz_dt holds all state variables (in our case 6)
  real[] dll_dt(real t, 
               real[] z_ll, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l = z_ll[1];

    real gamma = theta_ll[1];
    
    real dl_dt = gamma*(1-l);

    return { dl_dt };
  }
}
data {
  // length data
  real ll_init[1]; // initial length value 
  
  // reproduction data
  real ts[22]; // time points
  real r_y[22]; // cumulative reproduction  
  real l_y_obs[22];

}

parameters {
  real<lower = 0> theta_ll[1]; // gamma & l
  real<lower = 0> Lp;
  real<lower = 0> Rm;
  real<lower = 0> Lm;
  real<lower = 0> tau_l;
  real<lower = 0> tau_r;
  
  //real<lower = 0> sigma[2];   // error scale
}
transformed parameters {
  
  real z_ll[22,1] = 
      integrate_ode_rk45(dll_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
}
model {
  // definitions of values
  real R0[22];
  real eq0[22];
  real L0[22];

  // priors
  theta_ll[1] ~ normal(0.11, 0.009); //gamma
  Lp ~ normal(0.49, 0.14); // length at puberty
  Rm ~ normal(10.74, 13.1); // max reproduction
  Lm ~ normal(4.77, 1.2);
  tau_l ~ gamma(0.001, 0.01);
  tau_r ~ gamma(0.001, 0.01);
  // priors will be added for concentration 
  
 // theta[{2}] ~ normal(0.5,0.5); // cq
  //theta[{3}] ~ uniform(0, 500); // NEC
  //theta[{4}] ~ uniform(0.5, 5); // ke
  //sigma ~ lognormal(-1, 1);
  //z_init ~ lognormal(log(140), 1)

  
   // do the length estimation now
  for(i in 1:22){ // days
      
      real z_ll_temp = z_ll[i,1];
      // get the theoretical length
      L0[i] = Lm*z_ll_temp;
      
      // fit observed length
      l_y_obs[i] ~ normal(L0[i], tau_l);
      
    }
  // do the reproduction estimation 
  R0[1] = 0; // initialization of cumulated reproduction for the control
  eq0[1] = 0; // helper value/equation
    
  for(y in 2:22){ // every day from 2 to 21
      
    real z_temp = z_ll[y,1];// value from the ode solver 
      // equation that is either 0 or 1, 1 if the scaled length is > Lp
    if (z_temp <= Lp)
      eq0[y] = 0;
    else
      eq0[y] = 1;
      
      // cumulative reproduction at each time step
    R0[y] = R0[y-1] + eq0[y]*(Rm/(1-(Lp^3)))*((1*((z_temp)^2)) *
                ((1+z_temp)/(1+1))-(Lp^3));
      // fit R
    r_y[y] ~ normal(R0[y], tau_r); // estimation step 
  }
}
// Uncertainty due to parameter estimation is rolled into the values of z_init
// z, and sigma. The uncertainty due to unexplained variation and measurement 
// error is captured through the use of the normal pseudorandom number generator
// *normal_rng*. The additional noise in the measurements y over that of the 
// underlying population predictions is visualized in plots
generated quantities {
  
  //real z_ll_rep[22,1];
  real R0_rep[22];
  real eq0_rep[22];
  real L0_rep[22];
  real r_y_rep[22];
  real l_y_rep[22];
  
   //z_ll_rep
    // = integrate_ode_rk45(dll_dt, ll_init, 0, ts, theta_ll, x_r, x_i);
  
  // do the reproduction estimation 
    
  R0_rep[1] = 0; // initialization of cumulated reproduction for the control
  eq0_rep[1] = 0;
  r_y_rep[1]= 0;
    
  for(y in 2:22){ // every day from 2 to 21
      
    real z_temp_rep = z_ll[y,1];
      // equation that is either 0 or 1, 1 if the scaled length is > Lp
    if (z_temp_rep <= Lp)
      eq0_rep[y] = 0;
    else
      eq0_rep[y] = 1;
      // cumulative reproduction at each time step
    R0_rep[y] = R0_rep[y-1] + eq0_rep[y]*(Rm/(1-(Lp^3)))*((1*((z_temp_rep)^2)) *
                ((1+z_temp_rep)/(1+1))-(Lp^3));
      // fit R
    r_y_rep[y] = normal_rng(R0_rep[y], tau_r);
  }
  
  // do the length estimation now

    
  for(i in 1:22){
      
    real z_ll_temp_rep = z_ll[i,1];
      // get the theoretical length
    L0_rep[i] = Lm*z_ll_temp_rep;
      
      // fit observed length
    l_y_rep[i] = normal_rng(L0_rep[i], tau_l);
      
    }
}


