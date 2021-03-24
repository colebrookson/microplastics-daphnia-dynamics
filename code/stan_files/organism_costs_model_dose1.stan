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
    //real Lm = theta[1];
    //real ke = theta[2];
    //real gamma = theta[3];
    //real R_m = theta[4];
    //real lp = theta[5];
    //real cstar = theta[6];  
    //real NEC = theta[7];
// define the functions first 
functions { // dz_dt holds all state variables (in our case 6)
  real[] dcq_dt(real t, 
               real[] z_cq, // specifying the output   
               real[] theta_cq, 
               real[] x_r,  
               int[] x_i) {
    real cq = z_cq[1];
    //real l = z_cq[2];
    //real R = z[3];
    
    real ke = theta_cq[1];
   // real gamma = theta[2];

    
    //real l = LL/Lm;
    real d_cq = ke*(0-cq);
    //real d_l = gamma*(1-l);
    //real R = (R_m/(1-(lp^3)))*(1*(l^2)*((1+l)/(1+1))-(lp^3))*((1+((cstar^-1)*fmax(0, (cq-NEC))))^-1);
    
    return {d_cq};
  }
  real[] dll_dt(real t, 
               real[] z_ll, // specifying the output   
               real[] theta_ll, 
               real[] x_r,  
               int[] x_i) {
    //real cq = z[1];
    real l = z_ll[1];
    //real R = z[3];
    
    //real ke = theta[1];
    real gamma = theta_ll[1];

    
    //real l = LL/Lm;
    //real d_cq = ke*(0-cq);
    real d_l = gamma*(1-l);
    //real R = (R_m/(1-(lp^3)))*(1*(l^2)*((1+l)/(1+1))-(lp^3))*((1+((cstar^-1)*fmax(0, (cq-NEC))))^-1);
    
    return {d_l};
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  // length data
  real ll_init[1]; // initial length value 
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  // the following two arrays contain the indexes of the final array l_y which
  // is the observed data (coded as a data vector l_y_obs) and the missing data
  // (coded as a parameter vector )
  int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs]; 
  int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];
  real l_y_obs[N_obs];
  
  // reproduction data
  real ts[21]; // time points
  real cq_init[1]; // initial value for cq equation (0)
  int<lower = 0> rep_r; // number of replicates in the reproduction trials
  real r_y[21]; // cumulative reproduction  
  
  ///////// BEGIN NOTE /////////////////////////////////////////////////////////
  // In a complex system you have to understand both the parameter values 
  // and abundances while still getting the dynamics wrong because you start in
  // a different place. So passing it the initial value is helping the model 
  // find the correct place to start, from which you'll have a much easier time
  // getting the model to converge. You can then test both ways
  ///////// END NOTE ///////////////////////////////////////////////////////////
}
transformed data {
  
  int<lower = 0> N = N_obs + N_mis;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower = 0> theta_cq[1];   // ke
  real<lower = 0> theta_ll[1]; // gamma
  real<lower = 0> z_init_cq[1];  // initial values
  real<lower = 0> z_init_ll[1];
  real l_y_mis[N_mis];
  real<lower = 0> Lp;
  real<lower = 0> Rm;
  real<lower = 0> cstar;
  real<lower = 0> nec;
  real<lower = 0> Lm;
  real tau_l;
  real tau_r;
  
  //real<lower = 0> sigma[2];   // error scale
}
transformed parameters {
  real z_cq[N]
    = integrate_ode_rk45(dcq_dt, z_init_cq, 0, ts, theta_cq,
                          rep_array(0.0, 0), rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll[N]
    = integrate_ode_rk45(dll_dt, z_init_ll, 0, ts, theta_ll,
                          rep_array(0.0, 0), rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real l_y[N]
  l_y[ii_obs] = l_y_obs;
  l_y[ii_mis] = l_y_mis;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  theta_cq[{1}] ~ uniform(0.5, 5); // ke
  theta_ll[{1}] ~ normal(0.11, 0.009); //gamma
  Lp ~ normal(0.49, 0.049); // length at puberty
  Rm ~ normal(10.74, 13.1044); // max reproduction
  cstar ~ uniform(0,10000);
  nec ~ uniform(0,6000);
  Lm ~ normal(4.77, 1.98);
  tau_l ~ gamma(0.001, 0.001);
  tau_r ~ gamma(0.001, 0.001)
  
 // theta[{2}] ~ normal(0.5,0.5); // cq
  //theta[{3}] ~ uniform(0, 500); // NEC
  //theta[{4}] ~ uniform(0.5, 5); // ke
  //sigma ~ lognormal(-1, 1);
  //z_init ~ lognormal(log(140), 1);
  
  real R0[N, rep_r];
  real eq0[N, rep_r];
  real L0[N, 10];
  
  // do the reproduction estimation 
  
  for(x in 1:rep_r){ // for every replicate
    
    R0[1,x] = 0; // initialization of cumulated reproduction for the control
    
    for(y in 2:N){ // every day from 2 to 21
    
      // equation that is either 0 or 1, 1 if the scaled length is > Lp
      eq0[y,x] = fmax(0,(z_ll[y])-Lp)/sqrt((z_ll[y]-Lp)^2);
      
      // cumulative reproduction at each time step
      R0[y,x] = R0[y-1,x] + eq0[y,x]*(Rm/(1-(Lp^3))*((1*(z_ll[y]^2)) *
                ((1+z_ll[y])/(1+1))-(Lp^3))*(1+((cstar^-1)*(fmax(0,(z_cq[y]-nec))))^-1);
      // fit R
      Rep0[y,x] ~ normal(R0[y,x], tau_r);
    }
    
  }
  
  // do the length estimation now
  
  for(i in 1:N){
    
    for(j in 1:10){
      
      // get the theoretical length
      L0[i,j] = Lm*z_ll[i];
      
      // fit observed length
      l_y[i,j] ~ normal(L0[i,j], tau_g);
      
    }
  }
}
// Uncertainty due to parameter estimation is rolled into the values of z_init
// z, and sigma. The uncertainty due to unexplained variation and measurement 
// error is captured through the use of the normal pseudorandom number generator
// *normal_rng*. The additional noise in the measurements y over that of the 
// underlying population predictions is visualized in plots
generated quantities {
  
  real R0_rep[N, rep_r];
  real eq0_rep[N, rep_r];
  real L0_rep[N, 10];
  
  // do the reproduction estimation 
  
  for(x in 1:rep_r){ // for every replicate
    
    R0_rep[1,x] = 0; // initialization of cumulated reproduction for the control
    
    for(y in 2:N){ // every day from 2 to 21
    
      // equation that is either 0 or 1, 1 if the scaled length is > Lp
      eq0_rep[y,x] = fmax(0,(z_ll[y])-Lp)/sqrt((z_ll[y]-Lp)^2);
      
      // cumulative reproduction at each time step
      R0_rep[y,x] = R0_rep[y-1,x] + eq0_rep[y,x]*(Rm/(1-(Lp^3))*((1*(z_ll[y]^2)) *
                ((1+z_ll[y])/(1+1))-(Lp^3))*(1+((cstar^-1)*(fmax(0,(z_cq[y]-nec))))^-1);
      // fit R
      Rep0[y,x] ~ normal(R0_rep[y,x], tau_r);
    }
    
  }
  
  // do the length estimation now
  
  for(i in 1:N){
    
    for(j in 1:10){
      
      // get the theoretical length
      L0_rep[i,j] = Lm*z_ll[i];
      
      // fit observed length
      l_y[i,j] ~ normal(L0_rep[i,j], tau_g);
      
    }
  }
}


