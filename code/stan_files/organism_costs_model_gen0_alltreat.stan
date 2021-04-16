//////////////////// 
////////////////////
// This code contains the stan code for the individual component as a
// part of the larger analysis of the effect of microplastics on daphnia
////////////////////
////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2021-04-06
////////////////////
////////////////////
functions { // dz_dt holds all state variables (in our case 6)
  real[] dll_con_dt(real t, 
               real[] z_ll_con, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_con = z_ll_con[1];

    real gamma = theta_ll[1];
    
    real dl_con_dt = gamma*(1-l_con);

    return { dl_con_dt };
  }
  real[] dll_400_dt(real t, 
               real[] z_ll_400, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_400 = z_ll_400[1];

    real gamma = theta_ll[1];
    
    real dl_400_dt = gamma*(1-l_400);

    return { dl_400_dt };
  }
  real[] dll_2000_dt(real t, 
               real[] z_ll_2000, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_2000 = z_ll_2000[1];

    real gamma = theta_ll[1];
    
    real dl_2000_dt = gamma*(1-l_2000);

    return { dl_2000_dt };
  }
  real[] dll_10000_dt(real t, 
               real[] z_ll_10000, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_10000 = z_ll_10000[1];

    real gamma = theta_ll[1];
    
    real dl_10000_dt = gamma*(1-l_10000);

    return { dl_10000_dt };
  }
  real[] dcq_con_dt(real t, 
              real[] z_con_cq, // specifying the output   
              real[] theta_cq, 
              real[] x_r,  
              int[] x_i) {
   real cq_con = z_con_cq[1];
   
   real ke = theta_cq[1];
   
   real d_con_cq = ke*(400-cq_con); // NOTE - constant value changes w/ concentration
   
   return { d_con_cq };
 }
  real[] dcq_400_dt(real t, 
              real[] z_400_cq, // specifying the output   
              real[] theta_cq, 
              real[] x_r,  
              int[] x_i) {
   real cq_400 = z_400_cq[1];
   
   real ke = theta_cq[1];
   
   real d_400_cq = ke*(400-cq_400); // NOTE - constant value changes w/ concentration
   
   return { d_400_cq };
 }
  real[] dcq_2000_dt(real t, 
              real[] z_2000_cq, // specifying the output   
              real[] theta_cq, 
              real[] x_r,  
              int[] x_i) {
   real cq_2000 = z_2000_cq[1];
   
   real ke = theta_cq[1];
   
   real d_2000_cq = ke*(2000-cq_2000); // NOTE - constant value changes w/ concentration
   
   return { d_2000_cq };
 }
  real[] dcq_10000_dt(real t, 
              real[] z_10000_cq, // specifying the output   
              real[] theta_cq, 
              real[] x_r,  
              int[] x_i) {
   real cq_10000 = z_10000_cq[1];
   
   real ke = theta_cq[1];
   
   real d_10000_cq = ke*(10000-cq_10000); // NOTE - constant value changes w/ concentration
   
   return { d_10000_cq };
 }
}
data {

  real ll_init[1]; // initial length value 
  real cq_init[1]; // initial concentration value
  // real cq_init[1];

  // reproduction data
  real ts[22]; // time points
  real r_y_con[22]; // cumulative reproduction  
  real r_y_400[22]; // cumulative reproduction  
  real r_y_2000[22]; // cumulative reproduction  
  real r_y_10000[22]; // cumulative reproduction  
  
  real l_y_obs_con[22];
  real l_y_obs_400[22];
  real l_y_obs_2000[22];
  real l_y_obs_10000[22];
  //real cq[22];

}
transformed data {
  
  //int<lower = 0> N = N_obs + N_mis;
  //real x_r[1];
  //int x_i[0];
  
  
}
parameters {
  real<lower = 0> theta_ll[1]; // gamma & l
  real<lower = 0> theta_cq[1]; // ke
  real<lower = 0, upper = 10000> cstar;
  real<lower = 0> NEC;
  real<lower = 0> Lp;
  real<lower = 0> Rm;
  real<lower = 0> Lm;
  real<lower = 0> tau_l;
  real<lower = 0> tau_r;
  
  //real<lower = 0> sigma[2];   // error scale
}
transformed parameters {
  
  real z_ll_con[22,1] = 
      integrate_ode_rk45(dll_con_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_400[22,1] = 
      integrate_ode_rk45(dll_400_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_2000[22,1] = 
      integrate_ode_rk45(dll_2000_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_10000[22,1] = 
      integrate_ode_rk45(dll_10000_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_400_cq[22,1] = 
      integrate_ode_rk45(dcq_400_dt, // function (defined above)
                         cq_init, // initial cq value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_cq, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                        1e-5, 1e-3, 5e2);
  real z_2000_cq[22,1] = 
      integrate_ode_rk45(dcq_2000_dt, // function (defined above)
                         cq_init, // initial cq value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_cq, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                        1e-5, 1e-3, 5e2);
  real z_10000_cq[22,1] = 
      integrate_ode_rk45(dcq_10000_dt, // function (defined above)
                         cq_init, // initial cq value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_cq, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                        1e-5, 1e-3, 5e2);
    
}
model {
  // definitions of values
  real R0[22];
  real eq0[22];
  real L0[22];
  real R1[22];
  real eq1[22];
  real L1[22];
  real R2[22];
  real eq2[22];
  real L2[22];
  real R3[22];
  real eq3[22];
  real L3[22];

  // priors
  theta_ll[1] ~ normal(0.11, 0.009); //gamma
  theta_cq[1] ~ normal(0, 2); //ke
  log(cstar) ~ normal(0, 5000); // tolerance concentration
  log(NEC) ~ normal(0, 2000); // no effect concentration 
  Lp ~ normal(0.49, 0.049); // length at puberty
  Rm ~ normal(10.74, 13.1044); // max reproduction
  Lm ~ normal(4.77, 1.98);
  tau_l ~ gamma(0.001, 0.001);
  tau_r ~ gamma(0.001, 0.001);

  // do the reproduction estimation 
    
  R0[1] = 0; // initialization of cumulated reproduction for the control
  eq0[1] = 0; // helper value/equation
  R1[1] = 0; 
  eq1[1] = 0;
  R2[1] = 0; 
  eq2[1] = 0;
  R3[1] = 0; 
  eq3[1] = 0;
  
  for(y in 2:22){ // every day from 2 to 21
      
    real z_ll_con_temp = z_ll_con[y,1];// value from the ode solver  for length
    //real z_cq_con_temp = z_con_cq[y,1]; // value from ode solver for cq
    real z_ll_400_temp = z_ll_400[y,1];
    real z_cq_400_temp = z_400_cq[y,1]; 
    real z_ll_2000_temp = z_ll_2000[y,1];
    real z_cq_2000_temp = z_2000_cq[y,1]; 
    real z_ll_10000_temp = z_ll_10000[y,1];
    real z_cq_10000_temp = z_10000_cq[y,1]; 
    
    real s_cq_con = 0; 
    real s_cq_400 = (cstar^(-1))*(fmax(0, (z_cq_400_temp-NEC)));
    real s_cq_2000 = (cstar^(-1))*(fmax(0, (z_cq_2000_temp-NEC)));
    real s_cq_10000 = (cstar^(-1))*(fmax(0, (z_cq_10000_temp-NEC)));
    //print("s_cq_400: ", s_cq_400, "s_cq_2000: ", s_cq_2000, 
    //"s_cq_10000: ", s_cq_10000)
    // equation that is either 0 or 1, 1 if the scaled length is > Lp
    if (z_ll_con_temp <= Lp)
      eq0[y] = 0;
    else
      eq0[y] = 1;
      
    if (z_ll_400_temp <= Lp)
      eq1[y] = 0;
    else
      eq1[y] = 1;
      
    if (z_ll_2000_temp <= Lp)
      eq2[y] = 0;
    else
      eq2[y] = 1;
      
    if (z_ll_10000_temp <= Lp)
      eq3[y] = 0;
    else
      eq3[y] = 1;
     
    // cumulative reproduction at each time step
    R0[y] = R0[y-1] + eq0[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_con_temp)^2)) *
                ((1+z_ll_con_temp)/(1+1))-(Lp^3))*((1+s_cq_con)^-1);
                
    R1[y] = R1[y-1] + eq1[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_400_temp)^2)) *
                ((1+z_ll_400_temp)/(1+1))-(Lp^3))*((1+s_cq_400)^-1);
                
    R2[y] = R2[y-1] + eq2[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_2000_temp)^2)) *
                ((1+z_ll_2000_temp)/(1+1))-(Lp^3))*((1+s_cq_2000)^-1);
                
    R3[y] = R3[y-1] + eq3[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_10000_temp)^2)) *
                ((1+z_ll_10000_temp)/(1+1))-(Lp^3))*((1+s_cq_10000)^-1);
                
    // fit R
    r_y_con[y] ~ normal(R0[y], tau_r); // estimation step 
    r_y_400[y] ~ normal(R1[y], tau_r);
    r_y_2000[y] ~ normal(R2[y], tau_r);
    r_y_10000[y] ~ normal(R3[y], tau_r);
  }

  
  // do the length estimation now
  
  for(i in 1:22){ // days
      
      real z_ll_con_temp = z_ll_con[i,1];
      real z_ll_400_temp = z_ll_400[i,1];
      real z_ll_2000_temp = z_ll_2000[i,1];
      real z_ll_10000_temp = z_ll_10000[i,1];
      // get the theoretical length
      L0[i] = Lm*z_ll_con_temp;
      L1[i] = Lm*z_ll_400_temp;
      L2[i] = Lm*z_ll_2000_temp;
      L3[i] = Lm*z_ll_10000_temp;
      
      // fit observed length
      l_y_obs_con[i] ~ normal(L0[i], tau_l);
      l_y_obs_400[i] ~ normal(L1[i], tau_l);
      l_y_obs_2000[i] ~ normal(L2[i], tau_l);
      l_y_obs_10000[i] ~ normal(L3[i], tau_l);
    }
}
// Uncertainty due to parameter estimation is rolled into the values of z_init
// z, and sigma. The uncertainty due to unexplained variation and measurement 
// error is captured through the use of the normal pseudorandom number generator
// *normal_rng*. The additional noise in the measurements y over that of the 
// underlying population predictions is visualized in plots
generated quantities {
  
  real R0_rep[22];
  real eq0_rep[22];
  real L0_rep[22];
  real r_y_con_rep[22];
  real l_y_con_rep[22];
  real R1_rep[22];
  real eq1_rep[22];
  real L1_rep[22];
  real r_y_400_rep[22];
  real l_y_400_rep[22]; 
  real R2_rep[22];
  real eq2_rep[22];
  real L2_rep[22];
  real r_y_2000_rep[22];
  real l_y_2000_rep[22];
  real R3_rep[22];
  real eq3_rep[22];
  real L3_rep[22];
  real r_y_10000_rep[22];
  real l_y_10000_rep[22];
  // do the reproduction estimation 
    
  R0_rep[1] = 0; // initialization of cumulated reproduction for the control
  eq0_rep[1] = 0;
  r_y_con_rep[1]= 0;
  R1_rep[1] = 0; // initialization of cumulated reproduction for the control
  eq1_rep[1] = 0;
  r_y_400_rep[1]= 0;
  R2_rep[1] = 0; // initialization of cumulated reproduction for the control
  eq2_rep[1] = 0;
  r_y_2000_rep[1]= 0;
  R3_rep[1] = 0; // initialization of cumulated reproduction for the control
  eq3_rep[1] = 0;
  r_y_10000_rep[1]= 0;
    
  for(y in 2:22){ // every day from 2 to 21
      
    real z_ll_con_temp_rep = z_ll_con[y,1];
    //real z_cq_con_temp_rep = z_con_cq[y,1]; // value from ode solver for cq
    real s_cq_con_rep = 0; 
    real z_ll_400_temp_rep = z_ll_400[y,1];
    real z_cq_400_temp_rep = z_400_cq[y,1]; 
    real s_cq_400_rep = (cstar^(-1))*(fmax(0, (z_cq_400_temp_rep-NEC))); 
    real z_ll_2000_temp_rep = z_ll_2000[y,1];
    real z_cq_2000_temp_rep = z_2000_cq[y,1]; 
    real s_cq_2000_rep = (cstar^(-1))*(fmax(0, (z_cq_2000_temp_rep-NEC))); 
    real z_ll_10000_temp_rep = z_ll_10000[y,1];
    real z_cq_10000_temp_rep = z_10000_cq[y,1]; 
    real s_cq_10000_rep = (cstar^(-1))*(fmax(0, (z_cq_10000_temp_rep-NEC))); 
  
    
    // equation that is either 0 or 1, 1 if the scaled length is > Lp
    if (z_ll_con_temp_rep <= Lp)
      eq0_rep[y] = 0;
    else
      eq0_rep[y] = 1; 
      
    if (z_ll_400_temp_rep <= Lp)
      eq1_rep[y] = 0;
    else
      eq1_rep[y] = 1; 
      
    if (z_ll_2000_temp_rep <= Lp)
      eq2_rep[y] = 0;
    else
      eq2_rep[y] = 1; 
      
    if (z_ll_10000_temp_rep <= Lp)
      eq3_rep[y] = 0;
    else
      eq3_rep[y] = 1; 
            
    // cumulative reproduction at each time step
    R0_rep[y] = R0_rep[y-1] + eq0_rep[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_con_temp_rep)^2)) *
                ((1+z_ll_con_temp_rep)/(1+1))-(Lp^3))*((1+s_cq_con_rep)^-1);
    R1_rep[y] = R1_rep[y-1] + eq1_rep[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_400_temp_rep)^2)) *
                ((1+z_ll_400_temp_rep)/(1+1))-(Lp^3))*((1+s_cq_400_rep)^-1);
    R2_rep[y] = R2_rep[y-1] + eq2_rep[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_2000_temp_rep)^2)) *
                ((1+z_ll_2000_temp_rep)/(1+1))-(Lp^3))*((1+s_cq_2000_rep)^-1);
    R3_rep[y] = R3_rep[y-1] + eq3_rep[y]*(Rm/(1-(Lp^3)))*((1*((z_ll_10000_temp_rep)^2)) *
                ((1+z_ll_10000_temp_rep)/(1+1))-(Lp^3))*((1+s_cq_10000_rep)^-1);
                
    // fit R
    r_y_con_rep[y] = normal_rng(R0_rep[y], tau_r);
    r_y_400_rep[y] = normal_rng(R1_rep[y], tau_r);
    r_y_2000_rep[y] = normal_rng(R2_rep[y], tau_r);
    r_y_10000_rep[y] = normal_rng(R3_rep[y], tau_r);
  }
  
  // do the length estimation now

    
  for(i in 1:22){
      
    real z_ll_con_temp_rep = z_ll_con[i,1];
    real z_ll_400_temp_rep = z_ll_400[i,1];
    real z_ll_2000_temp_rep = z_ll_2000[i,1];
    real z_ll_10000_temp_rep = z_ll_10000[i,1];
    
    // get the theoretical length
    L0_rep[i] = Lm*z_ll_con_temp_rep;
    L1_rep[i] = Lm*z_ll_400_temp_rep;
    L2_rep[i] = Lm*z_ll_2000_temp_rep;
    L3_rep[i] = Lm*z_ll_10000_temp_rep;
      
    // fit observed length
    l_y_con_rep[i] = normal_rng(L0_rep[i], tau_l);
    l_y_400_rep[i] = normal_rng(L1_rep[i], tau_l);
    l_y_2000_rep[i] = normal_rng(L2_rep[i], tau_l);
    l_y_10000_rep[i] = normal_rng(L3_rep[i], tau_l);
      
    }
}

