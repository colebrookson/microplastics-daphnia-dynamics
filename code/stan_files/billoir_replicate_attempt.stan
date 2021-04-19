//////////////////// 
////////////////////
// This code contains the stan code that attempts to replicate the findings of 
// Billoir et al (2008) who performed their analysis in WinBUGS
////////////////////
////////////////////
// AUTHOR: Cole B. Brookson
// DATE OF CREATION: 2021-04-17
////////////////////
////////////////////

functions { // dz_dt holds all state variables (in our case 6)
  real[] dll_con_dt(real t, 
               real[] z_ll_con, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_con = z_ll_con[1];
    real cq_con = z_ll_con[2];
    

    real gamma = theta_ll[1];
    real cstar = theta_ll[2];
    real NEC = theta_ll[3];
    real ke = theta_ll[4];
    
    real dl_con_dt = gamma*(1-l_con*
                    (1+(cstar^(-1))*(fmax(0, (cq_con-NEC)))));
    real d_con_cq = ke*(0-cq_con);

    return { dl_con_dt, d_con_cq };
  }
  real[] dll_400_dt(real t, 
               real[] z_ll_400, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_400 = z_ll_400[1];
    real cq_400 = z_ll_400[2];
    

    real gamma = theta_ll[1];
    real cstar = theta_ll[2];
    real NEC = theta_ll[3];
    real ke = theta_ll[4];
    
    real dl_400_dt = gamma*(1-l_400*
                    (1+(cstar^(-1))*(fmax(0, (cq_400-NEC)))));
    real d_400_cq = ke*(0.074-cq_400);

    return { dl_400_dt, d_400_cq };
  }  real[] dll_2000_dt(real t, 
               real[] z_ll_2000, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_2000 = z_ll_2000[1];
    real cq_2000 = z_ll_2000[2];
    

    real gamma = theta_ll[1];
    real cstar = theta_ll[2];
    real NEC = theta_ll[3];
    real ke = theta_ll[4];
    
    real dl_2000_dt = gamma*(1-l_2000*
                    (1+(cstar^(-1))*(fmax(0, (cq_2000-NEC)))));
    real d_2000_cq = ke*(0.22-cq_2000);

    return { dl_2000_dt, d_2000_cq };
  }  real[] dll_10000_dt(real t, 
               real[] z_ll_10000, // specifying the output   
               real[] theta_ll, // parameters for this function 
               real[] x_r,  
               int[] x_i) {
    real l_10000 = z_ll_10000[1];
    real cq_10000 = z_ll_10000[2];
    

    real gamma = theta_ll[1];
    real cstar = theta_ll[2];
    real NEC = theta_ll[3];
    real ke = theta_ll[4];
    
    real dl_10000_dt = gamma*(1-l_10000*
                    (1+(cstar^(-1))*(fmax(0, (cq_10000-NEC)))));
    real d_10000_cq = ke*(0.66-cq_10000);

    return { dl_10000_dt, d_10000_cq };
  }
}
data {

  real ll_init[2]; // initial length value 
  //real cq_init[1]; // initial concentration value
  // real cq_init[1];

  // reproduction data
  real ts[21]; // time points
  real r_y_con[21]; // cumulative reproduction  
  real r_y_400[21]; // cumulative reproduction  
  real r_y_2000[21]; // cumulative reproduction  
  real r_y_10000[21]; // cumulative reproduction  
  
  real l_y_obs_con[21];
  real l_y_obs_400[21];
  real l_y_obs_2000[21];
  real l_y_obs_10000[21];
  //real cq[22];

}
transformed data {
  
  //int<lower = 0> N = N_obs + N_mis;
  //real x_r[1];
  //int x_i[0];
  
  
}
parameters {
  real<lower = 0> theta_ll[4]; // gamma 
  //real<lower = 0, upper = 20000> theta_ll[{2}]; // cstar
  //real<lower = 0> theta_ll[{3}]; // NEC
  //real<lower = 0> theta_ll[{4}]; // ke
  real<lower = 0.23> Lp;
  real<lower = 0> Rm;
  real<lower = 1> Lm;
  real<lower = 0> tau_l;
  real<lower = 0> tau_r;
  
  //real<lower = 0> sigma[2];   // error scale
}
transformed parameters {
  
  real z_ll_con[21,2] = 
      integrate_ode_rk45(dll_con_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_400[21,2] = 
      integrate_ode_rk45(dll_400_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_2000[21,2] = 
      integrate_ode_rk45(dll_2000_dt, // function (defined above)
                         ll_init, // initial length value
                         0, // initial time point 
                         ts, // time series to integrate over
                         theta_ll, // parameter vector
                         rep_array(0.0, 0), 
                         rep_array(0, 0),
                          1e-5, 1e-3, 5e2);
  real z_ll_10000[21,2] = 
      integrate_ode_rk45(dll_10000_dt, // function (defined above)
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
  real R0[21];
  real eq0[21];
  real L0[21];
  real R1[21];
  real eq1[21];
  real L1[21];
  real R2[21];
  real eq2[21];
  real L2[21];
  real R3[21];
  real eq3[21];
  real L3[21];

  // priors
  theta_ll[1] ~ normal(0.11, 0.009); //gamma
  theta_ll[2] ~ uniform(0,3); //cstar
  theta_ll[3] ~ uniform(0,0.66); // NEC
  theta_ll[4] ~ uniform(0.5, 5); //ke
  //cstar ~ uniform(0,3); // tolerance concentration
  //NEC ~ uniform(0,0.66); // no effect concentration 
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
  
  for(y in 2:21){ // every day from 2 to 21
      
    real z_ll_con_temp = z_ll_con[y,1];// value from the ode solver  for length
    //real z_cq_con_temp = z_con_cq[y,1]; // value from ode solver for cq
    real z_ll_400_temp = z_ll_400[y,1];
    real z_cq_400_temp = z_ll_400[y,2]; 
    real z_ll_2000_temp = z_ll_2000[y,1];
    real z_cq_2000_temp = z_ll_2000[y,2]; 
    real z_ll_10000_temp = z_ll_10000[y,1];
    real z_cq_10000_temp = z_ll_10000[y,2]; 
    
    real s_cq_con = 0; 
    real s_cq_400 = (theta_ll[2]^(-1))*(fmax(0, (z_cq_400_temp-theta_ll[3])));
    real s_cq_2000 = (theta_ll[2]^(-1))*(fmax(0, (z_cq_2000_temp-theta_ll[3])));
    real s_cq_10000 = (theta_ll[2]^(-1))*(fmax(0, (z_cq_10000_temp-theta_ll[3])));
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
    R0[y] = R0[y-1] + eq0[y] *
            (1+s_cq_con) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_con_temp^2))*
              ((1*((1+s_cq_con)^-1)+z_ll_con_temp)/
              (1+1)) - Lp^3
            );
    R1[y] = R1[y-1] + eq1[y] *
            (1+s_cq_400) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_400_temp^2))*
              ((1*((1+s_cq_400)^-1)+z_ll_400_temp)/
              (1+1)) - Lp^3
            );   
    R2[y] = R2[y-1] + eq2[y] *
            (1+s_cq_2000) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_2000_temp^2))*
              ((1*((1+s_cq_2000)^-1)+z_ll_2000_temp)/
              (1+1)) - Lp^3
            );
    R3[y] = R3[y-1] + eq3[y] *
            (1+s_cq_10000) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_10000_temp^2))*
              ((1*((1+s_cq_10000)^-1)+z_ll_10000_temp)/
              (1+1)) - Lp^3
            );
            
    // fit R
    r_y_con[y] ~ normal(R0[y], tau_r); // estimation step 
    r_y_400[y] ~ normal(R1[y], tau_r);
    r_y_2000[y] ~ normal(R2[y], tau_r);
    r_y_10000[y] ~ normal(R3[y], tau_r);
  }

  
  // do the length estimation now
  
  for(i in 1:21){ // days
      
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

generated quantities {
  
  real R0_rep[21];
  real eq0_rep[21];
  real L0_rep[21];
  real r_y_con_rep[21];
  real l_y_con_rep[21];
  real R1_rep[21];
  real eq1_rep[21];
  real L1_rep[21];
  real r_y_400_rep[21];
  real l_y_400_rep[21]; 
  real R2_rep[21];
  real eq2_rep[21];
  real L2_rep[22];
  real r_y_2000_rep[21];
  real l_y_2000_rep[21];
  real R3_rep[21];
  real eq3_rep[21];
  real L3_rep[21];
  real r_y_10000_rep[21];
  real l_y_10000_rep[21];
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
    
  for(y in 2:21){ // every day from 2 to 21
      
    real z_ll_con_temp_rep = z_ll_con[y,1];
    //real z_cq_con_temp_rep = z_con_cq[y,1]; // value from ode solver for cq
    real s_cq_con_rep = 0; 
    real z_ll_400_temp_rep = z_ll_400[y,1];
    real z_cq_400_temp_rep = z_ll_400[y,2]; 
    real s_cq_400_rep = (theta_ll[2]^(-1))*(fmax(0, (z_cq_400_temp_rep-theta_ll[3]))); 
    real z_ll_2000_temp_rep = z_ll_2000[y,1];
    real z_cq_2000_temp_rep = z_ll_2000[y,2]; 
    real s_cq_2000_rep = (theta_ll[2]^(-1))*(fmax(0, (z_cq_2000_temp_rep-theta_ll[3]))); 
    real z_ll_10000_temp_rep = z_ll_10000[y,1];
    real z_cq_10000_temp_rep = z_ll_10000[y,2]; 
    real s_cq_10000_rep = (theta_ll[2]^(-1))*(fmax(0, (z_cq_10000_temp_rep-theta_ll[3]))); 
  
    
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
    R0_rep[y] = R0_rep[y-1] + eq0_rep[y] *
            (1+s_cq_con_rep) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_con_temp_rep^2))*
              ((1*((1+s_cq_con_rep)^-1)+z_ll_con_temp_rep)/
              (1+1)) - Lp^3
            );
    R1_rep[y] = R1_rep[y-1] + eq1_rep[y] *
            (1+s_cq_400_rep) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_400_temp_rep^2))*
              ((1*((1+s_cq_400_rep)^-1)+z_ll_400_temp_rep)/
              (1+1)) - Lp^3
            );   
    R2_rep[y] = R2_rep[y-1] + eq2_rep[y] *
            (1+s_cq_2000_rep) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_2000_temp_rep^2))*
              ((1*((1+s_cq_2000_rep)^-1)+z_ll_2000_temp_rep)/
              (1+1)) - Lp^3
            );
    R3_rep[y] = R3_rep[y-1] + eq3_rep[y] *
            (1+s_cq_10000_rep) *
            (Rm/(1-(Lp^3)))*
            (
              (1*(z_ll_10000_temp_rep^2))*
              ((1*((1+s_cq_10000_rep)^-1)+z_ll_10000_temp_rep)/
              (1+1)) - Lp^3
            ); 
    // fit R
    r_y_con_rep[y] = normal_rng(R0_rep[y], tau_r);
    r_y_400_rep[y] = normal_rng(R1_rep[y], tau_r);
    r_y_2000_rep[y] = normal_rng(R2_rep[y], tau_r);
    r_y_10000_rep[y] = normal_rng(R3_rep[y], tau_r);
  }
  
  // do the length estimation now

    
  for(i in 1:21){
      
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
