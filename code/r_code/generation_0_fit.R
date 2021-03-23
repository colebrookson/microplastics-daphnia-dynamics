########## 
##########
# This code contains some data preparation and model fitting for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-03-22
##########
##########

# set-up =======================================================================

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)
library(parallel)
library(R2WinBUGS)
detectCores()


dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data/' #private data repo

growth_data = read_csv(paste0(dir, 'growth_data_gen0.csv'))
reproduction_data = read_csv(paste0(dir, 'offspring_data.csv'))

# keep only the LFC and generation 0
reproduction_data = reproduction_data %>% 
  select(-`X1`) %>% 
  filter(treatment == 'LFC') %>% 
  filter(generation == 0)
reproduction_data$offspring = 
  as.numeric(reproduction_data$offspring)

# make the length into mm
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment == 'LFC') %>% #keep only control treatment
  mutate(length_mm = length*0.001)  

# keep only the first replicate for each and only 21 for reproduction
reproduction_data_meanreps = reproduction_data %>% 
  group_by(day) %>% 
  summarize(mean_offspring = mean(offspring, na.rm = TRUE))
  
growth_data = growth_data %>% 
  filter(day < 22) %>% 
  group_by(day) %>% 
  summarize(mean_length = mean(length_mm, na.rm = TRUE))

growth_repro_data = left_join(reproduction_data_meanreps, 
                              growth_data, 
                              by = 'day')
growth_repro_data[is.na(growth_repro_data['mean_offspring']),
                  'mean_offspring']=0

# add concentration (0)
growth_data$conc = 0
# prep data ====================================================================

N = length(growth_data$day)
ts = 1:N
y_init = c(growth_data$conc[1], 
           growth_data$mean_length[1])
y = as.matrix(cbind(growth_data$conc,
                    growth_data$mean_length))
y = cbind(y[ , 1], y[ , 2]); 
gen0_fit_example = list(N = N, 
                        ts = ts,
                        y_init = y_init,
                        y = y)

# fit model ====================================================================
gen_0_fit = stan(file = here('./code/stan_files/organism_modeling.stan'),
               data = gen0_fit_example,
               chains = 8,
               cores = 8,
               warmup = 5000,
               iter = 20000,
               seed = 12,
               verbose = TRUE,
               #open_progress = TRUE,
               control = list(adapt_delta = 0.9999)); beep(3)
### diagnose model fit
gen_0_fit_summ = print(gen_0_fit, 
                     pars=c("theta", "sigma", "z_init"),
                     probs=c(0.1, 0.5, 0.9), digits = 3)

parms <- c("theta[1]","theta[2]","theta[3]","theta[4]")
gen_0_fit_output1 <- rstan::extract(gen_0_fit,permuted=TRUE,include=TRUE)

# Parms
nec_fit_Trace <- stan_trace(gen_0_fit,parms)
nec_fit_Dens <- mcmc_dens(gen_0_fit,parms)
nec_fit_Overlay <- mcmc_dens_overlay(gen_0_fit,parms)
nec_fit_Violin <- mcmc_violin(gen_0_fit,parms,probs = c(0.1, 0.5, 0.9))
nec_fit_Pairs <- mcmc_pairs(gen_0_fit,parms)
nec_fit_Trace
nec_fit_Dens
nec_fit_Overlay
nec_fit_Violin
nec_fit_Pairs

# try model in WinBUGS =========================================================
sink("gen_0_control_bugs.txt")
cat("
#...

model {
# system of differential equations specified via BUGS language : CQ corresponds to the
# scaled internal concentration, and LL to the scaled length. CQ0 and LL0 correspond to the
# control. CQ1 and LL1 correspond to the first exposure concentration, etc… cext is a vector
#containing the successive exposure concentrations.

solution[1:n.grid, 1:dim] <- ode(init[1:dim], grid[1:n.grid], D(C[1:dim], t), origin, tol)

D(C[CQ0], t) <- (cext[1] - C[CQ0] ) * ke
D(C[LL0], t) <- gam * (f - C[LL0] * (1+ pow(cm,-1) * max(0,(C[CQ0] - nec) ) ) )

# Initial conditions for the differential equations: Lb corresponds to the scaled length at birth

init[CQ0] <- 0; init[LL0] <- Lb

# Reproduction model: R0 corresponds to the theoretical cumulated number of offspring
# for controls, R1 for the first exposure concentration, etc..., N corresponds to the number of
# replicates, and h corresponds to the different times

  for(N in 1:1) {
  
  R0[1,N] <- 0
  
    for(h in 2:21){ # every day from 2 to 21
# eq = 0 if the scaled length LL is < lp, the scaled length at puberty, and eq = 1 if the
# scaled length LL is >= lp

    eq0[h,N] <- max(0,(solution[h,LL0]-Lp)/sqrt(pow(solution[h,LL0]-Lp,2)))

  R0[h,N] <- R0[h-1,N]+eq0[h,N]*(1+ pow(cm,-1) * max(0,(solution[h,CQ0] - nec) ))*Rm/(1-
pow(Lp,3))*((g+solution[h,LL0])/(g+f)*f*pow(solution[h,LL0],2)-pow(Lp,3))

  R0[h,N] <- R0[h-1,N]+eq0[h,N]*(1+ pow(cm,-1) * max(0,(solution[h,CQ0] - nec) ))*Rm/(1-
pow(Lp,3))*((g*pow(1+ pow(cm,-1) * max(0,(solution[h,CQ0] - nec) ),-
1)+solution[h,LL0])/(g+f)*f*pow(solution[h,LL0],2)-pow(Lp,3
    }
  }













    ", fill = TRUE)
sink()
