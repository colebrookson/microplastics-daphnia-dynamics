########## 
##########
# This code contains the stan code that attempts to replicate the findings of 
# Billoir et al (2008) who performed their analysis in WinBUGS
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-04-17
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
library(deSolve)
library(zoo)
library(PNWColors)
library(patchwork)
library(ggsci)
detectCores()

# assign values ================================================================

growth_data_0 = c(0.2,mean(2.61,2.82,2.73,2.58,2.81,2.56,2.52,2.81,2.58,2.80,2.64,2.80,2.59,2.81,2.60,2.79,2.78,2.58,2.65,2.71,
                         3.78,2.79,2.74,2.69,2.66,2.53,2.68,2.62,2.63,2.74,2.76,2.86,2.82,2.64,2.73,2.65,2.49,2.71,2.88,2.57,
                         2.84,2.63,2.79,2.56,2.66,2.73,2.59,2.76,2.54,2.59,2.57,2.74,2.65,2.76,2.73,2.65,2.73,2.63,2.54, na.rm = TRUE),
                mean(3.71,3.58,3.71,3.66,3.53,3.65,3.66,3.54,3.69,3.54,3.61,3.65,3.61,3.53,3.38,3.55,3.50,3.59,NA,NA,3.58,3.63,3.63,
                     3.71,3.66,3.71,3.78,3.66,3.65,3.81,3.82,3.69,3.70,3.81,3.67,3.77,NA,NA,NA,NA,3.47,3.59,3.70,3.64,
                     3.49,3.54,3.51,3.53,3.80,3.65,3.70,3.68,3.64,3.55,3.73,3.60,3.64,3.60,NA, na.rm = TRUE),
                mean(3.82,4.01,3.75,4.00,3.91,3.88,
                     3.79,3.63,3.95,4.01,3.96,3.90,3.93,3.79,3.82,3.72,3.79,3.96,NA,NA,4.00,3.89,3.90,3.89,3.66,3.85,
                     3.92,3.81,3.87,3.96,3.90,4.05,3.85,3.94,3.85,3.83,3.82,3.91,3.65,NA,3.70,3.79,3.82,3.76,3.83,3.75,4.01,
                     3.86,3.96,3.85,3.98,3.82,3.90,3.93,3.90,3.97,3.85,3.79,3.76, na.rm = TRUE))

days = data.frame(days = c(1:21))
growth_data_0 = data.frame(ts = c(1,8,15,21), 
                         grow = growth_data_0)
growth_data_0 = left_join(days, growth_data_0, by = c("days" = "ts"))
growth_data_0$grow = na.approx(growth_data_0$grow) # interpolate NAs for stan

growth_data_1 = c(0.2,mean(2.55,2.65,2.67,2.70,2.68,2.63,2.64,2.70,2.61,2.73,2.66,2.72,2.69,2.48,2.69,2.71,2.65,2.66,2.70,2.68,
                            2.57,2.60,2.83,2.51,2.62,2.62,2.47,2.59,2.70,2.73,2.73,2.49,2.66,2.71,2.71,2.53,2.73,2.75,2.69,2.63,2.78,
                            2.89,2.68,2.61,2.61,2.60,2.75,2.75,2.80,2.62,2.82,2.61,2.74,2.63,2.65,2.71,2.80,2.74,NA, na.rm = TRUE),
                      mean(3.65,3.64,
                           3.53,3.68,3.61,3.50,3.49,3.59,3.61,3.50,3.56,3.53,3.66,3.57,3.54,3.68,3.59,NA,NA,NA,3.61,3.68,3.55,
                           3.57,3.42,3.60,3.58,3.52,3.57,3.66,3.50,3.15,3.50,3.52,3.55,3.59,3.68,3.64,3.46,3.56,3.53,3.43,3.60,3.54,
                           3.52,3.57,3.47,3.66,3.61,3.75,3.53,3.52,3.59,3.57,3.58,3.64,3.45,3.71,NA, na.rm = TRUE),
                  mean(3.69,3.76,3.97,3.95,3.50,
                       3.77,3.80,3.77,3.85,3.75,3.80,3.97,4.03,3.77,3.85,3.64,3.76,3.60,3.81,NA,3.87,3.70,3.88,3.81,3.71,3.53,
                       3.57,3.66,3.82,3.67,3.87,3.88,3.89,3.85,3.83,3.84,3.91,3.88,3.68,NA,3.70,3.53,3.71,3.62,3.63,3.53,3.57,
                       3.65,3.60,3.70,4.00,3.62,3.69,3.43,3.57,3.53,3.70,3.64,NA, na.rm = TRUE))
growth_data_1 = data.frame(ts = c(1,8,15,21), 
                         grow = growth_data_1)
growth_data_1 = left_join(days, growth_data_1, by = c("days" = "ts"))
growth_data_1$grow = na.approx(growth_data_1$grow) # interpolate NAs for stan

growth_data_2 = c(0.2,mean(2.57,2.58,2.60,2.62,2.39,2.49,2.71,2.70,2.48,2.67,2.65,2.56,2.50,2.48,2.35,2.67,2.57,2.62,2.51,2.54,
                           2.62,2.62,2.61,2.54,2.62,2.62,2.54,2.41,2.62,2.58,2.76,2.39,2.58,2.70,2.46,2.65,2.36,2.45,2.60,2.47,2.55,
                           2.47,2.46,2.66,2.60,2.66,2.45,2.52,1.63,2.61,2.61,2.54,2.36,NA,NA,NA,NA,NA,NA, na.rm = TRUE),
                  mean(3.22,3.44,3.46,3.53,
                       3.58,3.38,3.45,3.28,3.54,3.42,3.55,3.59,3.40,3.31,3.55,3.53,3.48,3.45,NA,NA,
                       3.54,3.50,3.41,3.55,3.59,3.32,3.48,3.48,3.49,3.39,3.51,3.46,3.40,3.54,3.46,3.30,3.41,3.45,NA,NA,
                       3.41,3.62,3.57,3.55,3.19,3.59,3.45,3.57,3.18,3.54,3.34,3.59,3.50,3.26,3.58,3.42,3.50,NA,NA, na.rm = TRUE),
                  mean(3.68,3.58,3.20,3.22,3.57,2.71,3.60,3.70,3.72,3.64,3.49,3.75,3.58,3.66,3.62,3.59,NA,NA,NA,NA,
                       3.62,3.74,3.55,3.67,3.67,3.53,3.29,3.83,3.67,3.66,3.51,3.78,3.62,3.53,3.37,3.37,NA,NA,NA,NA,
                       3.64,3.66,3.52,3.57,3.54,3.69,3.63,3.50,3.54,3.67,3.69,3.63,3.37,3.66,3.50,3.69,3.61,3.57,NA, na.rm = TRUE))
growth_data_2 = data.frame(ts = c(1,8,15,21), 
                           grow = growth_data_2)
growth_data_2 = left_join(days, growth_data_2, by = c("days" = "ts"))
growth_data_2$grow = na.approx(growth_data_2$grow) # interpolate NAs for stan

growth_data_3 = c(0.2, mean(2.35,2.39,2.36,2.31,2.31,2.23,2.38,2.36,2.29,2.14,2.41,2.34,2.36,2.28,2.23,2.18,NA,NA,NA,NA,2.26,
                            2.44,2.53,1.86,2.28,2.26,2.32,2.41,2.30,2.33,2.33,2.20,2.50,2.15,2.35,NA,NA,NA,NA,NA,2.02,2.20,2.19,
                            2.25,2.36,2.24,1.87,2.42,1.92,2.18,2.28,2.37,2.39,2.43,2.39,2.14,1.52,2.27,NA, na.rm = TRUE),
                  mean(2.82,2.73,2.74,2.60,2.75,
                       2.79,2.96,2.93,2.81,2.60,2.78,2.61,3.01,2.74,NA,NA,NA,NA,NA,NA,
                       2.80,2.86,2.77,2.77,2.95,2.75,2.77,2.45,2.83,2.69,2.51,2.74,2.74,2.90,2.87,NA,NA,NA,NA,NA,
                       2.66,2.66,2.81,2.86,2.78,2.77,2.67,2.85,2.70,2.77,2.98,NA,NA,NA,NA,NA,NA,NA,NA, na.rm = TRUE),
                  mean(2.88,2.98,3.09,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       2.99,2.90,2.72,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       2.90,2.94,2.78,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, na.rm = TRUE))
growth_data_3 = data.frame(ts = c(1,8,15,21), 
                           grow = growth_data_3)
growth_data_3 = left_join(days, growth_data_3, by = c("days" = "ts"))
growth_data_3$grow = na.approx(growth_data_3$grow) # interpolate NAs for stan

rep_0 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.79,6.32,7.10,6.63,7.58,8.75,7.84,8.68,8.75,16.51,20.26,19.10,
               19.56,21.95,20.35,22.62,21.95,20.40,28.68,32.21,32.40,31.62,34.68,33.95,33.56,36.42,33.95,40.45,
               49.11,47.40,45.62,52.32,49.80,48.73,53.16,50.15,54.56,65.95,60.83,63.51,71.32,67.57),
               nrow = 3, ncol = 21)
rep_0 = t(rep_0)
rep_0 = rowMeans(rep_0)
rep_1 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.63,6.58,6.63,7.63,7.37,7.32,7.63,7.37,7.32,19.58,18.00,18.43,
                  19.58,18.00,19.09,19.58,19.89,19.93,33.37,33.11,33.70,35.68,33.11,34.82,35.68,34.68,34.82,43.68,
                  49.47,50.15,47.47,49.47,51.09,47.47,51.53,51.09,59.26,65.42,62.70,67.21,69.26,66.04),
                 nrow = 3, ncol = 21)
rep_1 = t(rep_1)
rep_1 = rowMeans(rep_1)
rep_2 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.47,4.94,3.45,6.37,6.61,5.35,6.37,6.61,5.85,14.32,14.67,15.00,
                 15.58,16.61,15.05,16.37,17.28,16.30,26.74,27.61,27.00,29.63,31.17,28.45,30.32,31.94,29.70,40.43,
                 41.94,41.00,44.87,45.82,42.80,45.70,46.82,44.12,51.29,55.76,51.43,57.42,59.01,58.26),
               nrow = 3, ncol = 21)
rep_2 = t(rep_2)
rep_2 = rowMeans(rep_2)
rep_3 = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00,0.00,0.00,0.00,0.00,0.00,1.06,0.00,0.00,3.81,0.15,0.90,
                 6.00,0.70,1.40,6.86,0.70,1.78,8.08,0.70,1.78,10.15,0.70,2.78,10.15,0.70,2.78,10.15,1.50,3.28,10.15,1.90,
                 3.28,10.15,1.90,3.28,10.15,3.23,3.28,10.15,4.90,3.28),
               nrow = 3, ncol = 21)
rep_3 = t(rep_3)
rep_3 = rowMeans((rep_3))



ll_init = 0.2
cq_init = 0
l_y_obs_con = growth_data_0[,2]
l_y_obs_400 = growth_data_1[,2]
l_y_obs_2000 = growth_data_2[,2]
l_y_obs_10000 = growth_data_3[,2]
r_y_con = rep_0
r_y_400 = rep_1
r_y_2000 = rep_2
r_y_10000 = rep_3
ts = 1:21

billoir_data = list(
  ll_init = c(ll_init, cq_init),
  #ll_init_con = array(ll_init),
  #cq_init = array(cq_init),
  #cq_init = array(cq_init),
  ts = ts,
  l_y_obs_con = l_y_obs_con,
  l_y_obs_400 = l_y_obs_400,
  l_y_obs_2000 = l_y_obs_2000, 
  l_y_obs_10000 = l_y_obs_10000, 
  r_y_con = r_y_con, 
  r_y_400 = r_y_400, 
  r_y_2000 = r_y_2000,
  r_y_10000 = r_y_10000
)
# fit model ====================================================================
warmups = 10000
total_iterations = 20000
max_treedepth = 11
adapt_delta = 0.999
n_cores = 4
n_chains = 4

billoir_fit = 
  stan(file = here('./code/stan_files/billoir_replicate_attempt.stan'),
       data = billoir_data,
       chains = n_chains,
       cores = n_cores,
       warmup = warmups,
       iter = total_iterations,
       #seed = 1,
       refresh = 500,
       #verbose = TRUE,
       #open_progress = TRUE,
       control = list(adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
  ); beep(3)

# diagnose model problems ======================================================
plot(billoir_fit)
check_divergences(billoir_fit)

gen0_diagnostics = get_sampler_params(billoir_fit) %>% 
  set_names(1:4) %>% 
  map_df(as_tibble, .id = "chain") %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= warmups)

percent_divergent_rungs = gen0_diagnostics %>% 
  group_by(warmup, chain) %>% 
  summarize(percent_divergent = mean(divergent__ > 0)) %>% 
  ggplot() + 
  geom_col(aes(chain, percent_divergent, fill = warmup), position = "dodge",
           colour = "black") +
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs") +
  scale_fill_npg()
tree_depth = gen0_diagnostics %>% 
  ggplot(aes(iteration, treedepth__, colour = chain)) + 
  geom_line() +
  geom_hline(aes(yintercept = 10), colour = 'red') +
  scale_colour_locuszoom()
step_size = gen0_diagnostics %>% 
  ggplot(aes(iteration, stepsize__, colour = chain)) +
  geom_line() +
  scale_colour_locuszoom()

# look at output/make plots ====================================================
print(billoir_fit, 
      pars=c("theta_ll[1]", "theta_ll[2]", "theta_ll[3]", "theta_ll[4]",
             "Lp", "Rm", "Lm", "tau_l", "tau_r"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
parms=c("theta_ll[1]", "theta_ll[2]", "theta_ll[3]", "theta_ll[4]",
       "Lp", "Rm", "Lm", "tau_l", "tau_r")
gen_0_alltreat_onerep_fit_Trace = stan_trace(billoir_fit,parms)
gen_0_alltreat_onerep_fit_Dens = mcmc_dens(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Overlay = mcmc_dens_overlay(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Violin = mcmc_violin(gen_0_alltreat_onerep_fit,parms,
                                               probs = c(0.1, 0.5, 0.9))

gen_0_alltreat_output = rstan::extract(billoir_fit,
                                       permuted=TRUE,include=TRUE)

# make reproduction data
rep_con = data.frame(gen_0_alltreat_output$r_y_con_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
rep_400 = data.frame(gen_0_alltreat_output$r_y_400_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
rep_2000 = data.frame(gen_0_alltreat_output$r_y_2000_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
rep_10000 = data.frame(gen_0_alltreat_output$r_y_10000_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()

concs = c(rep(0, 22), rep(400, 22), rep(2000, 22), rep(10000, 22))
rep_data = c(r_y_con, r_y_400, r_y_2000, r_y_10000)

repro_plot_data = data.frame(concs, 
                             rep_data,
                             predictions = c(as.vector(rep_con[,1]), 
                                             as.vector(rep_400[,1]), 
                                             as.vector(rep_2000[,1]), 
                                             as.vector(rep_10000[,1])),
                             time = rep(seq(1,22,1),4))
repro_plot_data$concs = factor(as.character(repro_plot_data$concs))
repro_plot_data$concs = factor(repro_plot_data$concs, 
                               levels = c('0', '400', '2000', '10000'))
# length dataframe 
len_con = data.frame(gen_0_alltreat_output$l_y_con_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
len_400 = data.frame(gen_0_alltreat_output$l_y_400_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
len_2000 = data.frame(gen_0_alltreat_output$l_y_2000_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()
len_10000 = data.frame(gen_0_alltreat_output$l_y_10000_rep) %>% 
  summarize(across(`X1`:`X22`, mean)) %>% 
  t()

len_data = c(l_y_obs_con, l_y_obs_400, l_y_obs_2000, l_y_obs_10000)

len_plot_data = data.frame(concs, 
                           len_data,
                           predictions = c(as.vector(len_con[,1]), 
                                           as.vector(len_400[,1]), 
                                           as.vector(len_2000[,1]), 
                                           as.vector(len_10000[,1])),
                           time = rep(seq(1,22,1),4))
len_plot_data$concs = factor(as.character(len_plot_data$concs))
len_plot_data$concs = factor(lenro_plot_data$concs, 
                             levels = c('0', '400', '2000', '10000'))

all_treat_reproduction = ggplot(data = repro_plot_data) +
  geom_point(aes(x = time, y = rep_data, 
                 colour = as.factor(as.character(concs))),
             size = 1.7) +
  geom_line(aes(x = time, y = predictions, 
                colour = as.factor(as.character(concs))),
            size = 1.3,
            position = position_dodge(width = 2)) + 
  labs(x = 'Time (days)', y = 'Cumulative Reproduction (neonates)') +
  scale_colour_manual("Concentrations (particles/mL)", 
                      values = pnw_palette('Starfish', 4, type = 'discrete')[1:4]) +
  theme_bw()

all_treat_length = ggplot(data = len_plot_data) +
  geom_point(aes(x = time, y = len_data, 
                 colour = as.factor(as.character(concs))),
             size = 1.7) +
  geom_line(aes(x = time, y = predictions, 
                colour = as.factor(as.character(concs))),
            size = 1.3,
            position = position_dodge(width = 2)) + 
  labs(x = 'Time (days)', y = 'Length (mm)') +
  scale_colour_manual("Concentrations (particles/mL)", 
                      values = pnw_palette('Starfish', 4, type = 'discrete')[1:4]) +
  theme_bw() +
  theme(
    legend.position = 'none'
  )

rep_len_plot = all_treat_reproduction + all_treat_length

