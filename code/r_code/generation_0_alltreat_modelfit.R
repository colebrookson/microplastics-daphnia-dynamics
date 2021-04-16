########## 
##########
# This code contains some data preparation and model fitting for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-04-06
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
 
ll_init = 0.2
cq_init = 0
l_y_obs_con = l_y_for_post_con[,2]
l_y_obs_400 = l_y_for_post_400[,2]
l_y_obs_2000 = l_y_for_post_2000[,2]
l_y_obs_10000 = l_y_for_post_10000[,2]
r_y_con = rowMeans(r_y_data_cumul_con)
r_y_400 = rowMeans(r_y_data_cumul_400)
r_y_2000 = rowMeans(r_y_data_cumul_2000)
r_y_10000 = r_y_data_10000
ts = 1:nrow(r_y_data_cumul_2000)

gen_0_alltreat_data = list(
  ll_init = array(ll_init),
  cq_init = array(cq_init),
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
warmups = 2000
total_iterations = 5000
max_treedepth = 10
adapt_delta = 0.99
n_cores = 4
n_chains = 4

gen_0_alltreat_onerep_fit = 
  stan(file = here('./code/stan_files/organism_costs_model_gen0_alltreat.stan'),
       data = gen_0_alltreat_data,
       chains = n_chains,
       cores = n_cores,
       warmup = warmups,
       iter = total_iterations,
       seed = 1,
       refresh = 500,
       #verbose = TRUE,
       # init = list(
       #   list(`theta_ll[1]` = 0.2,
       #        `theta_cq[1]` = 1.5,
       #        cstar = 5000,
       #        NEC = 2000,
       #        Lp = 0.9,
       #        Rm = 10,
       #        Lm = 4
       #        ),
       #   list(`theta_ll[1]` = 0.02,
       #        `theta_cq[1]` = 2,
       #        cstar = 6000,
       #        NEC = 1500,
       #        Lp = 0.2,
       #        Rm = 5,
       #        Lm = 2
       #        ),
       #   list(`theta_ll[1]` = 0.1,
       #        `theta_cq[1]` = 0.6,
       #        cstar = 7000,
       #        NEC = 6000,
       #        Lp = 0.1,
       #        Rm = 15,
       #        Lm = 1
       #        ),
       #   list(`theta_ll[1]` = 0.9,
       #        `theta_cq[1]` = 5,
       #        cstar = 4000,
       #        NEC = 3000,
       #        Lp = 1.6,
       #        Rm = 19,
       #        Lm = 41
       #        )
       # ),
       #open_progress = TRUE,
       control = list(adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
       ); beep(3)
saveRDS(gen_0_alltreat_onerep_fit, 
        here('/output/intermediate-objects/gen_0_alltreat_onerep_fit.RDS'))
parms = c("theta_ll[1]", 
          "theta_cq[1]", 
          "cstar", "NEC",
          "Lp", "Rm", "Lm", "tau_l", "tau_r")
# diagnose model problems ======================================================
plot(gen_0_alltreat_onerep_fit)
check_divergences(gen_0_alltreat_onerep_fit)

gen0_diagnostics = get_sampler_params(gen_0_alltreat_onerep_fit) %>% 
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
print(gen_0_alltreat_onerep_fit, 
      pars=c("theta_ll[1]", "theta_cq[1]", "cstar", "NEC",
             "Lp", "Rm", "Lm", "tau_l", "tau_r"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
gen_0_alltreat_onerep_fit_Trace = stan_trace(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Dens = mcmc_dens(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Overlay = mcmc_dens_overlay(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Violin = mcmc_violin(gen_0_alltreat_onerep_fit,parms,
                                             probs = c(0.1, 0.5, 0.9))

gen_0_alltreat_output = rstan::extract(gen_0_alltreat_onerep_fit,
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


#print("the z thing length: ", z_ll_400[i,1]);
#print("L1[i]: ", L1[i]);
#print("L2[i]: ", L2[i]);
#print("L3[i]: ", L3[i]);
