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
detectCores()


dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data/' #private data repo

growth_data = read_csv(paste0(dir, 'growth_data_gen0.csv'))
reproduction_data = read_csv(paste0(dir, 'offspring_data.csv'), 
                             guess_max = 15000)
unique(reproduction_data$treatment)
unique(growth_data$treatment)

# get the data required ========================================================

# keep only the LFC and generation 0
reproduction_data_con = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'LFC',
         generation == 0) 
reproduction_data_400 = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'PS_400',
         generation == 0) 
reproduction_data_2000 = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'PS_2000',
         generation == 0) 
reproduction_data_10000 = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'PS_10000',
         generation == 0) 

# make the length into mm
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment %in% c('LFC', 'MP400', 'MP2000', 'MP10000')) %>% #keep only ps2000 treatment
  filter(test == 1) %>% 
  filter(day < 23) %>%
  mutate(length_mm = length*0.001)  

# prep length data =============================================================

days = data.frame(ts = seq(1,22,1))
# create data matrix
l_y_obs_data_con = as.matrix(growth_data %>% 
                           filter(day < 22, 
                                  treatment == 'LFC') %>% 
                           select(day, length_mm, replicate) %>% 
                           pivot_wider(id_cols = c(replicate, day, length_mm),
                                       names_from = replicate,
                                       values_from = length_mm) %>% 
                           select_if(~ !any(is.na(.))) %>% 
                           select(-day))
l_y_obs_data_400 = as.matrix(growth_data %>% 
                               filter(day < 22, 
                                      treatment == 'MP400') %>% 
                               select(day, length_mm, replicate) %>% 
                               pivot_wider(id_cols = c(replicate, day, length_mm),
                                           names_from = replicate,
                                           values_from = length_mm) %>% 
                               select_if(~ !any(is.na(.))) %>% 
                               select(-day))
l_y_obs_data_2000 = as.matrix(growth_data %>% 
                               filter(day < 22, 
                                      treatment == 'MP2000') %>% 
                               select(day, length_mm, replicate) %>% 
                               pivot_wider(id_cols = c(replicate, day, length_mm),
                                           names_from = replicate,
                                           values_from = length_mm) %>% 
                               select_if(~ !any(is.na(.))) %>% 
                               select(-day))
l_y_obs_data_10000 = as.matrix(growth_data %>% 
                                filter(day < 22, 
                                       treatment == 'MP10000') %>% 
                                select(day, length_mm, replicate) %>% 
                                pivot_wider(id_cols = c(replicate, day, length_mm),
                                            names_from = replicate,
                                            values_from = length_mm) %>% 
                                select_if(~ !any(is.na(.))) %>% 
                                select(-day))
l_y_for_post_con = data.frame(growth_data %>% 
                            filter(day < 22,
                                   treatment == 'LFC') %>% 
                            select(day, length_mm, replicate) %>% 
                            pivot_wider(id_cols = c(replicate, day, length_mm),
                                        names_from = replicate,
                                        values_from = length_mm) %>% 
                            select_if(~ !any(is.na(.))))[,1:2]
l_y_for_post_con$day = l_y_for_post_con$day + 1
l_y_for_post_con = left_join(days, l_y_for_post_con, by = c('ts'= 'day'))
l_y_for_post_con$X1 = na.approx(l_y_for_post_con$X1)
l_y_for_post_400 = data.frame(growth_data %>% 
                            filter(day < 22,
                                   treatment == 'MP400') %>% 
                            select(day, length_mm, replicate) %>% 
                            pivot_wider(id_cols = c(replicate, day, length_mm),
                                        names_from = replicate,
                                        values_from = length_mm) %>% 
                            select_if(~ !any(is.na(.))))[,1:2]
l_y_for_post_400$day = l_y_for_post_400$day + 1
l_y_for_post_400 = left_join(days, l_y_for_post_400, by = c('ts'= 'day'))
l_y_for_post_400$X1 = na.approx(l_y_for_post_400$X1)
l_y_for_post_2000 = data.frame(growth_data %>% 
                                filter(day < 22,
                                       treatment == 'MP2000') %>% 
                                select(day, length_mm, replicate) %>% 
                                pivot_wider(id_cols = c(replicate, day, length_mm),
                                            names_from = replicate,
                                            values_from = length_mm) %>% 
                                select_if(~ !any(is.na(.))))[,1:2]
l_y_for_post_2000$day = l_y_for_post_2000$day + 1
l_y_for_post_2000 = left_join(days, l_y_for_post_2000, by = c('ts'= 'day'))
l_y_for_post_2000$X3 = na.approx(l_y_for_post_2000$X3)
l_y_for_post_10000 = data.frame(growth_data %>% 
                                 filter(day < 22,
                                        treatment == 'MP10000') %>% 
                                 select(day, length_mm, replicate) %>% 
                                 pivot_wider(id_cols = c(replicate, day, length_mm),
                                             names_from = replicate,
                                             values_from = length_mm) %>% 
                                 select_if(~ !any(is.na(.))))[,1:2]
l_y_for_post_10000$day = l_y_for_post_10000$day + 1
l_y_for_post_10000 = left_join(days, l_y_for_post_10000, by = c('ts'= 'day'))
l_y_for_post_10000$X7 = na.approx(l_y_for_post_10000$X7)

# prep reproduction data =======================================================

# identify replicates that didn't survive
reps_to_remove_con = reproduction_data_con %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')
unique(reps_to_remove_con$replicate)

# remove reps that didn't surviv
reproduction_data_x_removed_con = reproduction_data_con %>% 
  filter(!replicate %in% reps_to_remove_con$replicate)

# make numeric
reproduction_data_x_removed_con$offspring = 
  as.numeric(reproduction_data_x_removed_con$offspring)

# turn NA's into zeros
reproduction_data_x_removed_con$offspring[which(
  is.na(reproduction_data_x_removed_con$offspring))] = 0

# create data matrix
r_y_data_con = as.matrix(reproduction_data_x_removed_con %>% 
                       filter(day < 22) %>%  
                       pivot_wider(id_cols = c(replicate, day, offspring),
                                   names_from = replicate,
                                   values_from = offspring, 
                                   values_fn = max))
r_y_data_con = r_y_data_con[,2:ncol(r_y_data_con)]
r_y_data_cumul_con = matrix(nrow = nrow(r_y_data_con), 
                            ncol = ncol(r_y_data_con))
for(i in 1:nrow(r_y_data_cumul_con)) {
  for(j in 1:ncol(r_y_data_cumul_con)){
    
    r_y_data_cumul_con[i,j] = sum(r_y_data_con[1:i,j])
    
  }
}

# identify replicates that didn't survive
reps_to_remove_400 = reproduction_data_400 %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')
unique(reps_to_remove_400$replicate)

# remove reps that didn't surviv
reproduction_data_x_removed_400 = reproduction_data_400 %>% 
  filter(!replicate %in% reps_to_remove_400$replicate)

# make numeric
reproduction_data_x_removed_400$offspring = 
  as.numeric(reproduction_data_x_removed_400$offspring)

# turn NA's into zeros
reproduction_data_x_removed_400$offspring[which(
  is.na(reproduction_data_x_removed_400$offspring))] = 0

# create data matrix
r_y_data_400 = as.matrix(reproduction_data_x_removed_400 %>% 
                           filter(day < 22) %>%  
                           pivot_wider(id_cols = c(replicate, day, offspring),
                                       names_from = replicate,
                                       values_from = offspring, 
                                       values_fn = max))
r_y_data_400 = r_y_data_400[,2:ncol(r_y_data_400)]
r_y_data_cumul_400 = matrix(nrow = nrow(r_y_data_400), 
                            ncol = ncol(r_y_data_400))
for(i in 1:nrow(r_y_data_cumul_400)) {
  for(j in 1:ncol(r_y_data_cumul_400)){
    
    r_y_data_cumul_400[i,j] = sum(r_y_data_400[1:i,j])
    
  }
}

# identify replicates that didn't survive
reps_to_remove_2000 = reproduction_data_2000 %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')
unique(reps_to_remove_2000$replicate)

# remove reps that didn't surviv
reproduction_data_x_removed_2000 = reproduction_data_2000 %>% 
  filter(!replicate %in% reps_to_remove_2000$replicate)

# make numeric
reproduction_data_x_removed_2000$offspring = 
  as.numeric(reproduction_data_x_removed_2000$offspring)

# turn NA's into zeros
reproduction_data_x_removed_2000$offspring[which(
  is.na(reproduction_data_x_removed_2000$offspring))] = 0

# create data matrix
r_y_data_2000 = as.matrix(reproduction_data_x_removed_2000 %>% 
                           filter(day < 22) %>%  
                           pivot_wider(id_cols = c(replicate, day, offspring),
                                       names_from = replicate,
                                       values_from = offspring, 
                                       values_fn = max))
r_y_data_2000 = r_y_data_2000[,2:ncol(r_y_data_2000)]
r_y_data_cumul_2000 = matrix(nrow = nrow(r_y_data_2000), 
                            ncol = ncol(r_y_data_2000))
for(i in 1:nrow(r_y_data_cumul_2000)) {
  for(j in 1:ncol(r_y_data_cumul_2000)){
    
    r_y_data_cumul_2000[i,j] = sum(r_y_data_2000[1:i,j])
    
  }
}

# identify replicates that didn't survive
reps_to_remove_10000 = reproduction_data_10000 %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')
unique(reps_to_remove_10000$replicate)

# remove reps that didn't surviv
reproduction_data_x_removed_10000 = reproduction_data_10000 %>% 
  filter(!replicate %in% reps_to_remove_10000$replicate)

# make numeric
reproduction_data_x_removed_10000$offspring = 
  as.numeric(reproduction_data_x_removed_10000$offspring)

# turn NA's into zeros
reproduction_data_x_removed_10000$offspring[which(
  is.na(reproduction_data_x_removed_10000$offspring))] = 0

# create data matrix
r_y_data_10000 = as.matrix(reproduction_data_x_removed_10000 %>% 
                            filter(day < 22) %>%  
                            pivot_wider(id_cols = c(replicate, day, offspring),
                                        names_from = replicate,
                                        values_from = offspring, 
                                        values_fn = max))
r_y_data_10000 = r_y_data_10000[,2:ncol(r_y_data_10000)]


# assign values ================================================================

N_obs = 10
N_mis = 12
ii_obs = c(1, 3, 6, 8, 10, 13, 15, 17, 20, 22)
ii_mis = c(2, 4, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21) 
ll_init = 0.2
cq_init = 0.1
l_y_obs_con = l_y_for_post_con[,2]
l_y_obs_400 = l_y_for_post_400[,2]
l_y_obs_2000 = l_y_for_post_2000[,2]
l_y_obs_10000 = l_y_for_post_10000[,2]
r_y_con = rowMeans(r_y_data_cumul_con)
r_y_400 = rowMeans(r_y_data_cumul_400)
r_y_2000 = rowMeans(r_y_data_cumul_2000)
r_y_10000 = r_y_data_10000
ts = 1:nrow(r_y_data)

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
gen_0_alltreat_onerep_fit = 
  stan(file = here('./code/stan_files/organism_costs_model_gen0_alltreat.stan'),
       data = gen_0_alltreat_data,
       chains = 4,
       cores = 8,
       warmup = 2000,
       iter = 5000,
       seed = 1,
       verbose = TRUE,
       #open_progress = TRUE,
       control = list(adapt_delta = 0.999)); beep(3)
saveRDS(gen_0_alltreat_onerep_fit, 
        here('/output/intermediate-objects/gen_0_alltreat_onerep_fit.RDS'))
print(gen_0_alltreat_onerep_fit, 
      pars=c("theta_ll[1]", "theta_cq[1]", "cstar", "NEC",
             "Lp", "Rm", "Lm", "tau_l", "tau_r"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
gen_0_alltreat_onerep_fit_Trace = stan_trace(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Dens = mcmc_dens(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Overlay = mcmc_dens_overlay(gen_0_alltreat_onerep_fit,parms)
gen_0_alltreat_onerep_fit_Violin = mcmc_violin(gen_0_alltreat_onerep_fit,parms,
                                             probs = c(0.1, 0.5, 0.9))



# make reproduction plot
gen_0_alltreat_output = rstan::extract(gen_0_alltreat_onerep_fit,
                                     permuted=TRUE,include=TRUE)

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

ggplot(data = repro_plot_data) +
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




#print("the z thing length: ", z_ll_400[i,1]);
#print("L1[i]: ", L1[i]);
#print("L2[i]: ", L2[i]);
#print("L3[i]: ", L3[i]);
