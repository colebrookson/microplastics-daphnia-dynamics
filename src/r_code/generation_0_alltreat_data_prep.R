########## 
##########
# This code contains some data preparation for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-04-13
##########
##########

# set-up =======================================================================

library(here)
library(tidyverse)
library(zoo)

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

# make the length into mm, keep the appropriate tests and days, and keep only 
# treatments we want to use
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment %in% c('LFC', 'MP400', 
                          'MP2000', 'MP10000')) %>% 
  filter(test == 1) %>% 
  filter(day < 23) %>%
  mutate(length_mm = length*0.001)  

# prep length data =============================================================

days = data.frame(ts = seq(1,22,1))
# create data matrix
l_y_obs_data_con = as.matrix(growth_data %>% 
                               filter(day < 22, # keep days below 22
                                      treatment == 'LFC') %>% # keep treatment
                               select( # keep only the data we want to group by
                                 day, length_mm, replicate) %>% 
                               pivot_wider( # turn the data into wide format
                                 id_cols = c(replicate, day, 
                                                       length_mm),
                                           names_from = replicate,
                                           values_from = length_mm) %>% 
                               select_if( # get rid of any NA values
                                 ~ !any(is.na(.))) %>% 
                               select(-day))
l_y_obs_data_400 = as.matrix(growth_data %>% 
                               filter(day < 22, 
                                      treatment == 'MP400') %>% 
                               select(day, length_mm, replicate) %>% 
                               pivot_wider(id_cols = c(replicate, day, 
                                                       length_mm),
                                           names_from = replicate,
                                           values_from = length_mm) %>% 
                               select_if(~ !any(is.na(.))) %>% 
                               select(-day))
l_y_obs_data_2000 = as.matrix(growth_data %>% 
                                filter(day < 22, 
                                       treatment == 'MP2000') %>% 
                                select(day, length_mm, replicate) %>% 
                                pivot_wider(id_cols = c(replicate, day, 
                                                        length_mm),
                                            names_from = replicate,
                                            values_from = length_mm) %>% 
                                select_if(~ !any(is.na(.))) %>% 
                                select(-day))
l_y_obs_data_10000 = as.matrix(growth_data %>% 
                                 filter(day < 22, 
                                        treatment == 'MP10000') %>% 
                                 select(day, length_mm, replicate) %>% 
                                 pivot_wider(id_cols = c(replicate, day, 
                                                         length_mm),
                                             names_from = replicate,
                                             values_from = length_mm) %>% 
                                 select_if(~ !any(is.na(.))) %>% 
                                 select(-day))
l_y_for_post_con = data.frame(growth_data %>% 
                                filter(day < 22,
                                       treatment == 'LFC') %>% 
                                select(day, length_mm, replicate) %>% 
                                pivot_wider(id_cols = c(replicate, day, 
                                                        length_mm),
                                            names_from = replicate,
                                            values_from = length_mm) %>% 
                                select_if(~ !any(is.na(.))))[,1:2]
l_y_for_post_con$day = l_y_for_post_con$day + 1 # change days from 0-21 to 1-22
l_y_for_post_con = left_join( # join up days so we keep the NAs
  days, l_y_for_post_con, by = c('ts'= 'day'))
l_y_for_post_con$X1 = na.approx(l_y_for_post_con$X1) # interpolate NAs for stan
l_y_for_post_400 = data.frame(growth_data %>% 
                                filter(day < 22,
                                       treatment == 'MP400') %>% 
                                select(day, length_mm, replicate) %>% 
                                pivot_wider(id_cols = c(replicate, day, 
                                                        length_mm),
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
                                 pivot_wider(id_cols = c(replicate, day, 
                                                         length_mm),
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
                                  pivot_wider(id_cols = c(replicate, day, 
                                                          length_mm),
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


