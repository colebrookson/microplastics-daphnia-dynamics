########## 
##########
# This code contains the analysis for the fitting of the the individual component as a
# part of the larger analysis of the effect of microplastics on daphnia
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2020-08-25
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
detectCores()

source(here('./code/data_cleaning.R'))