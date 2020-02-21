##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 7.0 Probabilistic sensitivity analysis -- analysis of results
#### Copyright Isabelle Feldhaus
#### 27 January 2020
##################################


# 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(markovchain) 
library(readstata13)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')

## Load static parameters 
source('0.const_population.R') # Population characteristics
source('0.const_disability_weights.R') 
source('0.const_costs.R') 
source('1.costs.R')
source('1.model_states.R')
source('1.model_pop.R')
source('1.model_main.R')
source('1.outcomes.R')


# 1. Load results ---------------------------------------------------------

psa_results <- readr::read_csv('output/psa_results_1050L_2020-01-30_1356.csv')

ggplot(psa_results, aes(x = -incr_dm_dalys, y = incr_cost_total / n_hef_eligible, color = strategy)) + 
  geom_point() + 
  facet_grid(sex ~ strategy) + 
  theme_bw()

ggplot(psa_results, aes(x = -incr_che40, y = incr_cost_total / n_hef_eligible, color = strategy)) + 
  geom_point(alpha = 0.1) + 
  facet_grid(sex ~ strategy) + 
  theme_bw()


