##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 1.27 Markov model 
#### Copyright Isabelle Feldhaus
#### 19 January 2020
##################################

# 0. Setup ----------------------------------------------------------------

setwd("/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/diabetes_markov/")

OUTPUT_DIR = 'output'

library(tidyverse)
library(markovchain)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')


# I. Hypotheses -----------------------------------------------------------

source('0.const_population.R')


## 1. Disability weights ---------------------------------------------------

source('0.const_disability_weights.R')


## 3. Costs ----------------------------------------------------------------

source('0.const_costs.R')


## 4. Main parameters --------------------------------------------------

# Total population
TOTAL_SIMULATION_POP_SIZE <- 8e5

## HEF eligibility and enrollment
BASE_HEF_QUANTILE <- 0.3
hef_quantile <- 0.2

## Benefits of coverage
coverage <- 1 # 0.8 # If individual is under HEF, 80% of expenditures are covered.


# II. Simulate populations ----------------------------------------------------

source('1.model_states.R')
source('1.model_pop.R')

## Settings
n_population <- TOTAL_SIMULATION_POP_SIZE * BASE_HEF_QUANTILE
n_cycles = 45
create_new_population <- FALSE

if (create_new_population) {
  
  cat('Simulating population of', n_population, 'people')
  
  all_income <- rgamma(n = n_population / BASE_HEF_QUANTILE, shape = 0.5, scale = avg_income)
  income_vec <- all_income[all_income < quantile(all_income, BASE_HEF_QUANTILE)]
  
  ichar <- tibble(
    id = 1:n_population,
    age_start = sample(age_dist, n_population, replace = TRUE),
    sex = ifelse(rbinom(n = n_population, size = 1, p = p_female) == 0, "Male", "Female"), # Relabel for `look_up`
    income = income_vec, # Distribution derived from average income per month 2017, CSES (1960 thousand KHR / 481.93 USD 2019)
    disposable_income = income - subsistence,
    init_state = mapply(function(age, sex) sample(stateNames, size = 1, prob = initialStates_KHR(age, sex)), age_start, sex)
  ) %>% as.data.frame
  
  ## Matrix of state probabilities (pointers) for each individual and cycle
  person_cycle_x <- matrix(runif(n_population * (n_cycles * 2 + 2)), n_population, (n_cycles * 2 + 2))

  save_df_to_csv(ichar, "ichar")
  save_object_to_rdata(person_cycle_x)
  
  cat(' Done\n')
  
  } else {
    cat('Reading in existing population...\n')
    dir("output", sprintf("(ichar|person_cycle_x)_%dL.*", n_population))
    ichar <- data.frame(readr::read_csv("output/ichar_240000L_2020-01-30_0833.csv"))
    person_cycle_x <- load_object_from_rdata('output/person_cycle_x_240000L_2020-01-30_0833.RData')
}

population <- add_hef_to_pop(ichar, person_cycle_x, hef_quantile)


# Run Model ----------------------------------------------------------

source('1.model_main.R')

## Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(8)

all_results_final <- run_model(population, person_cycle_x)

## Export flat files
save_df_to_csv(all_results_final, "results")

