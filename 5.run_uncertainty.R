##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 5.0 Uncertainty analysis
#### Copyright Isabelle Feldhaus
#### 21 January 2020
##################################

# 0. Setup ----------------------------------------------------------------

setwd("/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/diabetes_markov/")

OUTPUT_DIR = 'output'

library(tidyverse)
library(markovchain) 
library(readstata13)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')
source('1.outcomes.R')


# 1. Bootstrap simulated population -------------------------------------------------

# Load population 
simulated_population <- readr::read_csv("output/ichar_500000L_2020-01-23_1757.csv")
simulated_population_results <- readr::read_csv('output/results_60000237L_2020-01-23_1849.csv') # 20% HEF threshold

n_population <- nrow(simulated_population)

# simulated_population <- readr::read_csv("output/ichar_200000L_2020-01-19_2306.csv")
# simulated_population_results <- readr::read_csv('output/results_24033414L_2020-01-22_1707.csv') # 30% HEF threshold

# Inflation factor based on simulated population size 
inf.fct <- 16e6 / n_population

# Sample population IDs (with replacement)
# bootstrapped_sample_size <- 200000
bootstrapped_sample_size <- 2e5
n_bootstraps <- 1000

# Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(4)

start_time <- Sys.time()
cat("Bootstrapping...")

nested_simulation_results_pop_per_id <- merge(simulated_population_results, simulated_population, by = 'id') %>% group_by(id) %>% nest
rm(simulated_population_results)
gc()

bootstrap_results <- foreach::foreach(bootstrap = 1:n_bootstraps) %dopar% {
  
  # Tracking progress 
  cat('.')
  
  # Bootstrap population IDs & select results
  bootstrapped_ids <- sample(1:n_population, bootstrapped_sample_size, replace = TRUE)
  bootstrapped_pop <- simulated_population[bootstrapped_ids, ]
  bootstrapped_results_pop <- nested_simulation_results_pop_per_id[bootstrapped_ids, ] %>% unnest()
  
  # Compute outcomes of interest and store
  return(compute_outcomes(bootstrapped_pop, bootstrapped_results_pop, hef_threshold = 0.2))
  
}

time_elapsed <- Sys.time() - start_time

# Combine results
bootstrapped_outcomes_total <- do.call(rbind, bootstrap_results)

# Save results
save(bootstrapped_outcomes_total, file = filename_for_saving("bootstrapped_outcomes_total", nlines = 20, extension = "Rdata")) # 20L = 20% threshold

cat("Done in", format(time_elapsed), '\n')
