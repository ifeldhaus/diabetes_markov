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
simulated_population <- readr::read_csv("output/ichar_200000L_2020-01-19_2041.csv")
simulated_population_results <- readr::read_csv('output/results_24034003L_2020-01-22_1646.csv') # 20% HEF threshold

# simulated_population <- readr::read_csv("output/ichar_200000L_2020-01-19_2306.csv")
# simulated_population_results <- readr::read_csv('output/results_24033414L_2020-01-22_1707.csv') # 30% HEF threshold

# Inflation factor based on simulated population size 
inf.fct <- 16e6 / nrow(simulated_population)

# Sample population IDs (with replacement)
bootstrapped_sample_size <- 10000
n_bootstraps <- 1000
bootstrapped_outcomes_total <- list()
bootstrapped_outcomes_sex <- list()

start_time <- Sys.time()
cat("Bootstrapping...")
for (bootstrap in 1:n_bootstraps) {
  
  # Tracking progress
  if (bootstrap %% 10 == 0) {
    cat('.')
  }
  
  # Bootstrap population IDs & select results
  bootstrapped_ids <- sample(simulated_population$id, bootstrapped_sample_size, replace = TRUE)
  bootstrapped_sample <- simulated_population[simulated_population$id %in% bootstrapped_ids, ]
  bootstrapped_sample_results <- simulated_population_results[simulated_population_results$id %in% bootstrapped_ids, ] # not working here
  
  # Compute outcomes of interest and store
  bootstrapped_outcomes_total[[bootstrap]] <- compute_outcomes(bootstrapped_sample, bootstrapped_sample_results, hef_threshold = 0.2) 
  bootstrapped_outcomes_sex[[bootstrap]] <- compute_outcomes_by_sex(bootstrapped_sample, bootstrapped_sample_results, hef_threshold = 0.2)
  
}

# Save results
save(bootstrapped_outcomes_total, file = filename_for_saving("bootstrapped_outcomes_total", nlines = 20, extension = "Rdata")) # 20L = 20% threshold
save(bootstrapped_outcomes_sex, file = filename_for_saving("bootstrapped_outcomes_sex", nlines = 20, extension = "Rdata"))

end_time <- Sys.time()
time <- end_time - start_time
cat("Done in", format(time), '\n')
