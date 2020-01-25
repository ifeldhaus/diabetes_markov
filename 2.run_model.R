##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 1.27 Markov model 
#### Copyright Isabelle Feldhaus
#### 19 January 2020
##################################

# 0. Setup ----------------------------------------------------------------

setwd("/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/diabetes_markov/")

OUTPUT_DIR = 'output'

# Install packages via SSH connection
# chmod 400 /Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/code/diabetes-markov.pem
# ssh -v -i "/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/code/diabetes-markov.pem" ubuntu@ec2-18-216-81-155.us-east-2.compute.amazonaws.com
# sudo -s
# sudo R -e "install.packages('heemod', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('markovchain', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('readstata13', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('foreach', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('doParallel', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('doMC', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('doSNOW', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('wesanderson', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('ggsci', repos='http://cran.us.r-project.org')"
# sudo R -e "install.packages('xtable', repos='http://cran.us.r-project.org')"

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


## 4. Coverage parameters --------------------------------------------------

## HEF eligibility and enrollment
hef_quantile <- 0.2
p_hef_enrollment <- 0.75

## Benefits of coverage
coverage <- 1 # 0.8 # If individual is under HEF, 80% of expenditures are covered.

## HEF utilization at point of care
p_hef_utilization = 1 # 0.16


## 6. States and strategies ------------------------------------------------

## Strategies to simulate
strategyNames <- c("base", "screen_only", "tx_only", "comp_only", "screen_tx", "tx_comp", "screen_tx_comp")


# II. Simulate populations ----------------------------------------------------

source('1.model_states.R')
source('1.model_pop.R')

## Settings
# n_population <- as.integer(5e5) # total: 16e6
n_population <- 1000
n_cycles = 45
create_new_population <- TRUE

if (create_new_population) {
  
  cat('Simulating population of', n_population, 'people')
  ichar <- tibble(
    id = 1:n_population,
    age_start = sample(age_dist, n_population, replace = TRUE),
    sex = ifelse(rbinom(n = n_population, size = 1, p = p_female) == 0, "Male", "Female"), # Relabel for `look_up`
    income = rgamma(n = n_population, shape = 0.5, scale = avg_income), # Distribution derived from average income per month 2017, CSES (1960 thousand KHR / 481.93 USD 2019)
    disposable_income = income - subsistence,
    init_state = mapply(function(age, sex) sample(stateNames, size = 1, prob = initialStates_KHR(age, sex)), age_start, sex)
  ) %>% as.data.frame
  
  ## Matrix of state probabilities (pointers) for each individual and cycle
  person_cycle_x <- matrix(runif(n_population * (n_cycles + 2)), n_population, (n_cycles + 2))

  save_df_to_csv(ichar, "ichar")
  save_object_to_rdata(person_cycle_x)
  
  cat(' Done\n')
  
  } else {
    cat('Reading in existing population...\n')
    dir("output", sprintf("(ichar|person_cycle_x)_%dL.*", n_population))
    ichar <- data.frame(readr::read_csv("output/ichar_500000L_2020-01-23_1757.csv"))
    person_cycle_x <- load_object_from_rdata('output/person_cycle_x_500000L_2020-01-23_1757.RData')
}

population <- add_hef_to_pop(ichar, person_cycle_x, hef_quantile)


# III. Run Model ----------------------------------------------------------

# Load required functions 
source("1.costs.R")
source("1.dalys.R")
source("1.FRP.R")
source("1.model_tmat_construction.R")

## Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(8)

# Determine sample groups for parallel execution, change first line only
indicative_parallel_samples <- 20
sample_size <- (n_population - 1) %/% indicative_parallel_samples + 1
n_parallel_samples <- (n_population - 1) %/% sample_size + 1

# Empty lists to store final results
all_results <- list()

## Simulation for each strategy
for (strategy in strategyNames) {
  
  cat("Starting strategy:", strategy)
  
  ## Prepare iteration
  # Set seed and time computation
  start_time <- Sys.time()
  
  cat(sprintf(', running %d samples of size %d... ', n_parallel_samples, sample_size))
  
  ## Perform iteration for each iteration group
  results_per_sample <- foreach::foreach(i_sample = 1:n_parallel_samples) %dopar% {
      
    ## Instantiate the list of vectors to store the results
    n_rows <- sample_size * n_cycles
    results <- list(
      strategy = character(n_rows),
      id = integer(n_rows),
      cycle = integer(n_rows),
      state = character(n_rows),
      cost_diagnosis = numeric(n_rows),
      cost_tx = numeric(n_rows),
      cost_complications = numeric(n_rows),
      cost_total = numeric(n_rows),
      dalys_dm = numeric(n_rows),
      dalys_other = numeric(n_rows),
      dalys_total = numeric(n_rows),
      oop_medical = numeric(n_rows),
      oop_nonmedical = numeric(n_rows),
      oop_total = numeric(n_rows),
      che10_result = integer(n_rows),
      che25_result = integer(n_rows),
      che40_result = integer(n_rows),
      che10_disposable_result = integer(n_rows),
      che25_disposable_result = integer(n_rows),
      che40_disposable_result = integer(n_rows),
      pov_result = integer(n_rows),
      pov_disposable_result = integer(n_rows)
    )
    
    ## Set counter for the sample, iterating along people and cycles
    j_counter <- 1
    
    ## Select population indices
    population_indices <- ((i_sample-1)*sample_size+1) : min(((i_sample)*sample_size), n_population)
    for (i_person in population_indices) {
      
      if (i_person %% 1000 == 0) {
        cat('.')
      }

      id <- population[i_person, 'id']
      age <- as.numeric(population[i_person, 'age_start'])
      sex <- population[i_person, 'sex']
      income <- as.numeric(population[i_person, 'income'])
      disposable_income <- as.numeric(population[i_person, 'disposable_income'])
      hef <- population[i_person, 'hef']
      hef_utilization <- population[i_person, 'hef_utilization']
      init_state <- population[i_person, 'init_state']
      
      #  Impact of strategy on Utilization:
      ## This only impacts cost computations
      is_strategy_comp <- strategy == "comp_only" | strategy == "tx_comp" | strategy == "screen_tx_comp"
      is_strategy_screen <- strategy == "screen_only" | strategy == "screen_tx" | strategy == "screen_tx_comp"
      e_strategy_util <- ifelse(is_strategy_screen & hef, e_screen_on_util, 1) * ifelse(is_strategy_comp & hef, e_comp_on_util, 1)
      p_op <- ifelse(hef, 0.172, 0.117) * e_strategy_util
      p_hosp <- 0.015 * e_strategy_util
      
      
      # Cycle through time horizon
      for (cycle in 1:n_cycles) {
        
        tmat <- build_tmat(strategy = strategy, hef = hef, age = age, sex = sex)
        
        # Advance one step in the chain
        current_state <- as.numeric(stateNames == init_state)
        state_probs <- current_state %*% tmat
        
        # Check the validity of the transition matrix
        if (cycle == 20 &  i_person %% 1000 == 0) {
          mc <- new("markovchain",
                    transitionMatrix = tmat)
        }
        
        # Monte Carlo simulation to determine end-cycle state
        end_state <- colnames(state_probs)[person_cycle_x[i_person, cycle] < cumsum(state_probs)][1]
        
        # Compute costs
        costs <- compute_costs(init_state, end_state, p_provider,
                               outpatient_costs, hospitalization_costs,
                               p_op, p_hosp, hef, hef_utilization, dr, cycle)
        
        # Compute DALYs
        dalys <- compute_dalys(end_state, age, sex)
        
        # Compute CHE
        che_result <- che(costs$oop_total, income)
        che_result_disposable <- che(costs$oop_total, disposable_income)
        
        # Compute POV
        pov_result <- pov(costs$oop_total, income)
        pov_result_disposable <- pov(costs$oop_total, disposable_income)
        
        # Store results in data frame
        results[[1]][j_counter] <- strategy
        results[[2]][j_counter] <- id
        results[[3]][j_counter] <- cycle
        results[[4]][j_counter] <- end_state
        results[[5]][j_counter] <- costs$cost_diagnosis
        results[[6]][j_counter] <- costs$cost_tx
        results[[7]][j_counter] <- costs$cost_complications
        results[[8]][j_counter] <- costs$cost_total
        results[[9]][j_counter] <- dalys$dalys_dm
        results[[10]][j_counter] <- dalys$dalys_other
        results[[11]][j_counter] <- dalys$dalys_total
        results[[12]][j_counter] <- costs$oop_medical
        results[[13]][j_counter] <- costs$oop_nonmedical
        results[[14]][j_counter] <- costs$oop_total
        results[[15]][j_counter] <- che_result$che10
        results[[16]][j_counter] <- che_result$che25
        results[[17]][j_counter] <- che_result$che40
        results[[18]][j_counter] <- che_result_disposable$che10
        results[[19]][j_counter] <- che_result_disposable$che25
        results[[20]][j_counter] <- che_result_disposable$che40
        results[[21]][j_counter] <- pov_result
        results[[22]][j_counter] <- pov_result_disposable
        
        # Re-assign `initial` state for next cycle
        init_state <- end_state
        
        # Re-assign `age` for next cycle
        age <- age + 1
        
        # Update counter
        j_counter <- j_counter + 1
        
        # Break out of loop under conditions
        if (age == 70) {
          break
        }
        
        if (end_state == "dm_death") {
          break
        }
        
        if (end_state == "other_death") {
          break
        } 
      }
    }
    return(results %>% as.tibble %>% filter(strategy != ''))
  }
  
  # End time computation
  end_time <- Sys.time()
  time <- end_time - start_time
  cat("Done in", format(time), '\n')
  
  results <- do.call(rbind, results_per_sample)
  all_results[[strategy]] <- results
}

# Combine results -- may be too large
all_results_final <- do.call(rbind, all_results) %>% as.data.frame

## Profiling
# Rprof(NULL)
# profvis::profvis(prof_input = 'profiling_model')

## Export flat files
save_df_to_csv(all_results_final, "results")

# Stop cluster
# parallel::stopCluster(cl)

