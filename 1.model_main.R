

## Strategies to simulate
strategyNames <- c("base", "screen_only", "tx_only", "comp_only", "screen_tx", "tx_comp", "screen_tx_comp")

# Load required functions 
source("1.costs.R")
source("1.dalys.R")
source("1.FRP.R")
source("1.model_tvec_construction.R")  # Includes definition of parameters

# Probabilities that we don't know where to put them anymore and it's too late to think about that.
p_op_hef <- 0.172
p_op_non_hef <- 0.117
p_hosp_base <- 0.015


run_model <- function(
  population,
  person_cycle_x
){
  
  n_population <- nrow(population)
  
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
    # results_per_sample <- lapply(1:n_parallel_samples, function(i_sample) {
      
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
        current_state <- population[i_person, 'init_state']
        
        #  Impact of strategy on Utilization:
        ## This only impacts cost computations
        is_strategy_comp <- strategy == "comp_only" | strategy == "tx_comp" | strategy == "screen_tx_comp"
        is_strategy_screen <- strategy == "screen_only" | strategy == "screen_tx" | strategy == "screen_tx_comp"
        e_strategy_util <- ifelse(is_strategy_screen & hef, e_screen_on_util, 1) * ifelse(is_strategy_comp & hef, e_comp_on_util, 1)
        p_op <- ifelse(hef, p_op_hef, p_op_non_hef) * e_strategy_util
        p_hosp <- p_hosp_base * e_strategy_util
        
        # Cycle through time horizon
        for (i_cycle in 1:n_cycles) {
          
          state_probs_base <- build_tvec(current_state = current_state, strategy = 'base', hef = hef, age = age, sex = sex)
          
          # Monte Carlo simulation to determine end-cycle state
          next_state <- draw_state_from_x(state_probs_base, person_cycle_x[i_person, i_cycle])
          next_state_vec <- stateNames == next_state
          
          if (strategy != 'base') {
            state_probs <- build_tvec(current_state = current_state, strategy = strategy, hef = hef, age = age, sex = sex)
            
            dp <- state_probs - state_probs_base
            
            if (dp[next_state_vec] < 0) {
              
              dp_plus <- (abs(dp) + dp) / 2
              p_to_go <- dp_plus / sum(dp_plus)  # Normalize by all the states you might go to
              
              dp_minus <- (abs(dp) - dp) / 2
              p_to_leave <- dp_minus / state_probs_base  # Normalize with the initial base prob to go to the state
              
              new_p <- next_state_vec * (1 - p_to_leave[next_state_vec]) + p_to_go * p_to_leave[next_state_vec]
              
              next_state <- draw_state_from_x(new_p, person_cycle_x[i_person, i_cycle * 2])
            }
          }
          
          # Compute costs
          costs <- compute_costs(strategy, current_state, next_state, p_provider,
                                 outpatient_costs, hospitalization_costs,
                                 p_op, p_hosp, hef, hef_utilization, dr, i_cycle)
          
          # Compute DALYs
          dalys <- compute_dalys(next_state, age, sex)
          
          # Compute CHE
          che_result <- che(costs$oop_total, income)
          che_result_disposable <- che(costs$oop_total, disposable_income)
          
          # Compute POV
          pov_result <- pov(costs$oop_total, income)
          pov_result_disposable <- pov(costs$oop_total, disposable_income)
          
          # Store results in data frame
          results[[1]][j_counter] <- strategy
          results[[2]][j_counter] <- id
          results[[3]][j_counter] <- i_cycle
          results[[4]][j_counter] <- next_state
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
          current_state <- next_state
          
          # Re-assign `age` for next cycle
          age <- age + 1
          
          # Update counter
          j_counter <- j_counter + 1
          
          # Break out of loop under conditions
          if (age == 70) {
            break
          }
          
          if (next_state == "dm_death") {
            break
          }
          
          if (next_state == "other_death") {
            break
          } 
        }
      }
      return(results %>% as.tibble %>% filter(strategy != ''))
    }
    # })
    
    # End time computation
    end_time <- Sys.time()
    time <- end_time - start_time
    cat("Done in", format(time), '\n')
    
    results <- do.call(rbind, results_per_sample)
    all_results[[strategy]] <- results
  }
  
  # Combine results -- may be too large
  return(do.call(rbind, all_results) %>% as.data.frame)
}