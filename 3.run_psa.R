##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 4.3 Probabilistic sensitivity analysis
#### Copyright Isabelle Feldhaus
#### 21 January 2020
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

# 1. Parameter Distributions ---------------------------------------------

# Define parameter distributions in dataframe
n_size = 10

e_oad <- function() rbeta(n_size, 18.98899,  8.935997)

modified_parameters_matrix <- tibble(
  
  # e_oad = rnorm(n_size, 0.50, 0.10),  ?->? rbeta(n_size, 26.00386, 8.667955), # Herman et al., 2017 (avg 0.75)  
  
  p_oad = rbeta(n_size, 15.35063, 53.17897),
  p_ins = rbeta(n_size, 11.34648, 656.0935),
  p_diet = (1 - p_oad - p_ins),
  
  p_nephro_to_retino = p_retino * rnorm(n_size, 5.01, 0.18),
  
  # Effects of therapies (RRs)
  e_diet = rnorm(n_size, 0.90, 0.05),
  
  e_oad_mi = rbeta(n_size, 9.068628, 5.797975),
  e_oad_stroke = rbeta(n_size, 3.760544, 2.613259),
  e_oad_nephro = 1 - rbeta(n_size, 3.987394, 1.668484),
  e_oad_retino = e_oad(),
  e_oad_neuro = rbeta(n_size, 2077.499, 132.6063),
  e_oad_ang = e_oad(),
  e_oad_pvd = rbeta(n_size, 3.337215, 1.172535),
  e_oad_fail = e_oad(),
  e_oad_dm_death = rbeta(n_size, 6.865432, 4.97152),
  e_oad_all_cause_death = rbeta(n_size, 10.06827, 5.6634),
  
  e_ins_mi = rbeta(n_size, 9.068628, 5.797975),
  e_ins_stroke = rbeta(n_size, 3.760544, 2.613259),
  e_ins_nephro = 1 - rbeta(n_size, 3.987394, 1.668484),
  e_ins_retino = e_oad(),
  e_ins_neuro = rbeta(n_size, 2077.499, 132.6063),
  e_ins_ang = e_oad(),
  e_ins_pvd = rbeta(n_size, 3.337215, 1.172535),
  e_ins_fail = e_oad(),
  e_ins_dm_death = rbeta(n_size, 6.865432, 4.97152),
  e_ins_all_cause_death = rbeta(n_size, 10.06827, 5.6634),
  
  # Probability of diabetes diagnosis
  p_diag_base = rbeta(n_size, 8.241926, 14.05766), # Assumption
  
  # Effects of strategies
  e_screen_on_diag = rnorm(n_size, 1.5, 0.4), # rr
  e_screen_on_util = rnorm(n_size, 2, 0.5), # rr
  e_comp_on_util = rnorm(n_size, 2, 0.5), # rr
  # p_ad_if_tx = rlnorm(n_size, -0.9877235, 0.3779755), # p
  p_ad_if_tx = rbeta(n_size, 3.866667, 5.8), # p
  
  # HEF policy 
  # p_hef_enrollment = rbeta(n_size, 5.5, 1.833333), 
  # coverage = rbeta(n_size, 17.1, 0.9), 
  # coverage = 1,
  # p_hef_utilization = 1,
  
  # Careseeking probabilities
  p_utilization_hc = rbeta(n_size, 9.755472, 107.7964),
  p_utilization_cpa1 = rbeta(n_size, 10.94302, 574.0505),
  p_utilization_cpa2 = rbeta(n_size, 11.86287, 557.8433),
  p_utilization_cpa3 = rbeta(n_size, 7.304759, 190.2188),

  p_op_hef = rbeta(n_size, 0.9166912, 4.412909),
  p_op_non_hef = rbeta(n_size, 0.4202172, 3.171383),
  p_hosp_base = rbeta(n_size, 8.85, 581.15),

  c_screening = rlnorm(n_size, -0.231791, 0.8088278),
  c_laboratory = rlnorm(n_size, -0.02394995,0.8573677),
  c_oad = rlnorm(n_size, 3.2449579, 0.4055468),
  c_ins = rlnorm(n_size, 4.81576532, 0.1945663),
  c_outpatient_cpa3 = rlnorm(n_size, 3.791499, 0.05900272),
  c_outpatient_cpa2 = rlnorm(n_size, 1.754165, 0.407934), 
  c_outpatient_cpa1 = rlnorm(n_size, 2.298063, 0.2684041),
  c_outpatient_hc = rlnorm(n_size, 1.299967, 0.4962679), 
  c_hospitalization_cpa3 = rlnorm(n_size, 3.703311, 0.1148535),
  c_hospitalization_cpa2 = rlnorm(n_size, 3.378472, 0.1148535), 
  c_hospitalization_cpa1 = rlnorm(n_size, 4.0832387, 0.1148535),
  c_hospitalization_hc = rlnorm(n_size, 1.479141, 0.407934), 
  c_nephro = rlnorm(n_size, 8.733453, 0.2191612), 
  c_retino = rlnorm(n_size, 5.7928000, 0.1148535),
  c_neuro_daily = rlnorm(n_size, 6.862624, 0.0299347),
  c_ang_daily = rlnorm(n_size, 7.855310, 0.0299347),
  c_pvd_daily = rlnorm(n_size, 7.855310, 0.0299347),
  c_mi_daily = rlnorm(n_size, 7.855310, 0.0299347),
  c_stroke_daily = rlnorm(n_size, 6.862624, 0.0299347),
  c_fail_daily = rlnorm(n_size, 7.855310, 0.0299347),
  c_transport_op = rlnorm(n_size, -0.1379559, 0.3303765),
  c_transport_hosp = rlnorm(n_size, 2.452990, 0.06806667),
  monthly_subsistence = 47.96339
) %>% as.matrix()


# 2. Get population -------------------------------------------------------

n_population <- 100000
hef_quantile <- 0.2

cat('Reading in existing population...\n')
ichar <- data.frame(readr::read_csv("output/ichar_200000L_2020-01-25_1855.csv"))[1:n_population,]
person_cycle_x <- load_object_from_rdata('output/person_cycle_x_200000L_2020-01-25_1855.RData')[1:n_population,]

n_cycles <- (dim(person_cycle_x)[2] - 2) / 2
population <- add_hef_to_pop(ichar, person_cycle_x, hef_quantile)


# 3. Run sensitivity analysis ---------------------------------------------

# Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(8)

# Run function for each parameter set
tt0 <- Sys.time()
i_run <- 1

### Function: `run_diabetes_markov`
# cbind(parameters, careseeking, costs, p_treatment, p_utilization)
outcomes_by_run <- apply(modified_parameters_matrix, 1, function(parameters) {
                         
  cat(sprintf('Run %d / %d\n', i_run, n_size))
  
  for (i in 1:length(parameters)) { assign(names(parameters)[i], parameters[[i]]) }
  
  outpatient_costs <- c(
    "HC" = c_outpatient_hc,
    "CPA1" = c_outpatient_cpa1,
    "CPA2" = c_outpatient_cpa2,
    "CPA3" = c_outpatient_cpa3
  )
  
  n_days <- 5
  hospitalization_costs <- c(
    "HC" = c_hospitalization_hc,
    "CPA1" = c_hospitalization_cpa1,
    "CPA2" = c_hospitalization_cpa2,
    "CPA3" = c_hospitalization_cpa3
  ) * n_days
  
  results <- run_model(population, person_cycle_x)
  print(dim(results))
  i_run <<- i_run + 1
  
  return(
    compute_outcomes(population, merge(results, population, by = 'id'), hef_threshold = hef_quantile)
  )
  
})

cat('Run in ', format(Sys.time() - tt0))

psa_results <- bind_rows(outcomes_by_run, .id = 'run')
save_df_to_csv(psa_results, 'psa_results')


