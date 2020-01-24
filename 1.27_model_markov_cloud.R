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
library(readstata13)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')


# I. Hypotheses -----------------------------------------------------------

## 1. Disability weights ---------------------------------------------------

source('0.disability_weights.R')


## 2. Population characteristics -------------------------------------------

source('0.population.R')


## 3. Costs ----------------------------------------------------------------

source('0.costs.R')


## 4. Coverage parameters --------------------------------------------------

## HEF eligibility and enrollment
hef_threshold <- quantile(income_dist, 0.3)[[1]] # Income threshold at which individuals are eligible for HEF (USD); ~ 20th percentile (30th percentile in scenario analysis) of (gamma) income distribution ($28-31 USD) used in the population income distribution parameter defined below
p_hef_enrollment <- 0.75

## Benefits of coverage
coverage <- 1 # 0.8 # If individual is under HEF, 80% of expenditures are covered.

## HEF utilization at point of care
p_hef_utilization = 1 # 0.16


## 5. Transition probabilities ---------------------------------------------

#  Complications; Various sources see Flessa & Zembok 2014 -- total: 0.1189
p_nephro <- 0.01
p_retino <- 0.02124
p_nephro_to_retino <- p_retino * 5.01 # Jeng et al., 2016 -- adjusted HRs for NPDR and PDR were 5.01 (95% confidence interval (CI) = 4.68–5.37) and 9.7 (95% CI = 8.15–11.5), respectively
p_neuro <- 0.04658
p_ang <- 0.00668
p_pvd <- 0.00847 
p_mi <- 0.01735 
p_stroke <- 0.00529
p_fail <- 0.00329
e_nephro_to_cvd <- 1.83 # RR, Gerstein et al., 2001
p_pvd_to_mi <- 0.17 # 0.34 2-year incidence; Rossi et al., 2002
p_pvd_to_ang <- 0.05 # 0.10 2-year incidence; Rossi et al., 2002
e_ang_to_stroke <- 1.14 # HR, Eisen et al., 2016 
e_ang_to_fail <- 1.17 # HR, Eisen et al., 2016
# Note: Other effects non-significant (see Eisen et al., 2016, Grenon et al., 2013)

#  Complication-related death
p_mi_death <- 0.7068 # see Flessa & Zembok 2014; Turner et al. 1998
p_stroke_death <- 0.38136 # see Flessa & Zembok 2014; Turner et al. 1998
p_nephro_death <- 0.311 # Afkarian et al. 2013
p_retino_death <- 0 # Frith & Loprinzi 2018; had to change because it increases beyond 1
p_neuro_death <- 0 
p_ang_death <- 0 
p_pvd_death <- 0.016 * 0.099 # combined probability of (1) amputation and (2) subsequent death within one year; Hoffstad et al. 2015
p_fail_death <- 0.327 # Bell et al. 2003

## Effects of treatments on transition probabilities (relative risks)
#  Diet (no treatment)
e_diet <- 0.90 # Assumption

#  OAD -- focus on metaformin (and glyburide/sulphonylureas) as most likely candidates for prescription in LMIC settings; probability of incidence (or mortality but used as incidence in transition matrix)
e_oad <- 0.5 # Assumption
e_oad_mi <- 0.39 # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS) [61% risk reduction? -> to check]
e_oad_stroke <- e_oad # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS) -- for metaformin
e_oad_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
e_oad_retino <- e_oad # Turner et al. 1998 (UKPDS)
e_oad_neuro <- 0.95 # Juster-Switlyk et al. 2016
e_oad_ang <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_pvd <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_fail <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)

#  Insulin -- same as above
e_ins <- 0.5 # Assumption
e_ins_mi <- 0.39 # Turner et al. 1998 (UKPDS) - did not differ from metaformin group; Herman et al. 2017; ranges to increased risk for CV events - see Table 1 (uses last registry entry for MACE)
e_ins_stroke <- e_ins # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS)
e_ins_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
e_ins_retino <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_neuro <- 0.95 # Turner et al. 1998 (UKPDS)
e_ins_ang <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_pvd <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_fail <- e_ins # Turner et al. 1998 (UKPDS)


## 6. States and strategies ------------------------------------------------

## States
stateNames <- c("healthy", "diabetes", "notx", "oad", "ins", "nephro", "retino", "neuro", "ang", "pvd", "mi", "stroke", "fail", "dm_death", "other_death")

## Distribution of Cambodian population
initialStates_KHR <- function(age, sex, p_diet, p_oad) {
  prevalence <- ifelse(sex == "Male", 0.057, 0.061) # WHO Fact Sheet, Cambodia
  prev_nephro <- ifelse(age > 52, 0.012, 0) # Thomas et al., 2014
  prev_retino <- ifelse(age > 52, 0.012, 0) # Sogbesan and Yutho, 2000
  prev_neuro <- ifelse(age > 52, 0.05, 0) # Lim et al., 2002
  p_oad <- 0.224
  p_ins <- 0.017
  p_diet <- (1 - p_oad - p_ins)
  init_states_KHR <- c((1 - prevalence * (1 + p_diet + p_oad + p_ins + prev_nephro + prev_retino + prev_neuro)), 
                       prevalence, prevalence * p_diet, prevalence * p_oad, prevalence * p_ins, 
                       prevalence * prev_nephro, prevalence * prev_retino, prevalence * prev_neuro, 
                       0, 0, 0, 0, 0, 0, 0) # CVD-related data (among diabetics) not available
  return(init_states_KHR)
}

## Strategies to simulate
strategyNames <- c("base", "screen_only", "tx_only", "comp_only", "screen_tx", "tx_comp", "screen_tx_comp")

# Effects of strategies
e_screen_on_diag <- 1.5
e_screen_on_util <- 2 # RR
e_comp_on_util <- 2 # RR
e_tx <- 0.40


# II. Simulate populations ----------------------------------------------------

n_population <- as.integer(5e5) # total: 16e6
# n_population <- 1000

# Time horizon
n_cycles = 45

# Load required functions 
source("1.22_model_cost.R")
source("1.20_model_dalys.R")
source("1.FRP.R")

create_new_population <- FALSE
modify_hef <- TRUE

if (create_new_population) {
  
  cat('Simulating population of', n_population, 'people')
  ichar <- tibble(
    id = 1:n_population,
    age_start = sample(age_dist, n_population, replace = TRUE),
    sex = ifelse(rbinom(n = n_population, size = 1, p = p_female) == 0, "Male", "Female"), # Relabel for `look_up`
    income = rgamma(n = n_population, shape = 0.5, scale = avg_income),
    disposable_income = income - subsistence,
    hef = ifelse(income < hef_threshold, rbinom(n = n_population, size = 1, p = p_hef_enrollment), 0),
    hef_utilization = ifelse(hef, rbinom(n = n_population, size = 1, p = p_hef_utilization), NA),
    init_state = mapply(function(age, sex) sample(stateNames, size = 1, prob = initialStates_KHR(age, sex)), age_start, sex)
  ) %>% as.data.frame
  
  ## Matrix of state probabilities (pointers) for each individual and cycle
  person_cycle_x <- matrix(runif(n_population * n_cycles), n_population, n_cycles)

  save_df_to_csv(ichar, "ichar")
  save_object_to_rdata(person_cycle_x)
  
  cat(' Done\n')
  
  } else {
    cat('Reading in existing population...\n')
    dir("output", sprintf("(ichar|person_cycle_x)_%dL.*", n_population))
    ichar <- data.frame(readr::read_csv("output/ichar_500000L_2020-01-23_1757.csv"))
    person_cycle_x <- load_object_from_rdata('output/person_cycle_x_500000L_2020-01-23_1757.RData')
    if (modify_hef) {
      cat('Mutating HEF status.\n')
      ichar <- ichar %>% mutate(
        hef = ifelse(income < hef_threshold, rbinom(n = 1, size = 1, p = p_hef_enrollment), 0),
        hef_utilization = ifelse(hef, rbinom(n = 1, size = 1, p = p_hef_utilization), NA)
      )
      save_df_to_csv(ichar, "ichar")
    }
}



# III. Run Model ----------------------------------------------------------

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
      
    ## Instantiate the data frame to store the results
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
    
    ## Initialize the transition matrix
    tmat <- matrix(NA, nrow = 15, ncol = 15)
    colnames(tmat) <- stateNames
    rownames(tmat) <- stateNames
    
    ## Select population indices
    population_indices <- ((i_sample-1)*sample_size+1) : min(((i_sample)*sample_size), n_population)
    for (i_person in population_indices) {
      
      if (i_person %% 1000 == 0) {
        cat('.')
      }

      id <- ichar[i_person, 'id']
      age <- as.numeric(ichar[i_person, 'age_start'])
      sex <- ichar[i_person, 'sex']
      income <- as.numeric(ichar[i_person, 'income'])
      disposable_income <- as.numeric(ichar[i_person, 'disposable_income'])
      hef <- ichar[i_person, 'hef']
      hef_utilization <- ichar[i_person, 'hef_utilization']
      init_state <- ichar[i_person, 'init_state']
      
      ## Set relevant probabilities
      is_strategy_screen <- strategy == "screen_only" | strategy == "screen_tx" | strategy == "screen_tx_comp"
      is_strategy_tx <- !(strategy == "base" | strategy == "screen_only" | strategy == "comp_only")
      is_strategy_comp <- strategy == "comp_only" | strategy == "tx_comp" | strategy == "screen_tx_comp"

      
      #  Diagnosis; Source: Flessa & Zembok 2014
      p_diag <- 0.3696 * ifelse(is_strategy_screen & hef, e_screen_on_diag, 1)
      
      #  Treatment
      p_oad <- 0.224 # Taniguchi et al., 2016; differs significantly from King et al., 2005
      p_ins <- 0.017 # Van Olmen et al. 2015
      p_diet <- (1 - p_oad - p_ins)
      p_oad_to_ins <- 0.04 # Ringborg et al. 2010
      p_oad_ad <- ifelse(is_strategy_tx & hef, e_tx, 0.125)
      p_ins_ad <- ifelse(is_strategy_tx & hef, e_tx, 0.125)
      
      #  Utilization
      e_strategy_util <- ifelse(is_strategy_screen & hef, e_screen_on_util, 1) * ifelse(is_strategy_comp & hef, e_comp_on_util, 1)
      p_op <- ifelse(hef, 0.172, 0.117) * e_strategy_util
      p_hosp <- 0.015 * e_strategy_util
      
      # Cycle through time horizon
      for (cycle in 1:n_cycles) {
        
        # Assign parameters based on individual characteristics
        p_other_death <- all_cause_mortality[ all_cause_mortality$age == age & all_cause_mortality$sex == sex, "mr"]
        p_dm_death <- dm_mortality[ dm_mortality$age == age & dm_mortality$sex_name == sex, "val"]
        p_dm <- dm_incidence[dm_incidence$age == age & dm_incidence$sex_name == sex, "val"]
        
        ## Define all values of tmat
        diag(tmat) <- 0
        
        tmat["healthy", "diabetes"] <- p_dm * (1 - p_diag)
        tmat["healthy", "notx"] <- p_dm * p_diag * p_diet
        tmat["healthy", "oad"] <- p_dm * p_diag * p_oad
        tmat["healthy", "ins"] <- p_dm * p_diag * p_ins
        tmat["healthy", "nephro"] <- p_dm * p_nephro
        tmat["healthy", "retino"] <- p_dm * p_retino
        tmat["healthy", "neuro"] <- p_dm * p_neuro
        tmat["healthy", "ang"] <- p_dm * p_ang
        tmat["healthy", "pvd"] <- p_dm * p_pvd
        tmat["healthy", "mi"] <- p_dm * p_mi
        tmat["healthy", "stroke"] <- p_dm * p_stroke
        tmat["healthy", "fail"] <- p_dm * p_fail
        tmat["healthy", "dm_death"] <- p_dm * p_dm_death
        tmat["healthy", "other_death"] <- p_other_death
        tmat["healthy", "healthy"] <- 1 - sum(tmat["healthy",], na.rm = T)
        
        tmat["diabetes", "healthy"] <- 0
        tmat["diabetes", "notx"] <- p_diag * p_diet
        tmat["diabetes", "oad"] <- p_diag * p_oad
        tmat["diabetes", "ins"] <- p_diag * p_ins
        tmat["diabetes", "nephro"] <- p_nephro 
        tmat["diabetes", "retino"] <- p_retino 
        tmat["diabetes", "neuro"] <- p_neuro 
        tmat["diabetes", "ang"] <- p_ang 
        tmat["diabetes", "pvd"] <- p_pvd 
        tmat["diabetes", "mi"] <- p_mi 
        tmat["diabetes", "stroke"] <- p_stroke 
        tmat["diabetes", "fail"] <- p_fail 
        tmat["diabetes", "dm_death"] <- p_dm_death
        tmat["diabetes", "other_death"] <- p_other_death
        tmat["diabetes", "diabetes"] <- 1 - sum(tmat["diabetes",], na.rm = T)
        
        tmat["notx", "healthy"] <- 0
        tmat["notx", "diabetes"] <- 0
        tmat["notx", "oad"] <- p_oad 
        tmat["notx", "ins"] <- p_ins
        tmat["notx", "nephro"] <- p_nephro * e_diet
        tmat["notx", "retino"] <- p_retino * e_diet
        tmat["notx", "neuro"] <- p_neuro * e_diet
        tmat["notx", "ang"] <- p_ang * e_diet
        tmat["notx", "pvd"] <- p_pvd * e_diet
        tmat["notx", "mi"] <- p_mi * e_diet
        tmat["notx", "stroke"] <- p_stroke * e_diet
        tmat["notx", "fail"] <- p_fail * e_diet
        tmat["notx", "dm_death"] <- p_dm_death
        tmat["notx", "other_death"] <- p_other_death
        tmat["notx", "notx"] <- 1 - sum(tmat["notx",], na.rm = T)
        
        tmat["oad", "healthy"] <- 0
        tmat["oad", "diabetes"] <- 0
        tmat["oad", "notx"] <- 0
        tmat["oad", "ins"] <- p_oad_to_ins
        tmat["oad", "nephro"] <- p_nephro * (1 - p_oad_ad * e_oad_nephro)
        tmat["oad", "retino"] <- p_retino * (1 - p_oad_ad * e_oad_retino)
        tmat["oad", "neuro"] <- p_neuro * (1 - p_oad_ad * e_oad_neuro)
        tmat["oad", "ang"] <- p_ang * (1 - p_oad_ad * e_oad_ang)
        tmat["oad", "pvd"] <- p_pvd * (1 - p_oad_ad * e_oad_pvd)
        tmat["oad", "mi"] <- p_mi * (1 - p_oad_ad * e_oad_mi)
        tmat["oad", "stroke"] <- p_stroke * (1 - p_oad_ad * e_oad_stroke)
        tmat["oad", "fail"] <- p_fail * (1 - p_oad_ad * e_oad_fail)
        tmat["oad", "dm_death"] <- p_dm_death
        tmat["oad", "other_death"] <- p_other_death
        tmat["oad", "oad"] <- 1 - sum(tmat["oad",], na.rm = T)
        
        tmat["ins", "healthy"] <- 0
        tmat["ins", "diabetes"] <- 0
        tmat["ins", "notx"] <- 0
        tmat["ins", "oad"] <- 0
        tmat["ins", "nephro"] <- p_nephro * (1 - p_ins_ad * e_ins_nephro)
        tmat["ins", "retino"] <- p_retino * (1 - p_ins_ad * e_ins_retino)
        tmat["ins", "neuro"] <- p_neuro * (1 - p_ins_ad * e_ins_neuro)
        tmat["ins", "ang"] <- p_ang * (1 - p_ins_ad * e_ins_ang)
        tmat["ins", "pvd"] <- p_pvd * (1 - p_ins_ad * e_ins_pvd)
        tmat["ins", "mi"] <- p_mi * (1 - p_ins_ad * e_ins_mi)
        tmat["ins", "stroke"] <- p_stroke * (1 - p_ins_ad * e_ins_stroke)
        tmat["ins", "fail"] <- p_fail * (1 - p_ins_ad * e_ins_fail)
        tmat["ins", "dm_death"] <- p_dm_death
        tmat["ins", "other_death"] <- p_other_death
        tmat["ins", "ins"] <- 1 - sum(tmat["ins",], na.rm = T)
        
        tmat["nephro", "healthy"] <- 0
        tmat["nephro", "diabetes"] <- 0
        tmat["nephro", "notx"] <- 0
        tmat["nephro", "oad"] <- 0
        tmat["nephro", "ins"] <- 0
        tmat["nephro", "retino"] <- p_nephro_to_retino
        tmat["nephro", "neuro"] <- p_neuro
        tmat["nephro", "ang"] <- p_ang * e_nephro_to_cvd
        tmat["nephro", "pvd"] <- p_pvd * e_nephro_to_cvd
        tmat["nephro", "mi"] <- p_mi * e_nephro_to_cvd
        tmat["nephro", "stroke"] <- p_stroke
        tmat["nephro", "fail"] <- p_fail * e_nephro_to_cvd
        tmat["nephro", "dm_death"] <- p_nephro_death
        tmat["nephro", "other_death"] <- p_other_death
        tmat["nephro", "nephro"] <- 1 - sum(tmat["nephro",], na.rm = T)
        
        tmat["retino", "healthy"] <- 0
        tmat["retino", "diabetes"] <- 0
        tmat["retino", "notx"] <- 0
        tmat["retino", "oad"] <- 0
        tmat["retino", "ins"] <- 0
        tmat["retino", "nephro"] <- p_nephro
        tmat["retino", "neuro"] <- p_neuro
        tmat["retino", "ang"] <- p_ang
        tmat["retino", "pvd"] <- p_pvd
        tmat["retino", "mi"] <- p_mi
        tmat["retino", "stroke"] <- p_stroke
        tmat["retino", "fail"] <- p_fail
        tmat["retino", "dm_death"] <- p_retino_death
        tmat["retino", "other_death"] <- p_other_death
        tmat["retino", "retino"] <- 1 - sum(tmat["retino",], na.rm = T)
        
        tmat["neuro", "healthy"] <- 0
        tmat["neuro", "diabetes"] <- 0
        tmat["neuro", "notx"] <- 0
        tmat["neuro", "oad"] <- 0
        tmat["neuro", "ins"] <- 0
        tmat["neuro", "nephro"] <- p_nephro
        tmat["neuro", "retino"] <- p_retino
        tmat["neuro", "ang"] <- p_ang
        tmat["neuro", "pvd"] <- p_pvd
        tmat["neuro", "mi"] <- p_mi
        tmat["neuro", "stroke"] <- p_stroke
        tmat["neuro", "fail"] <- p_fail
        tmat["neuro", "dm_death"] <- p_neuro_death
        tmat["neuro", "other_death"] <- p_other_death
        tmat["neuro", "neuro"] <- 1 - sum(tmat["neuro",], na.rm = T)
        
        tmat["ang", "healthy"] <- 0
        tmat["ang", "diabetes"] <- 0
        tmat["ang", "notx"] <- 0
        tmat["ang", "oad"] <- 0
        tmat["ang", "ins"] <- 0
        tmat["ang", "nephro"] <- p_nephro
        tmat["ang", "retino"] <- p_retino
        tmat["ang", "neuro"] <- p_neuro
        tmat["ang", "pvd"] <- p_pvd
        tmat["ang", "mi"] <- p_mi
        tmat["ang", "stroke"] <- p_stroke * e_ang_to_stroke
        tmat["ang", "fail"] <- p_fail * e_ang_to_fail
        tmat["ang", "dm_death"] <- p_ang_death
        tmat["ang", "other_death"] <- p_other_death
        tmat["ang", "ang"] <- 1 - sum(tmat["ang",], na.rm = T)
        
        tmat["pvd", "healthy"] <- 0
        tmat["pvd", "diabetes"] <- 0
        tmat["pvd", "notx"] <- 0
        tmat["pvd", "oad"] <- 0
        tmat["pvd", "ins"] <- 0
        tmat["pvd", "nephro"] <- p_nephro
        tmat["pvd", "retino"] <- p_retino
        tmat["pvd", "neuro"] <- p_neuro
        tmat["pvd", "ang"] <- p_pvd_to_ang
        tmat["pvd", "mi"] <- p_pvd_to_mi
        tmat["pvd", "stroke"] <- p_stroke
        tmat["pvd", "fail"] <- p_fail
        tmat["pvd", "dm_death"] <- p_pvd_death
        tmat["pvd", "other_death"] <- p_other_death
        tmat["pvd", "pvd"] <- 1 - sum(tmat["pvd",], na.rm = T)
        
        tmat["mi", "healthy"] <- 0
        tmat["mi", "diabetes"] <- 0
        tmat["mi", "nephro"] <- p_nephro
        tmat["mi", "retino"] <- p_retino
        tmat["mi", "neuro"] <- p_neuro
        tmat["mi", "ang"] <- p_ang
        tmat["mi", "pvd"] <- p_pvd
        tmat["mi", "mi"] <- p_mi
        tmat["mi", "stroke"] <- p_stroke
        tmat["mi", "fail"] <- p_fail
        tmat["mi", "dm_death"] <- p_mi_death
        tmat["mi", "other_death"] <- p_other_death
        p_mi_diagnosed <- (1 - sum(tmat["mi", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
        tmat["mi", "notx"] <- p_mi_diagnosed * p_diet
        tmat["mi", "oad"] <- p_mi_diagnosed * p_oad
        tmat["mi", "ins"] <- p_mi_diagnosed * p_ins
        
        tmat["stroke", "healthy"] <- 0
        tmat["stroke", "diabetes"] <- 0
        tmat["stroke", "nephro"] <- p_nephro
        tmat["stroke", "retino"] <- p_retino
        tmat["stroke", "neuro"] <- p_neuro
        tmat["stroke", "ang"] <- p_ang
        tmat["stroke", "pvd"] <- p_pvd
        tmat["stroke", "mi"] <- p_mi
        tmat["stroke", "stroke"] <- p_stroke
        tmat["stroke", "fail"] <- p_fail
        tmat["stroke", "dm_death"] <- p_mi_death
        tmat["stroke", "other_death"] <- p_other_death
        p_stroke_diagnosed <- (1 - sum(tmat["stroke", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
        tmat["stroke", "notx"] <- p_stroke_diagnosed * p_diet
        tmat["stroke", "oad"] <- p_stroke_diagnosed * p_oad
        tmat["stroke", "ins"] <- p_stroke_diagnosed * p_ins
        
        tmat["fail", "healthy"] <- 0
        tmat["fail", "diabetes"] <- 0
        tmat["fail", "nephro"] <- p_nephro
        tmat["fail", "retino"] <- p_retino
        tmat["fail", "neuro"] <- p_neuro
        tmat["fail", "ang"] <- p_ang
        tmat["fail", "pvd"] <- p_pvd
        tmat["fail", "mi"] <- p_mi
        tmat["fail", "stroke"] <- p_stroke
        tmat["fail", "fail"] <- p_fail
        tmat["fail", "dm_death"] <- p_mi_death
        tmat["fail", "other_death"] <- p_other_death
        p_fail_diagnosed <- (1 - sum(tmat["fail", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
        tmat["fail", "notx"] <- p_fail_diagnosed * p_diet
        tmat["fail", "oad"] <- p_fail_diagnosed * p_oad
        tmat["fail", "ins"] <- p_fail_diagnosed * p_ins
        
        tmat["dm_death", "healthy"] <- 0
        tmat["dm_death", "diabetes"] <- 0
        tmat["dm_death", "notx"] <- 0
        tmat["dm_death", "oad"] <- 0
        tmat["dm_death", "ins"] <- 0
        tmat["dm_death", "nephro"] <- 0
        tmat["dm_death", "retino"] <- 0
        tmat["dm_death", "neuro"] <- 0
        tmat["dm_death", "ang"] <- 0
        tmat["dm_death", "pvd"] <- 0
        tmat["dm_death", "mi"] <- 0
        tmat["dm_death", "stroke"] <- 0
        tmat["dm_death", "fail"] <- 0
        tmat["dm_death", "dm_death"] <- 1
        tmat["dm_death", "other_death"] <- 0
        
        tmat["other_death", "healthy"] <- 0
        tmat["other_death", "diabetes"] <- 0
        tmat["other_death", "notx"] <- 0
        tmat["other_death", "oad"] <- 0
        tmat["other_death", "ins"] <- 0
        tmat["other_death", "nephro"] <- 0
        tmat["other_death", "retino"] <- 0
        tmat["other_death", "neuro"] <- 0
        tmat["other_death", "ang"] <- 0
        tmat["other_death", "pvd"] <- 0
        tmat["other_death", "mi"] <- 0
        tmat["other_death", "stroke"] <- 0
        tmat["other_death", "fail"] <- 0
        tmat["other_death", "dm_death"] <- 0
        tmat["other_death", "other_death"] <- 1
        
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
  
  # Combine results into single dataframes (results, ichar)
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

