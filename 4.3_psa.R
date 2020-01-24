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
source("1.22_model_cost.R")
source("1.20_model_dalys.R")
source("1.FRP.R")
source("1.outcomes.R")

## Load static parameters 
source('0.population.R') # Population characteristics

## Strategies to simulate
strategyNames <- c("base", "screen_only", "tx_only", "comp_only", "screen_tx", "tx_comp", "screen_tx_comp")

# I. Model Function ----------------------------------------------------------

### Function: `run_diabetes_markov`

run_diabetes_markov <- function(create_new_population = FALSE, modify_hef = FALSE,
                                n_population = 20000, n_cycles = 45, 
                                costs, parameters, careseeking = NULL, p_treatment = NULL, p_utilization = NULL, 
                                person_cycle_x = NULL) {
  
  ## Define parameters from inputs
  #  Complications
  p_nephro <- parameters[['p_nephro']]
  p_retino <- parameters[['p_retino']]
  p_nephro_to_retino <- parameters[['p_retino']] * parameters[['rr_nephro_to_retino']] # Jeng et al., 2016 -- adjusted HRs for NPDR and PDR were 5.01 (95% confidence interval (CI) = 4.68–5.37) and 9.7 (95% CI = 8.15–11.5), respectively
  p_neuro <- parameters[['p_neuro']]
  p_ang <- parameters[['p_ang']]
  p_pvd <- parameters[['p_pvd']] 
  p_mi <- parameters[['p_mi']] 
  p_stroke <- parameters[['p_stroke']]
  p_fail <- parameters[['p_fail']]
  e_nephro_to_cvd <- parameters[['rr_nephro_to_cvd']] # RR, Gerstein et al., 2001
  p_pvd_to_mi <- parameters[['p_pvd_to_mi']] # 0.34 2-year incidence; Rossi et al., 2002
  p_pvd_to_ang <- parameters[['p_pvd_to_ang']] # 0.10 2-year incidence; Rossi et al., 2002
  e_ang_to_stroke <- parameters[['rr_ang_to_stroke']] # HR, Eisen et al., 2016 
  e_ang_to_fail <- parameters[['rr_ang_to_fail']] # HR, Eisen et al., 2016
  
  #  Complication-related death
  p_mi_death <- parameters[['p_mi_death']] # see Flessa & Zembok 2014; Turner et al. 1998
  p_stroke_death <- parameters[['p_stroke_death']] # see Flessa & Zembok 2014; Turner et al. 1998
  p_nephro_death <- parameters[['p_nephro_death']] # Afkarian et al. 2013
  p_retino_death <- parameters[['p_retino_death']] # Frith & Loprinzi 2018; had to change because it increases beyond 1
  p_neuro_death <- parameters[['p_neuro_death']] 
  p_ang_death <- parameters[['p_ang_death']]
  p_pvd_death <- parameters[['p_pvd_death']] # combined probability of (1) amputation and (2) subsequent death within one year; Hoffstad et al. 2015
  p_fail_death <- parameters[['p_fail_death']] # Bell et al. 2003
  
  # Effects of treatments on transition probabilities (relative risks)
  #  Diet (no treatment)
  e_diet <- parameters[['rr_diet']] # Assumption
  
  #  OAD -- focus on metaformin (and glyburide/sulphonylureas) as most likely candidates for prescription in LMIC settings; probability of incidence (or mortality but used as incidence in transition matrix)
  e_oad <- parameters[['rr_oad']] # Assumption
  e_oad_mi <- parameters[['rr_oad_mi']] # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS) [61% risk reduction? -> to check]
  e_oad_stroke <- parameters[['rr_oad_stroke']] # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS) -- for metaformin
  e_oad_nephro <- parameters[['rr_oad_nephro']] # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
  e_oad_retino <- parameters[['rr_oad_retino']] # Turner et al. 1998 (UKPDS)
  e_oad_neuro <- parameters[['rr_oad_neuro']] # Juster-Switlyk et al. 2016
  e_oad_ang <- parameters[['rr_oad_ang']] # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
  e_oad_pvd <- parameters[['rr_oad_pvd']] # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
  e_oad_fail <- parameters[['rr_oad_fail']] # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
  
  #  Insulin -- same as above
  e_ins <- parameters[['rr_ins']] # Assumption
  e_ins_mi <- parameters[['rr_ins_mi']] # Turner et al. 1998 (UKPDS) - did not differ from metaformin group; Herman et al. 2017; ranges to increased risk for CV events - see Table 1 (uses last registry entry for MACE)
  e_ins_stroke <- parameters[['rr_ins_stroke']] # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS)
  e_ins_nephro <- parameters[['rr_ins_nephro']] # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
  e_ins_retino <- parameters[['rr_ins_retino']] # Turner et al. 1998 (UKPDS)
  e_ins_neuro <- parameters[['rr_ins_neuro']] # Turner et al. 1998 (UKPDS)
  e_ins_ang <- parameters[['rr_ins_ang']] # Turner et al. 1998 (UKPDS)
  e_ins_pvd <- parameters[['rr_ins_pvd']] # Turner et al. 1998 (UKPDS)
  e_ins_fail <- parameters[['rr_ins_fail']] # Turner et al. 1998 (UKPDS)
  
  # Careseeking probabilities (at which facility, scaled based on overall probability of care utilization) (DHS 2014, Care seeking for first treatment, total population)
  if (is.null(careseeking)) {
    p_utilization_hc <- 0.114
    p_utilization_cpa1 <- 0.025
    p_utilization_cpa2 <- 0.030
    p_utilization_cpa3 <- 0.042
    p_provider <- c(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3) / sum(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3)
  } else {
    p_utilization_hc <- careseeking[['p_utilization_hc']]
    p_utilization_cpa1 <- careseeking[['p_utilization_cpa1']]
    p_utilization_cpa2 <- careseeking[['p_utilization_cpa2']]
    p_utilization_cpa3 <- careseeking[['p_utilization_cpa3']]
    p_provider <- c(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3) / sum(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3)
  }
  
  # Effects of strategies on transition probabilities
  e_screen_on_diag <- parameters[['e_screen_on_diag']] # rr
  e_screen_on_util <- parameters[['e_screen_on_util']] # rr
  e_comp_on_util <- parameters[['e_comp_on_util']] # rr
  e_tx <- parameters[['e_tx']] # probability
  
  # Probabilities of being prescribed treatment
  if (is.null(p_treatment)) {
    p_oad <- 0.224 # Taniguchi et al., 2016; differs significantly from King et al., 2005
    p_ins <- 0.017 # Van Olmen et al. 2015
  } else {
    p_oad <- p_treatment[['p_oad']]
    p_ins <- p_treatment[['p_ins']]
  }
  
  # Probabilities of utilization
  if (is.null(p_utilization)) {
    p_op_hef <- 0.172
    p_op_nohef <- 0.117
    p_hospitalization <- 0.015
  } else {
    p_op_hef <- p_utilization[['p_op_hef']]
    p_op_nohef <- p_utilization[['p_op_nohef']]
    p_hospitalization <- p_utilization[['p_hospitalization']]
  }
  
  ## Costs
  # Direct medical costs; Source: Flessa & Zembok, 2014; adjusted to 2019 USD
  c_screening <- costs[['c_screening']] # i.e. screening; unit cost FPG-test
  c_laboratory <- costs[['c_laboratory']] # diagnosis cost
  c_oad <- costs[['c_oad']]
  c_ins <- costs[['c_ins']]
  c_outpatient_cpa3 <- costs[['c_outpatient_cpa3']] # Maximum prices -- Costing of Health Care Services in Three Provinces of Cambodia (2018)
  c_outpatient_cpa2 <- costs[['c_outpatient_cpa2']] # Secondary hospital
  c_outpatient_cpa1 <- costs[['c_outpatient_cpa1']] # Primary hospital
  c_outpatient_hc <- costs[['c_outpatient_hc']] 
  c_hospitalization_cpa3 <- costs[['c_hospitalization_cpa3']]
  c_hospitalization_cpa2 <- costs[['c_hospitalization_cpa2']]
  c_hospitalization_cpa1 <- costs[['c_hospitalization_cpa1']]
  c_hospitalization_hc <- costs[['c_hospitalization_hc']]
  
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
  
  # Direct medical costs associated with each complication
  c_nephro <- costs[['c_nephro']] # Hemodialysis, systematic review of LMICs; Int$ 3,424 to Int$ 42,785 overall, but used Sri Lanka numbers: Int$ 5,869–8,804; Mushi et al., 2015 -- adjusted for inflation
  c_retino <- costs[['c_retino']] # Cost of intravitreal injection in Indonesia, Sasongko et al., 2019
  c_neuro <- costs[['c_neuro_daily']] * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_ang <- costs[['c_ang_daily']] * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_pvd <- costs[['c_pvd_daily']] * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_mi <- costs[['c_mi_daily']] * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_stroke <- costs[['c_stroke_daily']] * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_fail <- costs[['c_fail_daily']] * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  
  ## Individual expenditures
  # Indirect (non-medical) costs
  oop_transport_op <- costs[['c_transport_op']]
  oop_transport_hosp <- costs[['c_transport_hosp']]
  
  # Discount rate
  dr <- 0.03
  
  # Subsistence expenditure
  subsistence <- costs[['monthly_subsistence']] * 12 # 196K KHR (CSES 2017) / 4086.45 KHR per USD (2 Oct 2019); average monthly value per capital, 2017
  
  ## Disability weights
  # Diagnosed, uncomplicated DM 
  dw_uncomplicated <- parameters[['uncomplicated']] # 0.049 (0.031 - 0.072)
  
  # Nephropathy
  dw_nephropathy_stage5 <- parameters[['dw_nephropathy']] # 0.569 (0.389 - 0.727)
  
  # Retinopathy
  dw_retinopathy <- parameters[['dw_retinopathy']] # 0.184 (0.125 - 0.258)
  
  # Distal symmetric neuropathy
  dw_neuropathy <- parameters[['dw_neuropathy']] # 0.133 (0.089 - 0.187)
  
  # Angina pectoris
  dw_angina_moderate <- parameters[['dw_angina']] # 0.080 (0.052 - 0.113)
  
  # Peripheral vascular disease
  dw_pvd <- parameters[['dw_pvd']] # 0.014 (0.007 - 0.025)
  
  # Myocardial infarction
  dw_mi <- parameters[['dw_mi']] # 0.432 (0.288 - 0.579)
  
  # Stroke
  dw_stroke_level5 <- parameters[['dw_stroke']] # 0.588 (0.411 - 0.744)
  
  # Heart failure
  dw_failure_severe <- parameters[['dw_failure']] # 0.179 (0.122 - 0.251)
  
  ## HEF
  p_hef_enrollment <- parameters[['p_enrollment']]
  p_hef_utilization <- parameters[['p_utilization']]
  hef_threshold <- parameters[['threshold']]
  coverage <- parameters[['coverage']]
  
  ## Designating matrix of state probabilities (pointers or 'fate') for each individual and cycle
  if (is.null(person_cycle_x)) {
    person_cycle_x <- load_object_from_rdata('output/person_cycle_x_500000L_2020-01-23_1757.RData')
  } else {
    person_cycle_x <- person_cycle_x
  }

  ## Define population
  if (create_new_population) {
    
    cat('Simulating population of', n_population, 'people')
    ichar <- tibble(
      id = 1:n_population,
      age_start = sample(age_dist, n_population, replace = TRUE),
      sex = ifelse(rbinom(n = n_population, size = 1, p = p_female) == 0, "Male", "Female"), 
      income = rgamma(n = n_population, shape = 0.5, scale = avg_income),
      disposable_income = income - subsistence,
      hef = ifelse(income < hef_threshold, rbinom(n = n_population, size = 1, p = p_hef_enrollment), 0),
      hef_utilization = ifelse(hef, rbinom(n = n_population, size = 1, p = p_hef_utilization), NA),
      init_state = mapply(function(age, sex) sample(stateNames, size = 1, prob = initialStates_KHR(age, sex)), age_start, sex)
    ) %>% as.data.frame
    
    save_df_to_csv(ichar, "ichar")
    
    cat(' Done\n')
    
  } else {
    cat('Reading in existing population...\n')
    dir("output", sprintf("(ichar|person_cycle_x)_%dL.*", n_population))
    ichar <- data.frame(readr::read_csv("output/ichar_500000L_2020-01-23_1757.csv"))[1:n_population,]
    
    if (modify_hef) {
      cat('Mutating HEF status.\n')
      ichar <- ichar %>% mutate(
        hef = ifelse(income < hef_threshold, rbinom(n = 1, size = 1, p = p_hef_enrollment), 0),
        hef_utilization = ifelse(hef, rbinom(n = 1, size = 1, p = p_hef_utilization), NA)
      )
      
      save_df_to_csv(ichar, "ichar")
      
    }
  }
  
  # Start parallel workers
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
        p_diag <- parameters[['p_diagnosis']] * ifelse(is_strategy_screen & hef, e_screen_on_diag, 1)
        
        #  Treatment
        p_oad <- p_oad # Taniguchi et al., 2016; differs significantly from King et al., 2005
        p_ins <- p_ins # Van Olmen et al. 2015
        p_diet <- (1 - p_oad - p_ins)
        p_oad_to_ins <- 0.04 # Ringborg et al. 2010
        p_oad_ad <- ifelse(is_strategy_tx & hef, e_tx, 0.125)
        p_ins_ad <- ifelse(is_strategy_tx & hef, e_tx, 0.125)
        
        #  Utilization
        e_strategy_util <- ifelse(is_strategy_screen & hef, e_screen_on_util, 1) * ifelse(is_strategy_comp & hef, e_comp_on_util, 1)
        p_op <- ifelse(hef, p_op_hef, p_op_nohef) * e_strategy_util
        p_hosp <- p_hospitalization * e_strategy_util
        
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
    all_results[[strategy]] <- do.call(rbind, results_per_sample)
  }
  
  # Combine results
  return(
    results = do.call(rbind, all_results) %>% as.data.frame,
    population = ichar
  )
}

run_and_compute_markov <- function(create_new_population = FALSE, modify_hef = FALSE,
                                   n_population = 20000, n_cycles = 45, 
                                   costs, parameters, careseeking = NULL, p_treatment = NULL, p_utilization = NULL, 
                                   person_cycle_x = NULL) {
  results_and_pop <- run_diabetes_markov(...)

  # Compute outcomes for each strategy (by income quintile and sex)
  outcomes <- compute_outcomes(results_and_pop$population, results_and_pop$results, hef_threshold = 0.2)

  # Return list of results
  return(list("outcomes" = outcomes, "outcomes_by_sex" = outcomes_by_sex))
}


# II. Parameter Distributions ---------------------------------------------

# Define parameter distributions in dataframe
n_size = 100
parameters <- data.frame(
  
  # Probabilities of onset of complications
  p_nephro = rep(0.01, n_size),
  p_retino = rep(0.02124, n_size),
  p_neuro = rep(0.04658, n_size),
  p_ang = rep(0.00668, n_size),
  p_pvd = rep(0.00847, n_size),
  p_mi = rep(0.01735, n_size), 
  p_stroke = rep(0.00529, n_size), 
  p_fail = rep(0.00329, n_size),
  e_nephro_to_cvd = rnorm(n_size, 1.83, 0.11) ,
  p_pvd_to_mi = rep(0.17, n_size),
  p_pvd_to_ang = rep(0.05, n_size),
  e_ang_to_stroke = rep(1.14, n_size), 
  e_ang_to_fail = rep(1.17, n_size), 
  
  # Probabilities of mortality
  p_mi_death = rep(0.7068, n_size),
  p_stroke_death = rep(0.38136, n_size), 
  p_nephro_death = rep(0.311, n_size), 
  p_retino_death = rep(0, n_size), 
  p_neuro_death = rep(0, n_size), 
  p_ang_death = rep(0, n_size), 
  p_pvd_death = rep(0.016 * 0.099, n_size),
  p_fail_death = rep(0.327, n_size),
  
  # Effects of therapies
  e_diet = rnorm(n_size, 0.90, 0.05), 
  e_oad = rnorm(n_size, 0.50, 0.10), 
  e_oad_mi = rep(0.39, n_size), 
  e_oad_stroke = rbeta(n_size, 26.00386, 8.667955), # Herman et al., 2017
  e_oad_nephro = rbeta(n_size, 942.1698, 413.4702),
  e_oad_neuro = rep(0.95, n_size), 
  e_oad_ang = rbeta(n_size, 26.00386, 8.667955), 
  e_oad_pvd = rbeta(n_size, 26.00386, 8.667955), 
  e_oad_fail = rbeta(n_size, 26.00386, 8.667955), 
  e_ins = rnorm(n_size, 0.50, 0.10), 
  e_ins_mi = 0.39, 
  e_ins_stroke = rbeta(n_size, 26.00386, 8.667955), 
  e_ins_nephro = rbeta(n_size, 235.0212, 103.1388),
  e_ins_neuro = 0.95, 
  e_ins_ang = rbeta(n_size, 26.00386, 8.667955), 
  e_ins_pvd = rbeta(n_size, 26.00386, 8.667955), 
  e_ins_fail = rbeta(n_size, 26.00386, 8.667955),
  
  # Effects of strategies
  e_screen_on_diag = rnorm(n_size, 1.5, 0.1), # rr
  e_screen_on_util = rnorm(n_size, 2, 0.1), # rr
  e_comp_on_util = rnorm(n_size, 2, 0.1), # rr
  e_tx = rlnorm(n_size, -0.9877235, 0.3779755), # p
  
  # Probability of diabetes diagnosis
  p_diagnosis = rbeta(n_size, 8.241926, 14.05766), # Assumption
  
  # HEF policy 
  p_hef_enrollment = rbeta(n_size, 5.5, 1.833333), 
  p_hef_utilization = rep(1, n_size), 
  hef_threshold = rnorm(n_size, 0.25, 0.05), 
  coverage = rbeta(n_size, 17.1, 0.9), 
  
  # Disability weights 
  dw_uncomplicated = rnorm(n_size, mean = 0.049, sd = 0.0115), 
  dw_nephropathy = rnorm(n_size, mean = 0.569, sd = 0.079), 
  dw_retinopathy = rnorm(n_size, mean = 0.184, sd = 0.037), 
  dw_neuropathy = rnorm(n_size, mean = 0.133, sd = 0.027), 
  dw_angina = rnorm(n_size, mean = 0.080, sd = 0.033), 
  dw_pvd = rnorm(n_size, mean = 0.014, sd = 0.0055), 
  dw_mi = rnorm(n_size, mean = 0.432, sd = 0.0735), 
  dw_stroke = rnorm(n_size, mean = 0.588, sd = 0.078), 
  dw_failure = rnorm(n_size, mean = 0.179, sd = 0.036)
) %>% 
  mutate(
    p_nephro_to_retino = p_retino * rnorm(n_size, 5.01, 0.18),
    e_oad_retino = e_oad,
    e_ins_retino = e_ins
    
  )

careseeking <- data.frame(
  p_utilization_hc = rbeta(n_size, 9.755472, 107.7964),
  p_utilization_cpa1 = rbeta(n_size, 10.94302, 574.0505),
  p_utilization_cpa2 = rbeta(n_size, 11.86287, 557.8433),
  p_utilization_cpa3 = rbeta(n_size, 7.304759, 190.2188)
)

p_treatment <- data.frame(
  p_oad = rbeta(n_size, 15.35063, 53.17897),
  p_ins = rbeta(n_size, 11.34648, 656.0935)
)

p_utilization <- data.frame(
  p_op_hef = rbeta(n_size, 0.9166912, 4.412909),
  p_op_nohef = rbeta(n_size, 0.4202172, 3.171383),
  p_hospitalization = rbeta(n_size, 8.85, 581.15)
)

costs <- data.frame(
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
  monthly_subsistence = rep(47.96339, n_size)
) 


# 3. Run sensitivity analysis ---------------------------------------------

# Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(8)

# Run function for each parameter set
n_parameter_sets = n_size
psa_results <- foreach::foreach(i_parameters = 1:n_parameter_sets) %dopar% {
  run_and_compute_markov(
    costs = costs[i_parameters,],
    parameters = parameters[i_parameters,],
    careseeking = careseeking[i_parameters,],
    p_treatment = p_treatment[i_parameters,],
    p_utilization = p_utilization[i_parameters,]
  )
}
  






