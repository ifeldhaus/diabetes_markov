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

## Load static parameters 
source('0.const_population.R') # Population characteristics

# 1. Parameter Distributions ---------------------------------------------

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
  dw_uncomplicated <- rep(0.049, n_size), # 0.049 (0.031 - 0.072)
  dw_nephropathy_stage5 <- rep(0.569, n_size), # 0.569 (0.389 - 0.727)
  dw_retinopathy <- rep(0.184, n_size), # 0.184 (0.125 - 0.258)
  dw_neuropathy <- rep(0.133, n_size), # 0.133 (0.089 - 0.187)
  dw_angina <- rep(0.080, n_size), # 0.080 (0.052 - 0.113)
  dw_pvd <- rep(0.014, n_size), # 0.014 (0.007 - 0.025)
  dw_mi <- rep(0.432, n_size), # 0.432 (0.288 - 0.579)
  dw_stroke_level5 <- rep(0.588, n_size), # 0.588 (0.411 - 0.744)
  dw_failure_severe <- rep(0.179, n_size) # 0.179 (0.122 - 0.251)
  
  # dw_uncomplicated = rnorm(n_size, mean = 0.049, sd = 0.0115), 
  # dw_nephropathy = rnorm(n_size, mean = 0.569, sd = 0.079), 
  # dw_retinopathy = rnorm(n_size, mean = 0.184, sd = 0.037), 
  # dw_neuropathy = rnorm(n_size, mean = 0.133, sd = 0.027), 
  # dw_angina = rnorm(n_size, mean = 0.080, sd = 0.033), 
  # dw_pvd = rnorm(n_size, mean = 0.014, sd = 0.0055), 
  # dw_mi = rnorm(n_size, mean = 0.432, sd = 0.0735), 
  # dw_stroke = rnorm(n_size, mean = 0.588, sd = 0.078), 
  # dw_failure = rnorm(n_size, mean = 0.179, sd = 0.036)
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

### Function: `run_diabetes_markov`
# cbind(parameters, careseeking, costs, p_treatment, p_utilization)
apply(modified_parameters_matrix, 1, function(parameters) {
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
  
  results <- run_model(population, person_cycle_x)
  outcomes <- compute_outcomes(population, merge(results, population, by = 'id'), hef_threshold = hef_quantile)
  
}




