##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 1.25 Markov model 
#### Copyright Isabelle Feldhaus
#### 6 December 2019
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

library(tidyverse)
library(markovchain) 
library(readstata13)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')

# 1. Disability weights ---------------------------------------------------
# Source: GBD 2017

# Diagnosed, uncomplicated DM 
dw_uncomplicated <- 0.049 # 0.049 (0.031 - 0.072)

# Nephropathy
dw_nephropathy_stage5 <- 0.569 # 0.569 (0.389 - 0.727)
dw_nephropathy_stage4 <- 0.104 # 0.104 (0.070 - 0.147)
dw_nephropathy_stage3 <- 0.004 # 0.004 (0.001 - 0.008)

# Retinopathy
dw_retinopathy <- 0.184 # 0.184 (0.125 - 0.258)

# Distal symmetric neuropathy
dw_neuropathy <- 0.133 # 0.133 (0.089 - 0.187)

# Angina pectoris
dw_angina_mild <- 0.033 # 0.033 (0.020 - 0.052) 
dw_angina_moderate <- 0.080 # 0.080 (0.052 - 0.113)

# Peripheral vascular disease
dw_pvd <- 0.014 # 0.014 (0.007 - 0.025)

# Myocardial infarction
dw_mi <- 0.432 # 0.432 (0.288 - 0.579)

# Stroke
dw_stroke_level1 <- 0.019 # 0.019 (0.010 - 0.032)
dw_stroke_level2 <- 0.070 # 0.070 (0.046 - 0.099)
dw_stroke_level3 <- 0.316 # 0.316 (0.206 - 0.437)
dw_stroke_level4 <- 0.552 # 0.552 (0.377 - 0.707)
dw_stroke_level5 <- 0.588 # 0.588 (0.411 - 0.744)

# Heart failure
dw_failure_mild <- 0.041 # 0.041 (0.026 - 0.062)
dw_failure_moderate <- 0.072 # 0.072 (0.047 - 0.103)
dw_failure_severe <- 0.179 # 0.179 (0.122 - 0.251)


# 2. Population characteristics -------------------------------------------

### Demographic characteristics
## Age distribution
age_data <- read.dta13("input/dhskhr2014_age_distribution.dta") # Cambodia DHS, 2014
age_data <- age_data %>% filter(!is.na(age) & age >= 25 & age <=69)
age_dist <- age_data$age

## Sex distribution
p_female = 0.52

## Income distribution
avg_income <- 963.86 * 12
income_dist <- rgamma(n = 100000, shape = 0.5, scale = avg_income) # Distribution derived from average income per month 2017, CSES (1960 thousand KHR / 481.93 USD 2019)
income_quintile_thresholds <- qgamma(c(0.2, 0.4, 0.6, 0.8), shape = 0.5, scale = avg_income)

### Age- and sex-specific all-cause mortality
# Source: WHO GHO for Cambodia, 2016
age_specific_all_cause_mortality_male_2016 <- data.frame(
  mr = c(
    rep(0.012, 5), 
    rep(0.016, 5), 
    rep(0.021, 5), 
    rep(0.026, 5), 
    rep(0.032, 5),
    rep(0.040, 5),
    rep(0.063, 5),
    rep(0.111, 5),
    rep(0.163, 5)
  )
) %>% 
  mutate(age = 25:69,
         sex = "Male")
age_specific_all_cause_mortality_female_2016 <- data.frame(
  mr = c(
    rep(0.007, 5),
    rep(0.009, 5),
    rep(0.013, 5),
    rep(0.016, 5),
    rep(0.021, 5),
    rep(0.029, 5),
    rep(0.044, 5),
    rep(0.079, 5),
    rep(0.129, 5)
  )
) %>% 
  mutate(age = 25:69,
         sex = "Female")
all_cause_mortality <- rbind(age_specific_all_cause_mortality_male_2016, age_specific_all_cause_mortality_female_2016)

### Age- & sex-specific diabetes incidence
# Source: IHME GBD 2017 Cambodia
ihme_dm_khm <- read_csv("input/IHME-GBD_2017_DATA-0e7b15fe-1.csv")
age_specific_dm_incidence_male <- ihme_dm_khm %>% 
  filter(measure_name == "Incidence", sex_name == "Male", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_incidence_male <- rbind(age_specific_dm_incidence_male)
age_specific_dm_incidence_male$age <- c(25:84)
age_specific_dm_incidence_male <- data.frame(age_specific_dm_incidence_male)
age_specific_dm_incidence_female <- ihme_dm_khm %>% 
  filter(measure_name == "Incidence", sex_name == "Female", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_incidence_female <- rbind(age_specific_dm_incidence_female)
age_specific_dm_incidence_female$age <- c(25:84)
age_specific_dm_incidence_female <- data.frame(age_specific_dm_incidence_female)
dm_incidence <- rbind(age_specific_dm_incidence_male, age_specific_dm_incidence_female)

### Age- & sex-specific mortality due to diabetes
age_specific_dm_mortality_male <- ihme_dm_khm %>% 
  filter(measure_name == "Deaths", sex_name == "Male", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_mortality_male <- rbind(age_specific_dm_mortality_male)
age_specific_dm_mortality_male$age <- c(25:84)
age_specific_dm_mortality_male <- data.frame(age_specific_dm_mortality_male)
age_specific_dm_mortality_female <- ihme_dm_khm %>% 
  filter(measure_name == "Deaths", sex_name == "Female", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_mortality_female <- rbind(age_specific_dm_mortality_female)
age_specific_dm_mortality_female$age <- c(25:84)
age_specific_dm_mortality_female <- data.frame(age_specific_dm_mortality_female)
dm_mortality <- rbind(age_specific_dm_mortality_male, age_specific_dm_mortality_female)

### Careseeking probabilities (at which facility, scaled based on overall probability of care utilization) (DHS 2014, Care seeking for first treatment, total population)
p_utilization_hc <- 0.114
p_utilization_cpa1 <- 0.025
p_utilization_cpa2 <- 0.030
p_utilization_cpa3 <- 0.042
p_provider <- c(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3) / sum(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3)


# 3. Costs ----------------------------------------------------------------

## Lognormal computation: 
# m = exp(mu + sigma^2 / 2)
# v = exp(2*mu + sigma^2)(exp(sigma^2) - 1)
## Calculating lognormal distributions
lognormal_fxn <- function(m, z_alpha, alpha) {
  B <- -2 * qnorm(alpha)
  C <- 2*(log(z_alpha) - log(m))
  sigma <- (-B + c(-1,1)*sqrt(B^2 - 4*C)) / 2
  sigma <- sigma[sigma > 0 & !is.na(sigma)]
  mu <- log(m) - sigma^2 / 2
  return(cbind(mu=mu, sigma=sigma)[1,])
}

## Calculating beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## Discount rate
dr <- 0.03

### Financial risk protection
## Catastrophic health expenditure (CHE)
che_threshold = c(0.1, 0.25, 0.4)
# che <- function(oop_individual, income) {
#   if (income < 0) {
#     che_10 <- NA
#     che_25 <- NA
#     che_40 <- NA
#   } else {
#     che_10 <- ifelse(oop_individual > income * che_threshold[1], 1, 0)
#     che_25 <- ifelse(oop_individual > income * che_threshold[2], 1, 0)
#     che_40 <- ifelse(oop_individual > income * che_threshold[3], 1, 0)
#   }
#   
#   return(list("che10" = che_10, "che25" = che_25, "che40" = che_40))
# }

che <- function(oop_individual, income) {
  che_10 <- ifelse(income < 0, NA, 
                   ifelse(oop_individual > income * che_threshold[1], 1, 0))
  che_25 <- ifelse(income < 0, NA,
                   ifelse(oop_individual > income * che_threshold[2], 1, 0))
  che_40 <- ifelse(income < 0, NA, 
                   ifelse(oop_individual > income * che_threshold[3], 1, 0))
  
  return(list("che10" = che_10, "che25" = che_25, "che40" = che_40))
}

## Subsistence expenditure
subsistence <- 47.96339 * 12 # 196K KHR (CSES 2017) / 4086.45 KHR per USD (2 Oct 2019); average monthly value per capital, 2017

## Poverty cases
poverty_line = 52.21 * 12 # 7,112 KHR per day (* 30 days) / 4086.45 KHR per USD (2 Oct 2019)
pov <- function(oop_individual, income) {
  if (income > poverty_line) {
    ifelse(income - oop_individual <= poverty_line, 1, 0)
  } else {
    NA
  }
}

# 4. Coverage parameters --------------------------------------------------

## HEF eligibility and enrollment
hef_threshold <- quantile(income_dist, 0.2)[[1]] # Income threshold at which individuals are eligible for HEF (USD); ~ 20th percentile of (gamma) income distribution ($28-31 USD) used in the population income distribution parameter defined below
p_hef_enrollment <- 0.75

## Benefits of coverage
p_coverage <- c(0.8, 0.9, 1.0) # If individual is under HEF, 80% of expenditures are covered.
coverage <- p_coverage[1]


# 5. Transition probabilities ---------------------------------------------

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


# 6. States and strategies ------------------------------------------------

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
  init_states_KHR <- c((1 - prevalence * (1 + p_diet + p_oad + prev_nephro + prev_retino + prev_neuro)), 
                       prevalence, prevalence * p_diet, prevalence * p_oad, 0, 
                       prevalence * prev_nephro, prevalence * prev_retino, prevalence * prev_neuro, 
                       0, 0, 0, 0, 0, 0, 0) # CVD-related data (among diabetics) not available
  return(init_states_KHR)
}

## Strategies to simulate
strategyNames <- c("base", "screen_only", "tx_only", "comp_only", "screen_tx", "tx_comp", "screen_tx_comp")

# Naming strategies in final tables (alphabetical)
setNames <- function(x) {
  colnames(x) <- c("base", "comp_only", "screen_only", "screen_tx", "screen_tx_comp", "tx_comp", "tx_only")
  rownames(x) <- c("base", "comp_only", "screen_only", "screen_tx", "screen_tx_comp", "tx_comp", "tx_only")
  return(x)
}


# 8. Model implementation ------------------------------------------------

## Profiling (see below too)
# Rprof('profiling_model')

# Population / iterations
# n_population <- as.integer(1e5) # total: 16e6
n_population <- 1000

# Time horizon
n_cycles = 45

# Load required functions 
source("1.22_model_cost.R")
source("1.20_model_dalys.R")

### Model

## Set up parallel computing
foreach::getDoParWorkers()
doMC::registerDoMC(4)

# Determine sample groups for parallel execution, change first line only
indicative_parallel_samples <- 20
sample_size <- (n_population - 1) %/% indicative_parallel_samples + 1
n_parallel_samples <- (n_population - 1) %/% sample_size + 1

# Empty lists to store final results
all_results <- list()
all_ichar <- list()
population_results <- list()

## Defining the populations
# Set characteristics
ichar <- data.frame(matrix(ncol = 9))
colnames(ichar) <- c("id", "age_start", "sex", "income", "disposable_income", "income_quintile", "hef", "hef_utilization", "init_state")

cat('Simulating population of', n_population, 'people')
for (i_person in 1:n_population) {
  
  if (i_person %% 1000 == 0) {
    cat('.')
  }
  
  # Set individual characteristics
  age <- sample(age_dist, 1, replace = T)
  sex <- rbinom(n = 1, size = 1, p = p_female)
  sex <- ifelse(sex == 0, "Male", "Female") # Relabel for `look_up`
  income <- rgamma(n = 1, shape = 0.5, scale = avg_income) 
  disposable_income <- income - subsistence
  init_state <- sample(stateNames, size = 1, prob = initialStates_KHR(age, sex))
  
  # Determine HEF status & utilization behavior
  hef <- ifelse(income < hef_threshold, rbinom(n = 1, size = 1, p = p_hef_enrollment), 0) # 0 = non-HEF, 1 = HEF
  if (hef == 1) {
    hef_utilization <- rbinom(n = 1, size = 1, p = p_hef_utilization)
  } else {
    hef_utilization <- NA
  }
  
  income_quintile = ifelse(income > income_quintile_thresholds[4], 5,
                           ifelse(income <= income_quintile_thresholds[4] & income > income_quintile_thresholds[3], 4,
                                  ifelse(income <= income_quintile_thresholds[3] & income > income_quintile_thresholds[2], 3,
                                         ifelse(income <= income_quintile_thresholds[2] & income > income_quintile_thresholds[1], 2, 1))))
  
  # Store individual characteristics 
  ichar[i_person, 1] <- i_person
  ichar[i_person, 2] <- age
  ichar[i_person, 3] <- sex
  ichar[i_person, 4] <- income
  ichar[i_person, 5] <- disposable_income
  ichar[i_person, 6] <- income_quintile
  ichar[i_person, 7] <- hef
  ichar[i_person, 8] <- hef_utilization
  ichar[i_person, 9] <- init_state
  ichar <- data.frame(ichar)
}
cat(' Done\n')

n_populations <- 5
# Splitting into n_populations - won't do, in case we want different ichars for the psa populations
# number_of_people_per_sample <- nrow(ichar) / n_populations
# list_of_populations <- split(ichar, rep(1:n_populations, nrow(ichar) / n_populations))

## Defining parameters 

parameters_df <- data.frame(
  ## Direct medical costs; Source: Flessa & Zembok, 2014; adjusted to 2019 USD
  c_screening = rlnorm(n_populations, -0.231791, 0.8088278), # 1.10 # i.e. screening; unit cost FPG-test
  c_laboratory = rlnorm(n_populations, -0.02394995,0.8573677), # 1.41 # diagnosis cost
  c_oad = rlnorm(n_populations, 3.2449579, 0.4055468), # 27.86 # Herman, 2005: $50 annual cost
  c_ins = rlnorm(n_populations, 4.81576532, 0.1945663), # 125.80
  c_outpatient_cpa3 = rlnorm(n_populations, 3.791499, 0.05900272), # 44.40 # Maximum prices
  c_outpatient_cpa2 = rlnorm(n_populations, 1.754165, 0.407934), # 6.28 
  c_outpatient_cpa1 = rlnorm(n_populations, 2.298063, 0.2684041), # 10.32
  c_outpatient_hc = rlnorm(n_populations, 1.299967, 0.4962679), # 4.15 
  c_hospitalization_cpa3 = rlnorm(n_populations, 3.703311, 0.1148535), # 40.85
  c_hospitalization_cpa2 = rlnorm(n_populations, 3.378472, 0.1148535), # 29.52
  c_hospitalization_cpa1 = rlnorm(n_populations, 4.0832387, 0.1148535), # 59.73
  c_hospitalization_hc = rlnorm(n_populations, 1.479141, 0.407934), # 4.77
  
  ## Direct medical costs associated with each complication
  # See lognormal distribution above
  c_nephro = rlnorm(n_populations, 8.733453, 0.2191612), # 6358 # Hemodialysis, systematic review of LMICs; Int$ 3,424 to Int$ 42,785 overall, but used Sri Lanka numbers: Int$ 5,869–8,804; Mushi et al., 2015 -- adjusted for inflation
  c_retino = rlnorm(n_populations, 5.7928000, 0.1148535), # 330.1 # Cost of intravitreal injection in Indonesia, Sasongko et al., 2019
  c_neuro = rlnorm(n_populations, 6.862624, 0.0299347), # 2.62 * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_ang = rlnorm(n_populations, 7.855310, 0.0299347), # 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_pvd = rlnorm(n_populations, 7.855310, 0.0299347), # 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_mi = rlnorm(n_populations, 7.855310, 0.0299347), # 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_stroke = rlnorm(n_populations, 6.862624, 0.0299347), # 2.62 * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
  c_fail = rlnorm(n_populations, 7.855310, 0.0299347), # 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
  
  ### Individual expenditures
  ## Indirect (non-medical) costs
  # See lognormal distribution above
  oop_transport_op = rlnorm(n_populations, -0.1379559, 0.3303765), # 0.92
  oop_transport_hosp = rlnorm(n_populations, 2.452990, 0.06806667), # 11.65
  
  ## Effects of strategies
  e_screen_on_diag = rnorm(n_populations, 1.5, 0.1),
  e_screen_on_util = rnorm(n_populations, 2, 0.1),
  e_comp = rnorm(n_populations, 2, 0.1),
  e_tx = rlnorm(n_populations, -0.9877235, 0.3779755),
  
  ## HEF utilization at point of care
  p_hef_utilization = rlnorm(n_populations, -1.8333334, 0.15) # 0.16
)
parameters_df$e_tx <- ifelse(parameters_df$e_tx > 1, 1, parameters_df$e_tx) # Truncate effect of treatment at 1

## Matrix of state probabilities (pointers) for each individual and cycle
person_cycle_x <- matrix(runif(n_population * n_cycles), n_population, n_cycles)

## Simulation for each population
population_counter <- 1
populations_results <- data_frame()
for (i_population in 1:n_populations) {
  
  cat(sprintf('Starting population %d of %d populations: \n', population_counter, length(list_of_populations)))
  
  # Define parameters for the population
  c_screening <- parameters_df[population_counter, 'c_screening']
  c_laboratory <- parameters_df[population_counter, 'c_laboratory']
  c_oad <- parameters_df[population_counter, 'c_oad']
  c_nephro <- parameters_df[population_counter, 'c_nephro']
  c_retino <- parameters_df[population_counter, 'c_retino']
  c_neuro <- parameters_df[population_counter, 'c_neuro']
  c_ang <- parameters_df[population_counter, 'c_ang']
  c_pvd <- parameters_df[population_counter, 'c_pvd']
  c_mi <- parameters_df[population_counter, 'c_mi']
  c_stroke <- parameters_df[population_counter, 'c_stroke']
  c_fail <- parameters_df[population_counter, 'c_fail']
  oop_transport_op <- parameters_df[population_counter, 'oop_transport_op']
  oop_transport_hosp <- parameters_df[population_counter, 'oop_transport_hosp']
  
  outpatient_costs <- c(
    "HC" = parameters_df[population_counter, 'c_outpatient_hc'],
    "CPA1" = parameters_df[population_counter, 'c_outpatient_cpa1'],
    "CPA2" = parameters_df[population_counter, 'c_outpatient_cpa2'],
    "CPA3" = parameters_df[population_counter, 'c_outpatient_cpa3']
  )
  
  n_days <- 5
  hospitalization_costs <- c(
    "HC" = parameters_df[population_counter, 'c_hospitalization_hc'],
    "CPA1" = parameters_df[population_counter, 'c_hospitalization_cpa1'],
    "CPA2" = parameters_df[population_counter, 'c_hospitalization_cpa2'],
    "CPA3" = parameters_df[population_counter, 'c_hospitalization_cpa3']
  ) * n_days
  
  e_screen_on_diag <- parameters_df[population_counter, 'e_screen_on_diag']
  e_screen_on_util <- parameters_df[population_counter, 'e_screen_on_util']
  e_comp <- parameters_df[population_counter, 'e_comp']
  e_tx <- parameters_df[population_counter, 'e_tx']
  
  p_hef_utilization <- parameters_df[population_counter, 'p_hef_utilization']
  
  ## Simulation for each strategy
  for (strategy in strategyNames) {
    
    cat("Starting strategy:", strategy)
    
    ## Prepare iteration
    # Set seed and time computation
    set.seed(24)
    start_time <- Sys.time()
  
    ## Instantiate the data frame to store the results
    n_rows <- sample_size * n_cycles
    results <- data.frame(
      strategy = character(),
      id = integer(),
      cycle = integer(),
      state = character(),
      cost_diagnosis = numeric(),
      cost_tx = numeric(),
      cost_complications = numeric(),
      cost_total = numeric(),
      dalys_dm = numeric(),
      dalys_other = numeric(),
      dalys_total = numeric(),
      oop_medical = numeric(),
      oop_nonmedical = numeric(),
      oop_total = numeric(),
      che10_result = integer(),
      che25_result = integer(),
      che40_result = integer(),
      che10_disposable_result = integer(),
      che25_disposable_result = integer(),
      che40_disposable_result = integer(),
      pov_result = integer(),
      pov_disposable_result = integer(),
      stringsAsFactors = FALSE
    )
    # Instantiate all the rows with NA
    results[n_rows, 'id'] <- NA
    
    cat(sprintf(', running %d samples of size %d... ', length(list_of_populations), number_of_people_per_sample))
    
    ## Perform iteration for each iteration group
    model_results <- foreach::foreach(i_person = i_population$id) %dopar% {
      
      ## Set counter for the sample, iterating along people and cycles
      j_counter <- 1
      
      ## Initialize the transition matrix
      tmat <- matrix(NA, nrow = 15, ncol = 15)
      colnames(tmat) <- stateNames
      rownames(tmat) <- stateNames
      
      row <- i_population %>% filter(id == i_person)  
      id <- row$id
      age <- as.numeric(row$age_start)
      sex <- row$sex
      income <- as.numeric(row$income)
      disposable_income <- as.numeric(row$disposable_income)
      hef <- row$hef
      hef_utilization <- row$hef_utilization
      init_state <- row$init_state
      
      ## Set relevant probabilities
      is_strategy_screen <- strategy == "screen_only" | strategy == "screen_tx" | strategy == "screen_tx_comp"
      is_strategy_tx <- !(strategy == "base" | strategy == "screen_only" | strategy == "comp_only")
      is_strategy_comp <- strategy == "comp_only" | strategy == "tx_comp" | strategy == "screen_tx_comp"

      #  Diagnosis; Source: Flessa & Zembok 2014
      p_diag <- 0.3696 * ifelse(is_strategy_screen & hef, e_screen_on_diag, 1)
      p_undiag <- 1 - p_diag
      
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
        source("_tmat_construction.R", local = TRUE)
        
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
        results[j_counter, 1] <- strategy
        results[j_counter, 2] <- id
        results[j_counter, 3] <- cycle
        results[j_counter, 4] <- end_state
        results[j_counter, 5] <- costs$cost_diagnosis
        results[j_counter, 6] <- costs$cost_tx
        results[j_counter, 7] <- costs$cost_complications
        results[j_counter, 8] <- costs$cost_total
        results[j_counter, 9] <- dalys$dalys_dm
        results[j_counter, 10] <- dalys$dalys_other
        results[j_counter, 11] <- dalys$dalys_total
        results[j_counter, 12] <- costs$oop_medical
        results[j_counter, 13] <- costs$oop_nonmedical
        results[j_counter, 14] <- costs$oop_total
        results[j_counter, 15] <- che_result$che10
        results[j_counter, 16] <- che_result$che25
        results[j_counter, 17] <- che_result$che40
        results[j_counter, 18] <- che_result_disposable$che10
        results[j_counter, 19] <- che_result_disposable$che25
        results[j_counter, 20] <- che_result_disposable$che40
        results[j_counter, 21] <- pov_result
        results[j_counter, 22] <- pov_result_disposable
        
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
        
      return(list("results" = results %>% filter(!is.na(strategy))))
    }
    
    # End time computation
    end_time <- Sys.time()
    time <- end_time - start_time
    cat("Done in", time, '\n')
    
    # Combine results into single dataframes (results, ichar)
    results_per_sample <- lapply(1:length(list_of_populations), function(x) model_results[[x]]$results)
    results <- do.call(rbind, results_per_sample)
    all_results[[strategy]] <- results
    results_per_population <- do.call(rbind, all_results)
    
    # Compute outcomes for each strategy (by income quintile and sex)
    income_sex <- i_population %>% select(id, income_quintile, sex)
    results_per_population <- merge(results_per_population, income_sex, by = "id")
    
    costs <- results_per_population %>% 
      group_by(strategy, income_quintile, sex) %>% 
      summarize(
        cost_diagnosis = sum(cost_diagnosis),
        cost_tx = sum(cost_tx),
        cost_complications = sum(cost_complications), 
        cost_total = sum(cost_total)
      ) 
    
    health <- results_per_population %>% 
      group_by(strategy, income_quintile, sex) %>% 
      summarize(
        dm_deaths = sum(state == "dm_death"),
        dm_dalys = sum(dalys_dm)
      )
    
    frp <- results_per_population %>% 
      group_by(strategy, income_quintile, sex) %>% 
      summarize(
      che10 = sum(che10_result),
      che25 = sum(che25_result),
      che40 = sum(che40_result),
      pov = sum(pov_result, na.rm = T),
      che10_disposable = sum(che10_disposable_result),
      che25_disposable = sum(che25_disposable_result),
      che40_disposable = sum(che40_disposable_result),
      pov_disposable = sum(pov_disposable_result, na.rm = T)
    )
    
    outcomes <- list("costs" = costs, "health" = health, "frp" = frp)
    
    # Store results in data frame
    population_results[[population_counter]] <- outcomes
  }
  
  ### Analysis for population i_pop
  
  ### Saving results of population i_pop
  
  population_counter <- population_counter + 1
}

# Combine results -- may be too large
all_results_final <- do.call(rbind, all_results)

## Profiling
# Rprof(NULL)
# profvis::profvis(prof_input = 'profiling_model')

## Export flat files
save_df_to_csv <- function(df, title, output_dir = 'output') {
  filename <- sprintf('%s/%s_%dL_%s.csv', output_dir, title, nrow(df), format(Sys.time(), '%Y-%m-%d_%H%M'))
  readr::write_csv(df, filename)
}
save_df_to_csv(all_results_final, "results")
save_df_to_csv(ichar, "ichar")

# Stop cluster
# parallel::stopCluster(cl)

