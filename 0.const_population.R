library(readstata13)


### Demographic characteristics
## Age distribution
age_data <- read.dta13("input/dhskhr2014_age_distribution.dta") # Cambodia DHS, 2014
age_data <- age_data %>% filter(!is.na(age) & age >= 25 & age <=69)
age_dist <- age_data$age

## Sex distribution
p_female = 0.52

## Income distribution
avg_income <- 963.86 * 12

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

### Age- & sex-specific diabetes prevalence
ihme_dm_khm_prev <- read_csv("input/IHME-GBD_2017_DATA-19f45038-1_prevalence.csv")
age_specific_dm_prevalence_male <- ihme_dm_khm_prev %>% 
  filter(measure_name == "Prevalence", sex_name == "Male", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_prevalence_male <- rbind(age_specific_dm_prevalence_male)
age_specific_dm_prevalence_male$age <- c(25:69)
age_specific_dm_prevalence_male <- data.frame(age_specific_dm_prevalence_male)
age_specific_dm_prevalence_female <- ihme_dm_khm_prev %>% 
  filter(measure_name == "Prevalence", sex_name == "Female", age_name != "All Ages") %>% 
  select(sex_name, age_name, val, upper, lower) %>% 
  slice(rep(1:n(), each = 5)) 
age_specific_dm_prevalence_female <- rbind(age_specific_dm_prevalence_female)
age_specific_dm_prevalence_female$age <- c(25:69)
age_specific_dm_prevalence_female <- data.frame(age_specific_dm_prevalence_female)
dm_prevalence <- rbind(age_specific_dm_prevalence_male, age_specific_dm_prevalence_female)

### Careseeking probabilities (at which facility, scaled based on overall probability of care utilization) (DHS 2014, Care seeking for first treatment, total population)
p_utilization_hc <- 0.114
p_utilization_cpa1 <- 0.025
p_utilization_cpa2 <- 0.030
p_utilization_cpa3 <- 0.042
p_provider <- c(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3) / sum(p_utilization_hc, p_utilization_cpa1, p_utilization_cpa2, p_utilization_cpa3)

## Subsistence expenditure
subsistence <- 47.96339 * 12 # 196K KHR (CSES 2017) / 4086.45 KHR per USD (2 Oct 2019); average monthly value per capital, 2017

## HEF defaults
p_hef_enrollment <- 0.75
p_hef_utilization <- 1 # 0.16 ## HEF utilization at point of care
coverage <- 1

