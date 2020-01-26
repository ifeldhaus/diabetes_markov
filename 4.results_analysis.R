##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 2.6 Analysis of results
#### Copyright Isabelle Feldhaus
#### 24 January 2020
##################################


# 0. Preparation ----------------------------------------------------------

library(tidyverse)
library(scales)
source('helpers.R')
source('1.outcomes.R')
source('1.model_pop.R')
library(wesanderson)
library(ggsci)
library(xtable)

## Load results
results_threshold20_coverage80 <- readr::read_csv('output/results_23985420L_2020-01-25_1910.csv')
results_threshold30_coverage80 <- readr::read_csv('output/results_23985588L_2020-01-25_1927.csv')

# Load population and assign HEF status
population <- readr::read_csv("output/ichar_200000L_2020-01-25_1855.csv")
person_cycle_x <- load_object_from_rdata('output/person_cycle_x_200000L_2020-01-25_1855.RData')
population20 <- add_hef_to_pop(population, person_cycle_x, hef_quantile = 0.2) %>% mutate(threshold = "20%")
population30 <- add_hef_to_pop(population, person_cycle_x, hef_quantile = 0.3) %>% mutate(threshold = "30%")

## Define income quintile thresholds and add variable to results dataframe 
income_quintile_thresholds <- quantile(population$income, c(0.2, 0.4, 0.6, 0.8))
population20 <- population20 %>% mutate(
  income_quintile = ifelse(income > income_quintile_thresholds[4], 5,
                           ifelse(income <= income_quintile_thresholds[4] & income > income_quintile_thresholds[3], 4,
                                  ifelse(income <= income_quintile_thresholds[3] & income > income_quintile_thresholds[2], 3,
                                         ifelse(income <= income_quintile_thresholds[2] & income > income_quintile_thresholds[1], 2, 1))))
)
population20$income_quintile <- ordered(population20$income_quintile, levels = 1:5)

population30 <- population30 %>% mutate(
  income_quintile = ifelse(income > income_quintile_thresholds[4], 5,
                           ifelse(income <= income_quintile_thresholds[4] & income > income_quintile_thresholds[3], 4,
                                  ifelse(income <= income_quintile_thresholds[3] & income > income_quintile_thresholds[2], 3,
                                         ifelse(income <= income_quintile_thresholds[2] & income > income_quintile_thresholds[1], 2, 1))))
)
population30$income_quintile <- ordered(population30$income_quintile, levels = 1:5)

## HEF eligibility income thresholds 
hef_threshold20 <- quantile(population$income, 0.2)[[1]]
hef_threshold30 <- quantile(population$income, 0.3)[[1]]

## Add sex and income quintile to results data frame and filter by HEF eligible population
results20_80 <- merge(results_threshold20_coverage80, population20 %>% select(id, sex, income, disposable_income, income_quintile, hef, threshold), by = "id") %>% 
  filter(income < hef_threshold20) %>% 
  mutate(coverage = "80%")
results30_80 <- merge(results_threshold30_coverage80, population30 %>% select(id, sex, income, disposable_income, income_quintile, hef, threshold), by = "id") %>% 
  filter(income < hef_threshold30) %>% 
  mutate(coverage = "80%")

## Setting HEF coverage to 100%
results20_100 <- results20_80 %>%
  mutate(oop_medical = ifelse(strategy != "base" & hef == 1, 0, oop_medical),
         oop_total = oop_medical + oop_nonmedical,
         che10_result = che(oop_total, income)$che10,
         che25_result = che(oop_total, income)$che25,
         che40_result = che(oop_total, income)$che40,
         che10_disposable_result = che(oop_total, disposable_income)$che10,
         che25_disposable_result = che(oop_total, disposable_income)$che25,
         che40_disposable_result = che(oop_total, disposable_income)$che40,
         pov_result = pov(oop_total, income),
         pov_disposable_result = pov(oop_total, disposable_income), 
         coverage = "100%")

results30_100 <- results30_80 %>%
  mutate(oop_medical = ifelse(strategy != "base" & hef == 1, 0, oop_medical),
         oop_total = oop_medical + oop_nonmedical,
         che10_result = che(oop_total, income)$che10,
         che25_result = che(oop_total, income)$che25,
         che40_result = che(oop_total, income)$che40,
         che10_disposable_result = che(oop_total, disposable_income)$che10,
         che25_disposable_result = che(oop_total, disposable_income)$che25,
         che40_disposable_result = che(oop_total, disposable_income)$che40,
         pov_result = pov(oop_total, income),
         pov_disposable_result = pov(oop_total, disposable_income), 
         coverage = "100%")

# Inflation factor
inf.fct <- 16e6 / nrow(population)

## Quick look
# Prevalence of diabetes over 45 years among HEF eligible population
n_diabetes <- results20_80 %>% filter(state != "healthy") %>% filter(state != "other_death") %>% group_by(id) %>% summarise(n = n()) %>% nrow() * inf.fct
n_diabetes_sex <- results20_80 %>% filter(state != "healthy") %>% filter(state != "other_death") %>% group_by(id, sex) %>% summarise(n = n()) %>% 
  ungroup %>% group_by(sex) %>% summarise(n_diabetes = n() * inf.fct)
prevalence <- (n_diabetes / inf.fct) / nrow(population)

# All results
results <- rbind(results20_80, results20_100, results30_80, results30_100)

# HEF eligible by group
n_hef_eligible_income20 <- population20 %>% group_by(income_quintile) %>% summarise(n_hef_eligible = sum(hef) * inf.fct)
n_hef_eligible_income30 <- population30 %>% group_by(income_quintile) %>% summarise(n_hef_eligible = sum(hef) * inf.fct)
n_hef_eligible_sex20 <- population20 %>% group_by(sex) %>% summarise(n_hef_eligible = sum(hef) * inf.fct)
n_hef_eligible_sex30 <- population30 %>% group_by(sex) %>% summarise(n_hef_eligible = sum(hef) * inf.fct)

# Color palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 1. State duration -------------------------------------------------------

# States by strategy
statesByStrategy <- results20_80 %>% group_by(strategy, state, income_quintile, sex) %>% summarise(n = n()) %>% 
  filter(!(strategy %in% c("comp_only", "tx_comp", "screen_tx_comp"))) %>%
  ungroup() %>% group_by(state) %>% mutate(ratio_to_base = n / n[strategy == 'base'])
statesByStrategy$state <- ordered(statesByStrategy$state, 
                                  levels = c("healthy", "diabetes", "notx", "oad", "ins", "nephro", "retino", "neuro", "ang", "pvd", "mi", "stroke", "fail", "dm_death", "other_death"),
                                  labels = c("Healthy", "Undiagnosed diabetes", "No treatment", "OAD therapy", "Insulin therapy", 
                                             "Nephropathy", "Retinopathy", "Neuropathy", "Angina pectoris", "PVD", "Myocardial infarction", "Stroke", "Heart failure", "Diabetes-related death", "Other death"))
statesByStrategy$strategy <- ordered(statesByStrategy$strategy,
                                     levels = c("base", "screen_only", "screen_tx", "screen_tx_comp", "tx_only", "tx_comp", "comp_only"),
                                     labels = c("Current standard", "Diagnosis only", "Diagnosis + treatment", "Diagnosis + treatment + complications", "Treatment only", "Treatment + complications", "Complications only"))

# Plot number of life-years spent in each state by strategy per diabetic individual in population
# By income quintile
statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
  # filter(income_quintile %in% c(1, 2)) %>%
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = income_quintile)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(strategy ~ .) +
  theme_minimal() + 
  # scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Income quintile")

# By sex
statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
  # filter(income_quintile %in% c(1, 2)) %>%
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = sex)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(strategy ~ .) +
  theme_minimal() + 
  # scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Sex")

# 2. Results: Costs --------------------------------------------------------
# Total costs --------------------------------------------------------
total_costs_plot <- results %>% 
  group_by(strategy, threshold, coverage) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
    ) %>% 
  mutate(n_hef_eligible = ifelse(threshold == "20%", sum(population20$income < hef_threshold20) * inf.fct, sum(population30$income < hef_threshold30) * inf.fct)) %>% 
  gather(key = "cost", value = "value", -c(strategy, threshold, coverage, n_hef_eligible))

## Plot
total_costs_plot$cost <- ordered(total_costs_plot$cost,
                                 levels = c("cost_diagnosis", "cost_tx", "cost_complications", "cost_total"),
                                 labels = c("Screening + \nLaboratory", "Medications", "Service Utilization \nfor Complications", "Total"))
total_costs_plot$strategy <- ordered(total_costs_plot$strategy,
                                     levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                     labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
total_costs_plot$threshold <- factor(total_costs_plot$threshold,
                                     levels = c("20%", "30%"),
                                     labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
total_costs_plot$coverage <- factor(total_costs_plot$coverage, 
                                    levels = c("80%", "100%"),
                                    labels = c("HEF coverage: 80%", "HEF coverage: 100%"))

# Stacked bar plot of types of costs per capita by strategy
ggplot(data = total_costs_plot %>% filter(cost != "Total") %>% filter(coverage == "HEF coverage: 80%"), aes(x = strategy, y = value / n_hef_eligible, fill = cost)) +
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() + 
  facet_wrap(.~threshold) +
  # scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) + 
  scale_fill_manual(values = cbp1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Cost per person eligible for HEF (2019 USD)", fill = "Type of cost") + 
  scale_y_continuous(labels = scales::comma) 
ggsave("figures/cost_per_hef_eligible.png")

# Incremental costs -------------------------------------------------------
incr_costs <- results %>% 
  group_by(strategy, threshold, coverage) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    n_hef_eligible = ifelse(threshold == "20%", sum(population20$income < hef_threshold20) * inf.fct, sum(population30$income < hef_threshold30) * inf.fct)
  )

incr_costs$strategy <- ordered(incr_costs$strategy,
                               levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                               labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_costs$threshold <- factor(incr_costs$threshold,
                               levels = c("20%", "30%"),
                               labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_costs$coverage <- factor(incr_costs$coverage, 
                              levels = c("80%", "100%"),
                              labels = c("HEF coverage: 80%", "HEF coverage: 100%"))

#### *************************************************************************************************************************
# Incremental costs per HEF eligible
ggplot(data = incr_costs %>% filter(coverage == "HEF coverage: 80%") %>% filter(strategy != "Current standard"), 
       aes(x = strategy, y = incr_cost_total / n_hef_eligible)) +
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() + 
  facet_wrap(.~threshold) +
  # scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Incremental cost per person eligible for HEF (USD 2019)", fill = "Type of cost") + 
  scale_y_continuous(labels = scales::comma) 
ggsave("figures/incr_costs_per_hef.png")

incr_costs_by_type <- incr_costs %>% 
  gather(key = "cost", value = "value", -c("strategy", "threshold", "coverage", "n_hef_eligible")) %>% 
  filter(cost %in% c("incr_cost_diagnosis", "incr_cost_tx", "incr_cost_complications")) 
incr_costs_by_type$cost <- ordered(incr_costs_by_type$cost,
                                   levels = c("incr_cost_diagnosis", "incr_cost_tx", "incr_cost_complications"),
                                   labels = c("Screening + \nLaboratory", "Medications", "Service Utilization \nfor Complications"))
ggplot(data = incr_costs_by_type %>% 
         filter(coverage == "HEF coverage: 80%") %>% 
         filter(strategy != "Current standard"), 
       aes(x = strategy, y = value / n_hef_eligible, fill = cost)) +
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() + 
  facet_wrap(.~threshold) +
  # scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Incremental cost per person eligible for HEF (USD 2019)", fill = "Type of cost") + 
  scale_y_continuous(labels = scales::comma) 
ggsave("figures/incr_costs_per_hef_by_type.png")

# Costs by income quintile ------------------------------------------------
incr_costs_by_income <- results %>% 
  group_by(strategy, threshold, coverage, income_quintile) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  group_by(threshold, coverage, income_quintile) %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base']
  ) 


incr_costs_by_income$n_hef_eligible <- ifelse(incr_costs_by_income$threshold == "20%" & incr_costs_by_income$income_quintile == 1, n_hef_eligible_income20[1, 2][[1]], 
                                              ifelse(incr_costs_by_income$threshold == "30%" & incr_costs_by_income$income_quintile == 1, n_hef_eligible_income30[1, 2][[1]], n_hef_eligible_income30[2, 2][[1]]))

# Bar plot by strategy
incr_costs_by_income$strategy <- ordered(incr_costs_by_income$strategy,
                                         levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

# Cost per person eligible for HEF by income
ggplot(data = incr_costs_by_income,
       aes(x = income_quintile, y = incr_cost_total / n_hef_eligible, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  # scale_fill_manual(values = cbp1) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Income quintile", y = "Cost per person eligible for HEF (2019 USD)", fill = "Strategy") + 
  scale_y_continuous(labels = scales::comma)

# Costs by sex ------------------------------------------------------------
incr_costs_by_sex <- results %>% 
  group_by(strategy, threshold, coverage, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  group_by(sex) %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base']
  )

incr_costs_by_sex$n_hef_eligible <- ifelse(incr_costs_by_sex$threshold == "20%" & incr_costs_by_sex$sex == "Female", n_hef_eligible_sex20[1, 2][[1]], 
                                           ifelse(incr_costs_by_sex$threshold == "20%" & incr_costs_by_sex$sex == "Male", n_hef_eligible_sex20[2, 2][[1]],
                                                  ifelse(incr_costs_by_sex$threshold == "30%" & incr_costs_by_sex$sex == "Female", n_hef_eligible_sex30[1, 2][[1]], n_hef_eligible_sex30[2, 2][[1]])))


# Bar plot by strategy
incr_costs_by_sex$strategy <- ordered(incr_costs_by_sex$strategy,
                                      levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                      labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

ggplot(data = incr_costs_by_sex %>% filter(strategy != "Current standard"),
       aes(x = strategy, y = incr_cost_total / n_hef_eligible, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Strategy", y = "Cost per person eligible for HEF (2019 USD)", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)


# 3. Results: Health gains ------------------------------------------------
# Summary of health gains -------------------------------------------------
incr_health <- results %>% 
  group_by(strategy, threshold, coverage) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    n_hef_eligible = ifelse(threshold == "20%", sum(population20$income < hef_threshold20) * inf.fct, sum(population30$income < hef_threshold30) * inf.fct)
  )

# Bar plot by strategy
incr_health$strategy <- ordered(incr_health$strategy,
                                levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
ggplot(data = incr_health %>% filter(!(strategy %in% c("Current standard", "Complications only"))) %>% filter(coverage == "80%"), 
       aes(x = strategy, y = -incr_dm_dalys, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  scale_fill_manual(values = cbp1) + 
  facet_grid(. ~ threshold) +
  theme(axis.text.x = NULL,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = NULL, y = "Diabetes-related DALYs averted", fill = "Strategy")

incr_health_state_sex <- results %>% 
  group_by(strategy, threshold, coverage, state, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  group_by(sex) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base' & state == state & sex == sex],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'& state == state & sex == sex]
  )
incr_health_state_sex$n_hef_eligible <- ifelse(incr_health_state_sex$threshold == "20%" & incr_health_state_sex$sex == "Female", n_hef_eligible_sex20[1, 2][[1]], 
                                           ifelse(incr_health_state_sex$threshold == "20%" & incr_health_state_sex$sex == "Male", n_hef_eligible_sex20[2, 2][[1]],
                                                  ifelse(incr_health_state_sex$threshold == "30%" & incr_health_state_sex$sex == "Female", n_hef_eligible_sex30[1, 2][[1]], n_hef_eligible_sex30[2, 2][[1]])))


#### *************************************************************************************************************
# DALYs by state and sex
ggplot(data = incr_health_state_sex %>% filter(coverage == "80%"), 
       aes(x = strategy, y = -incr_dm_dalys, fill = state)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() +
  facet_grid(threshold~sex) +
  # scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") + 
  scale_y_continuous(labels = scales::comma) 

# Checking proportion of DALYs averted by state 
incr_health_state_sex %>% 
  mutate(dalys_averted = ifelse(incr_dm_dalys < 0, -incr_dm_dalys, 0)) %>% 
  group_by(strategy, sex, threshold) %>% 
  mutate(dalys_averted_sum = sum(dalys_averted)) %>% 
  filter(state == "neuro" & sex == "Male") %>% 
  mutate(prop_dalys_averted = dalys_averted / dalys_averted_sum) %>% 
  select(strategy, state, sex, threshold, dalys_averted, dalys_averted_sum, prop_dalys_averted) %>% View


# Health gains by income quintile -----------------------------------------
incr_health_by_income <- results %>% 
  group_by(strategy, threshold, coverage, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  group_by(income_quintile) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base']
  )

incr_health_by_income$n_hef_eligible <- ifelse(incr_health_by_income$threshold == "20%" & incr_health_by_income$income_quintile == 1, n_hef_eligible_income20[1, 2][[1]], 
                                              ifelse(incr_health_by_income$threshold == "30%" & incr_health_by_income$income_quintile == 1, n_hef_eligible_income30[1, 2][[1]], n_hef_eligible_income30[2, 2][[1]]))

# Plots of incremental health effects
incr_health_by_income$strategy <- ordered(incr_health_by_income$strategy,
                                          levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                          labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_health_by_income$income_quintile <- ordered(incr_health_by_income$income_quintile,
                                                 levels = c(1:5))

ggplot(data = incr_health_by_income %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = strategy, y = -incr_dm_dalys, fill = income_quintile)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_grid(.~threshold) +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") +
  scale_y_continuous(labels = scales::comma) 

# Disaggregating states
incr_health_by_income_and_state <- results %>% 
  group_by(strategy, threshold, coverage, income_quintile, state) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  group_by(income_quintile) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    n_hef_eligible = ifelse(threshold == "20%", sum(income < hef_threshold20) * inf.fct, sum(income < hef_threshold30) * inf.fct)
  )


# Health gains by sex -----------------------------------------------------
incr_health_by_sex <- results %>% 
  group_by(strategy, threshold, coverage, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  group_by(sex) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base']
  )

incr_health_by_sex$n_hef_eligible <- ifelse(incr_health_by_sex$threshold == "20%" & incr_health_by_sex$sex == "Female", n_hef_eligible_sex20[1, 2][[1]], 
                                           ifelse(incr_health_by_sex$threshold == "20%" & incr_health_by_sex$sex == "Male", n_hef_eligible_sex20[2, 2][[1]],
                                                  ifelse(incr_health_by_sex$threshold == "30%" & incr_health_by_sex$sex == "Female", n_hef_eligible_sex30[1, 2][[1]], n_hef_eligible_sex30[2, 2][[1]])))

# Plot
incr_health_by_sex$strategy <- ordered(incr_health_by_sex$strategy,
                                       levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                       labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

ggplot(data = incr_health_by_sex %>% filter(!(strategy %in% c("Current standard", "Complications only"))), 
       aes(x = strategy, y = -incr_dm_dalys / n_hef_eligible, fill = sex)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)

#### *************************************************************************************************************************
## Facet wrap for income quintile and sex -- which basically ends up being total DALYs averted and by sex
health_plot <- dplyr::bind_rows(incr_health, incr_health_by_sex) %>% 
  select(strategy, incr_dm_dalys, n_hef_eligible, threshold, sex) %>% 
  mutate(sex = ifelse(is.na(sex), "Total", sex))
health_plot$sex <- ordered(health_plot$sex,
                           levels = c("Male", "Female", "Total"))
health_plot$threshold <- factor(health_plot$threshold,
                                labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

ggplot(health_plot %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = sex, y = -incr_dm_dalys, fill = strategy)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  scale_fill_npg() +
  # scale_fill_manual(values = wes_palette("Darjeeling1", n = 3)) + 
  facet_grid(threshold ~ strategy) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.position = "none") + 
  labs(x = "Sex", y = "Diabetes-related DALYs averted") + 
  scale_y_continuous(labels = scales::comma)
ggsave("figures/dalys_by_sex.png")

# 4. Results: Individual expenditures ------------------------------------
# Summary of OOP ----------------------------------------------------------
incr_oop <- results %>% 
  group_by(strategy, threshold, coverage) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>%
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    n_hef_eligible = ifelse(threshold == "20%", sum(population20$income < hef_threshold20) * inf.fct, sum(population30$income < hef_threshold30) * inf.fct)
  )

# Bar plot by strategy
incr_oop$strategy <- ordered(incr_oop$strategy,
                             levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                             labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

ggplot(data = incr_oop %>% filter(strategy != "Current standard"), aes(x = strategy, y = oop_total / n_hef_eligible)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Type of OOP spending", y = "Amount per person eligible for HEF (2019 USD)", fill = "Strategy") 

# OOP by income quintile --------------------------------------------------
incr_oop_by_income <- results %>% 
  group_by(strategy, threshold, coverage, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  group_by(income_quintile) %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base']
  )

incr_oop_by_income$n_hef_eligible <- ifelse(incr_oop_by_income$threshold == "20%" & incr_oop_by_income$income_quintile == 1, n_hef_eligible_income20[1, 2][[1]], 
                                               ifelse(incr_oop_by_income$threshold == "30%" & incr_oop_by_income$income_quintile == 1, n_hef_eligible_income30[1, 2][[1]], n_hef_eligible_income30[2, 2][[1]]))

# Plot
incr_oop_by_income$strategy <- ordered(incr_oop_by_income$strategy,
                                       levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                       labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_oop_by_income$income_quintile <- ordered(incr_oop_by_income$income_quintile,
                                              levels = c(1:5))
incr_oop_by_income$threshold <- ordered(incr_oop_by_income$threshold, 
                                        levels = c("20%", "30%"),
                                        labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_oop_by_income$coverage <- ordered(incr_oop_by_income$coverage, levels = c("80%", "100%"),
                                       labels = c("HEF coverage: 80%", "HEF coverage: 100%"))

# OOP spending averted per HEF eligible
ggplot(data = incr_oop_by_income %>% filter(strategy != "Current standard"), 
       aes(x = strategy, y = -incr_oop_total, fill = income_quintile)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  # scale_fill_manual(values = cbp1, drop = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "OOP expenditures averted (USD 2019)", fill = "Income quintile")
ggsave("figures/oop_by_income.png")


# OOP by sex --------------------------------------------------------------
incr_oop_by_sex <- results %>% 
  group_by(strategy, threshold, coverage, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  group_by(sex) %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base']
  )

incr_oop_by_sex$n_hef_eligible <- ifelse(incr_oop_by_sex$threshold == "20%" & incr_oop_by_sex$sex == "Female", n_hef_eligible_sex20[1, 2][[1]], 
                                         ifelse(incr_oop_by_sex$threshold == "20%" & incr_oop_by_sex$sex == "Male", n_hef_eligible_sex20[2, 2][[1]],
                                                ifelse(incr_oop_by_sex$threshold == "30%" & incr_oop_by_sex$sex == "Female", n_hef_eligible_sex30[1, 2][[1]], n_hef_eligible_sex30[2, 2][[1]])))

# Plot
incr_oop_by_sex$strategy <- ordered(incr_oop_by_sex$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_oop_by_sex$threshold <- factor(incr_oop_by_sex$threshold, 
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_oop_by_sex$coverage <- ordered(incr_oop_by_sex$coverage, 
                                    levels = c("80%", "100%"),
                                    labels = c("HEF coverage: 80%", "HEF coverage: 100%"))

#### *************************************************************************************************************************
# OOP spending averted per HEF eligible by sex
ggplot(data = incr_oop_by_sex %>% filter(strategy != "Current standard"), 
       aes(x = strategy, y = -incr_oop_total / n_hef_eligible, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = wes_palette("Royal2", n = 2)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "OOP expenditures averted per person eligible for HEF (USD 2019)", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)
ggsave("figures/oop_by_sex.png")


# 4. Results: Financial risk protection ----------------------------------

incr_frp <- results %>% 
  group_by(strategy, threshold, coverage) %>% 
  summarise(
    che10 = sum(che10_result) * inf.fct,
    che25 = sum(che25_result) * inf.fct,
    che40 = sum(che40_result) * inf.fct,
    che10_disposable = sum(che10_disposable_result, na.rm = T) * inf.fct,
    che25_disposable = sum(che25_disposable_result, na.rm = T) * inf.fct,
    che40_disposable = sum(che40_disposable_result, na.rm = T) * inf.fct,
    pov = sum(pov_result, na.rm = T) * inf.fct,
    pov_disposable = sum(pov_disposable_result, na.rm = T)  * inf.fct
  ) %>% 
  ungroup %>% 
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
    n_hef_eligible = ifelse(threshold == "20%", sum(population20$income < hef_threshold20) * inf.fct, sum(population30$income < hef_threshold30) * inf.fct)
  )

incr_frp$strategy <- ordered(incr_frp$strategy,
                             levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                             labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_frp$threshold <- factor(incr_frp$threshold, 
                             levels = c("20%", "30%"),
                             labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))


# FRP by income quintile --------------------------------------------------
incr_frp_income <- results %>% 
  group_by(strategy, threshold, coverage, income_quintile) %>% 
  summarise(
    che10 = sum(che10_result) * inf.fct,
    che25 = sum(che25_result) * inf.fct,
    che40 = sum(che40_result) * inf.fct,
    che10_disposable = sum(che10_disposable_result, na.rm = T) * inf.fct,
    che25_disposable = sum(che25_disposable_result, na.rm = T) * inf.fct,
    che40_disposable = sum(che40_disposable_result, na.rm = T) * inf.fct,
    pov = sum(pov_result, na.rm = T) * inf.fct,
    pov_disposable = sum(pov_disposable_result, na.rm = T)  * inf.fct
  ) %>% 
  group_by(income_quintile) %>% 
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base']
  )

incr_frp_income$n_hef_eligible <- ifelse(incr_frp_income$threshold == "20%" & incr_frp_income$income_quintile == 1, n_hef_eligible_income20[1, 2][[1]], 
                                               ifelse(incr_frp_income$threshold == "30%" & incr_frp_income$income_quintile == 1, n_hef_eligible_income30[1, 2][[1]], n_hef_eligible_income30[2, 2][[1]]))

# Plot
incr_frp_income$strategy <- ordered(incr_frp_income$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_frp_income$income_quintile <- ordered(incr_frp_income$income_quintile,
                                           levels = c(1:5))
incr_frp_income$threshold <- ordered(incr_frp_income$threshold, 
                                     levels = c("20%", "30%"),
                                     labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_frp_income$coverage <- ordered(incr_frp_income$coverage, 
                                    levels = c("80%", "100%"),
                                    labels = c("HEF coverage: 80% of direct medical costs", "HEF coverage: 100% of direct medical costs"))
# incr_frp_income$frp <- ordered(incr_frp_income$frp, 
#                                levels = c("incr_che10", "incr_che25", "incr_che40", "incr_che10_disposable", "incr_che25_disposable", "incr_che40_disposable", "incr_pov", "incr_pov_disposable"),
#                                labels = c("CHE threshold: 10%", "CHE threshold: 25%", "CHE threshold: 40%", 
#                                           "CHE threshold: 10% \n(disposable income)", "CHE threshold: 25% \n(disposable income)", "CHE threshold: 40% \n(disposable income)",
#                                           "Impoverishment", "Impoverishment \n(disposable income)"))

#### *************************************************************************************************************************
# CHE averted
ggplot(data = incr_frp_income %>% filter(strategy != "Current standard" & income_quintile == 1), 
       aes(x = strategy, y = -incr_che40)) +
  geom_bar(stat = "identity") + 
  facet_grid(threshold~coverage) +
  theme_bw() + 
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Cases of CHE averted per person eligible for HEF") + 
  scale_y_continuous(labels = scales::comma) 
ggsave("figures/che.png")

# FRP by sex --------------------------------------------------------------
incr_frp_by_sex <- results %>% 
  group_by(strategy, threshold, coverage, sex) %>% 
  summarise(
    che10 = sum(che10_result) * inf.fct,
    che25 = sum(che25_result) * inf.fct,
    che40 = sum(che40_result) * inf.fct,
    che10_disposable = sum(che10_disposable_result, na.rm = T) * inf.fct,
    che25_disposable = sum(che25_disposable_result, na.rm = T) * inf.fct,
    che40_disposable = sum(che40_disposable_result, na.rm = T) * inf.fct,
    pov = sum(pov_result, na.rm = T) * inf.fct,
    pov_disposable = sum(pov_disposable_result, na.rm = T)  * inf.fct
  ) %>% 
  group_by(sex) %>% 
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base']
  )

incr_frp_by_sex$n_hef_eligible <- ifelse(incr_frp_by_sex$threshold == "20%" & incr_frp_by_sex$sex == "Female", n_hef_eligible_sex20[1, 2][[1]], 
                                         ifelse(incr_frp_by_sex$threshold == "20%" & incr_frp_by_sex$sex == "Male", n_hef_eligible_sex20[2, 2][[1]],
                                                ifelse(incr_frp_by_sex$threshold == "30%" & incr_frp_by_sex$sex == "Female", n_hef_eligible_sex30[1, 2][[1]], n_hef_eligible_sex30[2, 2][[1]])))


# Plot
incr_frp_by_sex$strategy <- ordered(incr_frp_by_sex$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr_frp_by_sex$threshold <- factor(incr_frp_by_sex$threshold, 
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_frp_by_sex$coverage <- ordered(incr_frp_by_sex$coverage, 
                                    levels = c("80%", "100%"),
                                    labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
# incr_frp_by_sex$frp <- ordered(incr_frp_by_sex$frp, 
#                                  levels = c("incr_che10", "incr_che25", "incr_che40", "incr_che10_disposable", "incr_che25_disposable", "incr_che40_disposable", "incr_pov", "incr_pov_disposable"),
#                                  labels = c("CHE threshold: 10%", "CHE threshold: 25%", "CHE threshold: 40%", 
#                                             "CHE threshold: 10% \n(disposable income)", "CHE threshold: 25% \n(disposable income)", "CHE threshold: 40% \n(disposable income)",
#                                             "Impoverishment", "Impoverishment \n(disposable income)"))

#### *************************************************************************************************************************
# CHE averted by sex
ggplot(data = incr_frp_by_sex %>% filter(strategy != "Current standard"), 
       aes(x = strategy, y = -incr_che40, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = wes_palette("Royal2", n = 2)) + 
  # scale_fill_manual(values = cbp1, drop = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Cases of CHE averted", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)
ggsave("figures/che_sex.png")

# 5. ICERs ---------------------------------------------------------------

## Health gains
# DALYs
icer.dalys <- cbind(incr_costs, incr_health[4:15]) %>% 
  mutate(
    icer_dalys = incr_cost_total / -incr_dm_dalys
  )
icer.dalys$threshold <- factor(icer.dalys$threshold, 
                               levels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"),
                               labels = c("20%", "30%"))
icer.dalys$coverage <- ordered(icer.dalys$coverage,
                               levels = c("HEF coverage: 80%", "HEF coverage: 100%"),
                               labels = c("80%", "100%"))
icer.dalys.table <- icer.dalys %>% select(strategy, threshold, coverage, cost_total, incr_cost_total, dm_dalys, incr_dm_dalys, icer_dalys)
icer.dalys.table[,4:8] <- format(icer.dalys.table[,4:8], big.mark = ",")
print(xtable::xtable(icer.dalys.table, digits = 0), include.rownames = FALSE)

### ICER PLOT 
# DALYs
icer.dalys$strategy <- ordered(icer.dalys$strategy,
                               levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                               labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
icer.dalys$threshold <- factor(icer.dalys$threshold, 
                               levels = c("20%", "30%"),
                               labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
ggplot(icer.dalys %>% filter(coverage == "80%") %>% filter(!is.na(strategy)), 
       aes(x = -incr_dm_dalys, y = incr_cost_total / n_hef_eligible)) + 
  geom_point(aes(color = strategy), size = 4) + 
  facet_grid(threshold ~ .) + 
  theme_bw() + 
  scale_color_npg() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_shape_discrete(guide = FALSE) + 
  labs(x = "Diabetes-related DALYs averted", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Strategy") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0))
ggsave("figures/icer_dalys.png")

# By sex
icer.dalys <- cbind(incr_costs, incr_health[4:15]) %>% 
  mutate(
    icer_dalys = incr_cost_total / -incr_dm_dalys
  )

icer.dalys.sex <- merge(incr_costs_by_sex, incr_health_by_sex, by = c("strategy", "sex", "threshold", "coverage", "n_hef_eligible")) %>% 
  mutate(
    icer_dalys_sex = incr_cost_total / -incr_dm_dalys
  )
icer.dalys.sex$threshold <- ordered(icer.dalys.sex$threshold, 
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

icer.dalys.sex.plot <- dplyr::bind_rows(icer.dalys, icer.dalys.sex) %>% 
  select(strategy, sex, threshold, incr_cost_total, incr_dm_dalys, n_hef_eligible) %>% 
  mutate(
    sex = ifelse(is.na(sex), "Total", sex)
  )

ggplot(icer.dalys.sex.plot %>% filter(!is.na(strategy)), aes(x = -incr_dm_dalys, y = incr_cost_total / n_hef_eligible)) + 
  geom_point(aes(color = strategy), size = 4) + 
  facet_grid(sex ~ threshold) + 
  theme_bw() + 
  scale_color_npg() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_shape_discrete(guide = FALSE) + 
  labs(x = "Diabetes-related DALYs averted", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Strategy") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0))
ggsave("figures/icer_dalys.png")

# # OOP
incr_oop$threshold <- ordered(incr_oop$threshold, 
                              levels = c("20%", "30%"),
                              labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_oop$coverage <- incr_oop$coverage <- ordered(incr_oop$coverage, 
                                                  levels = c("80%", "100%"),
                                                  labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
icer.oop <- merge(incr_costs, incr_oop, by = c("strategy", "threshold", "coverage", "n_hef_eligible"))

icer.oop.sex <- merge(incr_costs_by_sex, incr_oop_by_sex, by = c("strategy", "sex", "threshold", "coverage", "n_hef_eligible")) %>% 
  mutate(
    icer_oop_sex = incr_cost_total / incr_oop_total
  ) 

icer.oop.plot <- dplyr::bind_rows(icer.oop, icer.oop.sex) %>% 
  mutate(
    sex = ifelse(is.na(sex), "Total", sex)
  )  
icer.oop.plot$coverage <- icer.oop.plot$coverage <- ordered(icer.oop.plot$coverage, 
                                                  levels = c("HEF coverage: 100%", "HEF coverage: 80%"),
                                                  labels = c("100%", "80%"))

ggplot(icer.oop.plot, aes(x = -incr_oop_total / 1e6, y = incr_cost_total / n_hef_eligible)) +
  geom_point(aes(color = strategy, shape = coverage), size = 4) +
  facet_grid(sex ~ threshold) +
  theme_bw() +
  scale_color_npg(drop = FALSE) +
  scale_shape_manual(values = c(16, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "OOP expenditures averted (millions, USD 2019)", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Coverage strategy", shape = "HEF coverage") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0))
ggsave("figures/icer_oop.png")

## Financial risk protection
incr_frp$coverage <- ordered(incr_frp$coverage,
                             levels = c("80%", "100%"),
                             labels = c("HEF coverage: 80%", "HEF coverage: 100%"))

icer.che <- merge(incr_costs, incr_frp, by = c("strategy", "threshold", "coverage", "n_hef_eligible")) %>% 
  mutate(
    icer_che = incr_cost_total / -incr_che40
  )

incr_costs_by_sex$threshold <- ordered(incr_costs_by_sex$threshold, 
                                       levels = c("20%", "30%"),
                                       labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr_costs_by_sex$coverage <- ordered(incr_costs_by_sex$coverage,
                                      levels = c("80%", "100%"),
                                      labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
icer.frp.sex <-  merge(incr_costs_by_sex, incr_frp_by_sex, by = c("strategy", "sex", "threshold", "coverage", "n_hef_eligible")) %>% 
  mutate(
    icer.che.sex = incr_cost_total / -incr_che40
  )

icer.frp.plot <- dplyr::bind_rows(icer.che, icer.frp.sex) %>% 
  mutate(
    sex = ifelse(is.na(sex), "Total", sex)
  )  

icer.frp.plot$coverage <- ordered(icer.frp.plot$coverage,
                                  levels = c("HEF coverage: 100%", "HEF coverage: 80%"),
                                  labels = c("100%", "80%"))

ggplot(icer.frp.plot %>%filter(!is.na(strategy)), aes(x = -incr_che40 / 1e3, y = incr_cost_total / n_hef_eligible)) + 
  geom_point(aes(color = strategy, shape = coverage), size = 4) + 
  facet_grid(sex ~ threshold) + 
  theme_bw() + 
  scale_color_npg(drop = FALSE) + 
  scale_shape_manual(values = c(16, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Cases of CHE averted (thousands)", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Strategy", shape = "HEF Coverage") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0))
ggsave("figures/icer_frp.png")

ggplot(icer.frp.plot %>%filter(!is.na(strategy)) %>% filter(threshold == "Poorest 30% eligible for HEF"), aes(x = -incr_pov / 1e3, y = incr_cost_total / n_hef_eligible)) + 
  geom_point(aes(color = strategy, shape = coverage), size = 4) + 
  facet_grid(sex ~ .) + 
  theme_bw() + 
  scale_color_npg(drop = FALSE) + 
  scale_shape_manual(values = c(16, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Cases of poverty averted (thousands)", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Strategy", shape = "HEF Coverage") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0))
ggsave("figures/icer_pov.png")

# Print table
icer.che$threshold <- factor(icer.che$threshold, 
                               levels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"),
                               labels = c("20%", "30%"))
icer.che$coverage <- ordered(icer.che$coverage,
                               levels = c("HEF coverage: 80%", "HEF coverage: 100%"),
                               labels = c("80%", "100%"))
icer.che.table <- icer.che %>% select(strategy, threshold, coverage, cost_total, incr_cost_total, che40, incr_che40, icer_che)
icer.che.table[,4:8] <- format(icer.che.table[,4:8], big.mark = ",")
print(xtable::xtable(icer.che.table, digits = 0), include.rownames = FALSE)


# Inspect model runs ------------------------------------------------------

checks <- results20_80 %>% 
  group_by(id, strategy) %>% 
  arrange(cycle) %>% 
  summarise(
    n_cycles = last(cycle),
    last_age = last(age_start + cycle),
    has_stopped_at_70 = last_age == 70,
    stopped_reason = ifelse(has_stopped_at_70, 'hit_70', last(state))
  )

checks %>% ggplot(data = ., aes(x = n_cycles, fill = stopped_reason)) +
    geom_bar() + facet_grid(strategy ~ .)

checks %>% ggplot(data = ., aes(x = last_age, fill = stopped_reason)) +
  geom_bar() + facet_grid(strategy ~ .)


