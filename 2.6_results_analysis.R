##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 2.6 Analysis of results
#### Copyright Isabelle Feldhaus
#### 24 January 2020
##################################


# 0. Preparation ----------------------------------------------------------

library(wesanderson)
library(ggsci)
library(xtable)

## Designating appropriate date types in results
results_threshold20_coverage80 <- readr::read_csv('output/results_60000237L_2020-01-23_1849.csv')
results_threshold30_coverage80 <- readr::read_csv('output/results_60002004L_2020-01-24_1156.csv')
population20 <- readr::read_csv("output/ichar_500000L_2020-01-23_1757.csv")
population30 <- readr::read_csv("output/ichar_500000L_2020-01-24_1059.csv")

## Define income quintile thresholds and add variable to results dataframe 
# avg_income <- 963.86 * 12
# income_quintile_thresholds <- qgamma(c(0.2, 0.4, 0.6, 0.8), shape = 0.5, scale = avg_income)
income_quintile_thresholds <- quantile(income_dist, c(0.2, 0.4, 0.6, 0.8))
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
hef_threshold20 <- quantile(income_dist, 0.2)[[1]]
hef_threshold30 <- quantile(income_dist, 0.3)[[1]]

## Add sex and income quintile to results data frame
results20_80 <- merge(results_threshold20_coverage80, population20 %>% select(id, income_quintile, sex, hef), by = "id")
results30_80 <- merge(results_threshold30_coverage80, population30 %>% select(id, income_quintile, sex, hef), by = "id")

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
         pov_disposable_result = pov(oop_total, disposable_income))

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
         pov_disposable_result = pov(oop_total, disposable_income))

# Inflation factor
inf.fct <- 16e6 / nrow(population20)

## Quick look
# Prevalence of diabetes over 45 years
n_diabetes <- results20_80 %>% filter(state != "healthy") %>% filter(state != "other_death") %>% group_by(id) %>% summarise(n = n()) %>% nrow() * inf.fct
n_diabetes_sex <- results20_80 %>% filter(state != "healthy") %>% filter(state != "other_death") %>% group_by(id, sex) %>% summarise(n = n()) %>% 
  ungroup %>% group_by(sex) %>% summarise(n_diabetes = n() * inf.fct)
prevalence <- (n_diabetes / inf.fct) / nrow(population20)

# Number of HEF beneficiaries
# Threshold 20, Coverage 80/100 
n_hef_beneficiaries20_80 <- results20_80 %>% 
  group_by(strategy, id) %>%
  summarise(n_cycles = sum(hef)) %>% 
  filter(n_cycles != 0) %>% 
  ungroup %>% 
  group_by(strategy) %>% 
  summarise(n_hef = n())
n_hef_beneficiaries20_80 <- n_hef_beneficiaries20_80[1, 2][[1]] 
n_hef_beneficiaries20_80 <- n_hef_beneficiaries20_80 * inf.fct

# Threshold 30, Coverage 80/100
n_hef_beneficiaries30_80 <- results30_80 %>% 
  group_by(strategy, id) %>%
  summarise(n_cycles = sum(hef)) %>% 
  filter(n_cycles != 0) %>% 
  ungroup %>% 
  group_by(strategy) %>% 
  summarise(n_hef = n())
n_hef_beneficiaries30_80 <- n_hef_beneficiaries30_80[1, 2][[1]] 
n_hef_beneficiaries30_80 <- n_hef_beneficiaries30_80 * inf.fct

# Count of HEF benficiaries by sex
n_hef_beneficiaries20_80.sex <- results20_80 %>% 
  filter(strategy == "base") %>% 
  group_by(id, sex) %>%
  summarise(n_cycles = sum(hef)) %>% 
  filter(n_cycles != 0) %>% 
  ungroup %>% 
  group_by(sex) %>% 
  summarise(n_hef = n() * inf.fct)

n_hef_beneficiaries20_100.sex <- results20_100 %>% 
  filter(strategy == "base") %>% 
  group_by(id, sex) %>%
  summarise(n_cycles = sum(hef)) %>% 
  filter(n_cycles != 0) %>% 
  ungroup %>% 
  group_by(sex) %>% 
  summarise(n_hef = n() * inf.fct)

# Filter results to HEF eligible individuals 
results20_80 <- results20_80 %>% filter(income < hef_threshold20)
results20_100 <- results20_100 %>% filter(income < hef_threshold20)
results30_80 <- results30_80 %>% filter(income < hef_threshold30)
results30_100 <- results30_100 %>% filter(income < hef_threshold30)

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
statesByStrategy <- statesByStrategy %>% mutate(n_diabetes1 = ifelse(sex == "Female", n_diabetes_sex$n_diabetes[n_diabetes_sex$sex == "Female"], n_diabetes_sex$n_diabetes[n_diabetes_sex$sex == "Male"]))

# Look at first 2 income quintiles
statesByStrategy %>% select(strategy, state, income_quintile, n) %>% 
  filter(state != "Healthy" & state != "Other death") %>%  
  filter(income_quintile %in% c(1, 2)) %>% spread(strategy, n)

# Plot number of life-years spent in each state by strategy per diabetic individual in population
statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
ggplot(aes(x = state, y = n / n_diabetes, fill = strategy)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Strategy")
# ggsave("figures/state_by_strategy.png")

statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
  # filter(income_quintile %in% c(1, 2)) %>%
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = income_quintile)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(strategy ~ .) +
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Income quintile")

statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
  # filter(income_quintile %in% c(1, 2)) %>% 
  filter(strategy %in% c("Current standard", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = strategy)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  # facet_wrap(income_quintile ~ .) +
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Income quintile")

statesByStrategy %>%
  filter(state != "Healthy" & state != "Other death") %>%  
  # filter(income_quintile %in% c(1, 2)) %>%
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = sex)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(strategy ~ .) +
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Sex")

statesByStrategy %>%
  filter(!(state %in% c("Healthy", "Other death"))) %>%  
  # filter(income_quintile %in% c(1, 2)) %>% 
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = sex, y = n / n_diabetes1, fill = strategy)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(state ~ .) +
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Strategy") 

statesByStrategy %>%
  filter(state %in% c("No treatment", "Insulin therapy")) %>%  
  # filter(income_quintile %in% c(1, 2)) %>% 
  filter(strategy %in% c("Current standard", "Diagnosis only", "Treatment only")) %>% 
ggplot(aes(x = state, y = n / n_diabetes, fill = strategy)) +
# ggplot(aes(x = state, y = ratio_to_base, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(sex ~ .) +
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "State", y = "Average life-years per diabetic patient", fill = "Strategy")

# results20_80 %>% 
#   filter(!(state %in% c("healthy", "other_death"))) %>% order(c('id', 'cycle')) %>% 
#   group_by(strategy, id) %>% 
#   summarise(duration = last(cycle),
#             undiagnosed = sum(state == "diabetes") / duration,
#             notx = sum(state == "notx"),
#             oad = sum(state == "oad"),
#             ins = sum(state == "ins"), 
#             nephro = sum(state == "nephro"),
#             retino = sum(state == "retino"),
#             neuro = sum(state == "neuro"),
#             ang = sum(state == "ang"), 
#             pvd = sum(state == "pvd"),
#             mi = sum(state == "mi"),
#             stroke = sum(state == "stroke"),
#             fail = sum(state == "fail"))


# 2. Results: Costs --------------------------------------------------------
# Summary of costs --------------------------------------------------------
# threshold 20, coverage 80
results.costs20_80 <- results20_80 %>% group_by(strategy) %>% summarise(
  cost_diagnosis = sum(cost_diagnosis) * inf.fct,
  cost_tx = sum(cost_tx) * inf.fct,
  cost_complications = sum(cost_complications) * inf.fct, 
  cost_total = sum(cost_total) * inf.fct, 
  threshold = "20%", 
  coverage = "80%"
)
results.costs.plot20_80 <- results.costs20_80 %>% gather(key = "cost", value = "value", -c(strategy, threshold, coverage))

# threshold 30, coverage 80
results.costs30_80 <- results30_80 %>% group_by(strategy) %>% summarise(
  cost_diagnosis = sum(cost_diagnosis) * inf.fct,
  cost_tx = sum(cost_tx) * inf.fct,
  cost_complications = sum(cost_complications) * inf.fct, 
  cost_total = sum(cost_total) * inf.fct, 
  threshold = "30%", 
  coverage = "80%"
)
results.costs.plot30_80 <- results.costs30_80 %>% gather(key = "cost", value = "value", -c(strategy, threshold, coverage))

results.costs.plot <- rbind(results.costs.plot20_80, results.costs.plot30_80)

## Plot
results.costs.plot$cost <- ordered(results.costs.plot$cost,
                                   levels = c("cost_diagnosis", "cost_tx", "cost_complications", "cost_total"),
                                   labels = c("Screening + \nLaboratory", "Medications", "Service Utilization \nfor Complications", "Total"))
results.costs.plot$strategy <- ordered(results.costs.plot$strategy,
                                       levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                       labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
results.costs.plot$threshold <- factor(results.costs.plot$threshold,
                                       labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

# Stacked bar plot of types of costs per capita by strategy
ggplot(data = results.costs.plot %>% filter(cost != "Total"), aes(x = strategy, y = value / nrow(population), fill = cost)) +
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() + 
  facet_wrap(.~threshold) +
  # scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Cost per capita (2019 USD)", fill = "Type of cost") + 
  scale_y_continuous(labels = scales::comma) + 
  coord_cartesian(xlim = NULL, ylim = c(1e4, 2.3e4))


# Incremental costs -------------------------------------------------------
incr.costs20_80 <- results20_80 %>% 
  group_by(strategy) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "20%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.costs30_80 <- results30_80 %>% 
  group_by(strategy) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "30%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.costs <- rbind(incr.costs20_80, incr.costs30_80)
incr.costs.plot <- incr.costs %>% gather(key = "cost", value = "value", -c(strategy, threshold, n_hef_eligible))

incr.costs.plot$cost <- ordered(incr.costs.plot$cost,
                                levels = c("incr_cost_diagnosis", "incr_cost_tx", "incr_cost_complications", "incr_cost_total"),
                                labels = c("Screening + \nLaboratory", "Medications", "Service Utilization \nfor Complications", "Total"))
incr.costs.plot$strategy <- ordered(incr.costs.plot$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.costs.plot$threshold <- factor(incr.costs.plot$threshold,
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

#### *************************************************************************************************************************
# Incremental costs per HEF eligible
ggplot(data = incr.costs.plot %>% filter(cost != "Total") %>% filter(strategy != "Current standard"), 
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
ggsave("figures/incr_costs_per_hef.png")

### Summary table
summary.costs <- rbind(incr.costs20_80, incr.costs30_80)
total.costs <- summary.costs %>% select(strategy, threshold, cost_total, incr_cost_total)

# Costs by income quintile ------------------------------------------------
incr.costs.income.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    cost_total = cost_total - cost_total[strategy == 'base']
  )
incr.costs.income1.20_80$threshold <- "20%"
incr.costs.income1.20_80$coverage <- "80%"

# Bar plot by strategy
results.costs.income.plot <- results.costs.income %>% gather(key = "cost", value = "value", -c(strategy, income_quintile))
results.costs.income.plot$cost <- ordered(results.costs.income.plot$cost,
                                          levels = c("cost_diagnosis", "cost_tx", "cost_complications", "cost_total"),
                                          labels = c("Diagnosis", "Treatment", "Complications", "Total"))
results.costs.income.plot$strategy <- ordered(results.costs.income.plot$strategy,
                                              levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                       labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
results.costs.income.plot$income_quintile <- ordered(results.costs.income.plot$income_quintile,
                                                     levels = c(1:5))

# Cost per capita by income 
ggplot(data = results.costs.income.plot %>% filter(cost == "Total"),
       aes(x = income_quintile, y = value / nrow(population), fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  # scale_fill_manual(values = cbp1) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Income quintile", y = "Cost per capita (2019 USD)", fill = "Strategy") + 
  scale_y_continuous(labels = scales::comma)

# # Cost per capita of quintile by income -- unfinished, same plot as above
# ggplot(data = results.costs.income.plot %>% filter(cost == "Total"),
#        aes(x = income_quintile, y = value / nrow(population), fill = strategy)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_minimal() +
#   # scale_fill_manual(values = cbp1) +
#   scale_fill_brewer(palette = "Paired") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14)) +
#   labs(x = "Income quintile", y = "Cost per capita (2019 USD)", fill = "Strategy") + 
#   scale_y_continuous(labels = scales::comma)
# 
# # Cost per HEF beneficiary within quintiles by income -- unfinished, same plot as above
# ggplot(data = results.costs.income.plot %>% filter(cost == "Total"),
#        aes(x = income_quintile, y = value / nrow(population), fill = strategy)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_minimal() +
#   # scale_fill_manual(values = cbp1) +
#   scale_fill_brewer(palette = "Paired") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14)) +
#   labs(x = "Income quintile", y = "Cost per capita (2019 USD)", fill = "Strategy") + 
#   scale_y_continuous(labels = scales::comma)


# Costs by sex ------------------------------------------------------------
incr.costs.male.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%"
  )

incr.costs.female.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%"
  )

incr.costs.male.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%"
  )

incr.costs.female.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%"
  )

incr.costs.male.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%"
  )

incr.costs.female.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%"
  )

incr.costs.male.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%"
  )

incr.costs.female.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    cost_diagnosis = sum(cost_diagnosis) * inf.fct,
    cost_tx = sum(cost_tx) * inf.fct,
    cost_complications = sum(cost_complications) * inf.fct, 
    cost_total = sum(cost_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
    incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
    incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
    incr_cost_total = cost_total - cost_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%"
  )


incr.costs.sex <- rbind(incr.costs.male.20_80, incr.costs.female.20_80, incr.costs.male.20_100, incr.costs.female.20_100,
                        incr.costs.male.30_80, incr.costs.female.30_80, incr.costs.male.30_100, incr.costs.female.30_100)

# Bar plot by strategy
results.costs.sex.plot <- results.costs.sex %>% gather(key = "cost", value = "value", -c(strategy, sex))
results.costs.sex.plot$cost <- ordered(results.costs.sex.plot$cost,
                                       levels = c("cost_diagnosis", "cost_tx", "cost_complications", "cost_total"),
                                       labels = c("Diagnosis", "Treatment", "Complications", "Total"))
results.costs.sex.plot$strategy <- ordered(results.costs.sex.plot$strategy,
                                           levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                           labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

ggplot(data = results.costs.sex.plot %>% filter(cost == "Total"),
       aes(x = strategy, y = value / nrow(population), fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Strategy", y = "Cost per capita (2019 USD)", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)



# 3. Results: Health gains ------------------------------------------------
# Summary of health gains -------------------------------------------------
incr.health20_80 <- results20_80 %>% 
  group_by(strategy) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health30_80 <- results30_80 %>% 
  group_by(strategy) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "30%"
  )

incr.health <- rbind(incr.health20_80, incr.health30_80)
dalys.health <- incr.health %>% select(strategy, threshold, dm_dalys, incr_dm_dalys) # DALYs incurred (not averted) as denominator

# Bar plot by strategy
results.health.plot <- incr.health20_80 %>% select(strategy, dm_dalys) %>% gather(key = "dalys", value = "value", -strategy)
results.health.plot$strategy <- ordered(results.health.plot$strategy,
                                        levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                       labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
results.health.plot$value <- results.health.plot$value * -1
ggplot(data = results.health.plot %>% filter(!(strategy %in% c("Current standard", "Complications only"))), aes(x = dalys, y = value, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_minimal() + 
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = NULL,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = NULL, y = "DALYs averted", fill = "Coverage strategy")

incr.health20_80.state <- results20_80 %>% 
  group_by(strategy, state, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base' & state == state & sex == sex],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'& state == state & sex == sex],
    threshold = "20%"
  )

incr.health30_80.state <- results30_80 %>% 
  group_by(strategy, state, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base' & state == state & sex == sex],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base' & state == state & sex == sex], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base' & state == state & sex == sex], 
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'& state == state & sex == sex],
    threshold = "30%"
  )

health.state <- rbind(incr.health20_80.state, incr.health30_80.state)
#### *************************************************************************************************************
# DALYs by state and sex
ggplot(data = health.state, 
       aes(x = strategy, y = -incr_dm_dalys, fill = state)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() +
  facet_wrap(.~sex) +
  # scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") + 
  scale_y_continuous(labels = scales::comma) 

health.state %>% 
  mutate(dalys_averted = ifelse(incr_dm_dalys < 0, -incr_dm_dalys, 0)) %>% 
  group_by(strategy, sex) %>% 
  mutate(dalys_averted_sum = sum(dalys_averted)) %>% 
  filter(state == "neuro" & sex == "Male") %>% 
  mutate(prop_dalys_averted = dalys_averted / dalys_averted_sum) %>% 
  select(strategy, state, sex, dalys_averted, dalys_averted_sum, prop_dalys_averted)


# Health gains by income quintile -----------------------------------------
incr.health.income1.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.income2.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.income3.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 3) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.income4.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 4) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.income5.20_80 <- results20_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 5) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.income1.30_80 <- results30_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "30%"
  )

incr.health.income2.30_80 <- results30_80 %>% 
  group_by(strategy, income_quintile) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "30%"
  )

incr.health.income <- rbind(incr.health.income1.20_80, incr.health.income2.20_80,
                            incr.health.income1.30_80, incr.health.income2.30_80)
incr.health.income$coverage <- NA

# Bar plot of total health effects by strategy
# results.health.income.plot <- results.health.income %>% gather(key = "type", value = "value", -c(strategy, income_quintile))
# results.health.income.plot$strategy <- ordered(results.health.income.plot$strategy,
#                                                levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
#                                                labels = c("Current standard", "Diagnosis only", "Treatment only", "Diagnosis + treatment", "Complications only", "Treatment + complications", "Diagnosis + treatment + complications"))
# results.health.income.plot$income_quintile <- ordered(results.health.income.plot$income_quintile,
#                                                       levels = c(1:5))
# 
# ggplot(data = results.health.income.plot %>% filter(type == "dm_deaths"), 
#        aes(x = income_quintile, y = value, fill = strategy)) +
#   geom_bar(stat = "identity", position = position_dodge()) + 
#   theme_minimal() + 
#   scale_fill_manual(values = cbp1) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14)) + 
#   labs(x = "Income quintile", y = "Diabetes-related deaths", fill = "Strategy") + 
#   scale_y_continuous(labels = scales::comma) + 
#   coord_cartesian(ylim=c(1.8e5,2.1e5))
# 
# ggplot(data = results.health.income.plot %>% filter(type == "dm_dalys"), 
#        aes(x = income_quintile, y = value, fill = strategy)) +
#   geom_bar(stat = "identity", position = position_dodge()) + 
#   theme_minimal() + 
#   scale_fill_manual(values = cbp1) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14)) + 
#   labs(x = "Income quintile", y = "Diabetes-related DALYs", fill = "Strategy") + 
#   scale_y_continuous(labels = scales::comma) + 
#   coord_cartesian(ylim=c(4e6,4.75e6))


# Plots of incremental health effects
incr.health.income.plot <- incr.health.income %>% 
  select(strategy, income_quintile, incr_dm_dalys, threshold) %>% 
  gather(key = "incr_dm_dalys", value = "value", -c(strategy, income_quintile, threshold))
incr.health.income.plot$strategy <- ordered(incr.health.income.plot$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.health.income.plot$income_quintile <- ordered(incr.health.income.plot$income_quintile,
                                                      levels = c(1:5))
incr.health.income.plot$value <- incr.health.income.plot$value * -1 # to convert to DALYs _averted_

ggplot(data = incr.health.income.plot %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = strategy, y = value, fill = income_quintile)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  facet_grid(.~threshold) +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") + 
  scale_y_continuous(labels = scales::comma) 

# Per capita
ggplot(data = incr.health.income.plot %>% filter(!(strategy %in% c("Current standard", "Complications only"))), 
       aes(x = strategy, y = value / nrow(population), fill = income_quintile)) +
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted per capita", fill = "Income quintile") + 
  scale_y_continuous() 

# Per diabetic patient
incr.health.income.plot %>% 
  filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")) %>% 
  filter(income_quintile %in% c(1, 2)) %>% 
ggplot(aes(x = strategy, y = value / n_diabetes, fill = income_quintile)) +
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "top") + 
  labs(x = "Coverage strategy", y = "DALYs averted per diabetes case", fill = "Income quintile") + 
  scale_y_continuous() 

# Disaggregating states
results.health.income.state <- results20_80 %>% group_by(strategy, income_quintile, state) %>% summarise(
  dm_deaths = sum(state == "dm_death") * inf.fct,
  other_deaths = sum(state == "other_death") * inf.fct,
  total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
  dm_dalys = sum(dalys_dm) * inf.fct,
  other_dalys = sum(dalys_other) * inf.fct,
  total_dalys = sum(dalys_dm, dalys_other) * inf.fct
)

incr.health.income1.state <- results20_80 %>% 
  group_by(strategy, income_quintile, state) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    other_dalys = other_dalys - other_dalys[strategy == 'base'],
    total_dalys = total_dalys - total_dalys[strategy == 'base']
  )

incr.health.income1.state.prop <- incr.health.income1.state %>% 
  select(strategy, state, dm_dalys) %>% 
  filter(dm_dalys < 0) %>% 
  group_by(strategy) %>% 
  mutate(prop = dm_dalys / sum(dm_dalys))
incr.dalys.state.prop <- incr.health.income1.state.prop %>% ungroup %>% select(strategy, state, prop) %>% spread(key = "state", value = "prop")

incr.health.income1.state.plot <- incr.health.income1.state %>% select(strategy, income_quintile, state, dm_dalys) %>% gather(key = "dm_dalys", value = "value", -c(strategy, income_quintile, state)) %>% select(strategy, income_quintile, state, value)
incr.health.income1.state.plot$strategy <- ordered(incr.health.income1.state.plot$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.health.income1.state.plot$income_quintile <- ordered(incr.health.income1.state.plot$income_quintile,
                                                      levels = c(1:5))
incr.health.income1.state.plot$value <- incr.health.income1.state.plot$value * -1 # to convert to DALYs _averted_

ggplot(data = incr.health.income1.state.plot %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = strategy, y = value, fill = state)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  # scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") + 
  scale_y_continuous(labels = scales::comma) 


# Health gains by sex -----------------------------------------------------
incr.health.male20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.female20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "20%"
  )

incr.health.male30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "30%"
  )

incr.health.female30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    incr_other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    incr_total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    incr_other_dalys = other_dalys - other_dalys[strategy == 'base'],
    incr_total_dalys = total_dalys - total_dalys[strategy == 'base'],
    threshold = "30%"
  )

incr.health.sex <- rbind(incr.health.male20_80, incr.health.female20_80, incr.health.male30_80, incr.health.female30_80)

# Plot
incr.health.sex.plot <- incr.health.sex %>% select(strategy, sex, dm_dalys) %>% gather(key = "dm_dalys", value = "value", -c(strategy, sex)) %>% select(strategy, sex, value)
incr.health.sex.plot$strategy <- ordered(incr.health.sex.plot$strategy,
                                         levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.health.sex.plot$value <- incr.health.sex.plot$value * -1 # to convert to DALYs _averted_

ggplot(data = incr.health.sex.plot %>% filter(!(strategy %in% c("Current standard", "Complications only"))), 
       aes(x = strategy, y = value, fill = sex)) +
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)

# Per diabetic patient
incr.health.sex.plot %>% 
  filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")) %>% 
ggplot(aes(x = strategy, y = value / n_diabetes, fill = sex)) +
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Royal2", n = 2)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "top") + 
  labs(x = "Coverage strategy", y = "DALYs averted per diabetes case", fill = "Sex") + 
  scale_y_continuous(labels = scales::comma)

incr.health.male.state <- results20_80 %>% 
  group_by(strategy, sex, state) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    other_dalys = other_dalys - other_dalys[strategy == 'base'],
    total_dalys = total_dalys - total_dalys[strategy == 'base']
  )

incr.health.female.state <- results20_80 %>% 
  group_by(strategy, sex, state) %>% 
  summarise(
    dm_deaths = sum(state == "dm_death") * inf.fct,
    other_deaths = sum(state == "other_death") * inf.fct,
    total_deaths = sum(state == "dm_death" | state == "other_death") * inf.fct, 
    dm_dalys = sum(dalys_dm) * inf.fct,
    other_dalys = sum(dalys_other) * inf.fct,
    total_dalys = sum(dalys_dm, dalys_other) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
    other_deaths = other_deaths - other_deaths[strategy == 'base'], 
    total_deaths = total_deaths - total_deaths[strategy == 'base'], 
    dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
    other_dalys = other_dalys - other_dalys[strategy == 'base'],
    total_dalys = total_dalys - total_dalys[strategy == 'base']
  )

incr.health.sex.state <- rbind(incr.health.male.state, incr.health.female.state)
incr.health.sex.state.plot <- incr.health.sex.state %>% select(strategy, sex, state, dm_dalys) %>% gather(key = "dm_dalys", value = "value", -c(strategy, sex, state)) %>% select(strategy, sex, state, value)
incr.health.sex.state.plot$strategy <- ordered(incr.health.sex.state.plot$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.health.sex.state.plot$value <- incr.health.sex.state.plot$value * -1 # to convert to DALYs _averted_

ggplot(data = incr.health.sex.state.plot %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = strategy, y = value, fill = state)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() +
  facet_wrap(.~sex) +
  # scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "Diabetes-related DALYs averted", fill = "Income quintile") + 
  scale_y_continuous(labels = scales::comma) 

incr.health.sex.state.prop <- incr.health.sex.state %>% 
  select(strategy, state, sex, dm_dalys) %>% 
  filter(dm_dalys < 0) %>% 
  group_by(strategy) %>% 
  mutate(prop = dm_dalys / sum(dm_dalys))
incr.dalys.sex.state.prop <- incr.health.sex.state.prop %>% ungroup %>% select(strategy, state, sex, prop) %>% spread(key = "state", value = "prop")

#### *************************************************************************************************************************
## Facet wrap for income quintile and sex -- which basically ends up being total DALYs averted and by sex
health_plot <- dplyr::bind_rows(incr.health, incr.health.sex) %>% 
  select(strategy, incr_dm_dalys, threshold, sex) %>% 
  mutate(sex = ifelse(is.na(sex), "Both", sex))
health_plot$strategy <- ordered(health_plot$strategy,
                                levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
health_plot$sex <- ordered(health_plot$sex,
                           levels = c("Both", "Male", "Female"))
health_plot$threshold <- factor(health_plot$threshold,
                                labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
health_plot$incr_dm_dalys <- health_plot$incr_dm_dalys * -1 # to convert to health effects _averted_

ggplot(health_plot %>% filter(strategy %in% c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy")), 
       aes(x = sex, y = incr_dm_dalys, fill = strategy)) + 
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
incr.oop20_80 <- results20_80 %>% 
  group_by(strategy) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>%
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop20_100 <- results20_100 %>% 
  group_by(strategy) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>%
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop30_80 <- results30_80 %>% 
  group_by(strategy) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>%
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop30_100 <- results30_100 %>% 
  group_by(strategy) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>%
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop <- rbind(incr.oop20_80, incr.oop20_100, incr.oop30_80, incr.oop30_100)

# Table
summary.oop <- rbind(incr.oop20_80, incr.oop20_100, incr.oop30_80, incr.oop30_100)
total.oop <- summary.oop %>% select(strategy, threshold, coverage, oop_total, incr_oop_total) # OOP incurred (not averted) as denominator

# Bar plot by strategy
results.oop.plot <- incr.oop 
results.oop.plot$strategy <- ordered(results.oop.plot$strategy,
                                     levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))

ggplot(data = results.oop.plot %>% filter(strategy != "Current standard"), aes(x = strategy, y = oop_total / n_hef_eligible)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_minimal() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = cbp1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Type of OOP spending", y = "Amount (2019 USD)", fill = "Strategy") 

# OOP by income quintile --------------------------------------------------
incr.oop.income1.20_80 <- results20_80 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.income2.20_80 <- results20_80 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.income1.20_100 <- results20_100 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.income2.20_100 <- results20_100 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.income1.30_80 <- results30_80 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.income2.30_80 <- results30_80 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.income1.30_100 <- results30_100 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 1) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.income2.30_100 <- results30_100 %>% group_by(strategy, income_quintile) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(income_quintile == 2) %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'],
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.income <- rbind(incr.oop.income1.20_80, incr.oop.income1.20_100, incr.oop.income1.30_80, incr.oop.income1.30_100,
                         incr.oop.income2.20_80, incr.oop.income2.20_100, incr.oop.income2.30_80, incr.oop.income2.30_100)

# Plot
incr.oop.income.plot <- incr.oop.income %>% gather(key = "oop", value = "value", -c(strategy, income_quintile, threshold, coverage, n_hef_eligible))
incr.oop.income.plot$strategy <- ordered(incr.oop.income.plot$strategy,
                                         levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.oop.income.plot$income_quintile <- ordered(incr.oop.income.plot$income_quintile,
                                                levels = c(1:5))
incr.oop.income.plot$threshold <- factor(incr.oop.income.plot$threshold, 
                                         levels = c("20%", "30%"),
                                         labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr.oop.income.plot$coverage <- factor(incr.oop.income.plot$coverage, levels = c("80%", "100%"),
                                        labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
incr.oop.income.plot$value <- incr.oop.income.plot$value * -1 # convert to OOP _averted_

#### *************************************************************************************************************************
# OOP spending averted per HEF eligible
ggplot(data = incr.oop.income.plot %>% filter(strategy != "Current standard" & oop == "incr_oop_total"), 
       aes(x = strategy, y = value / n_hef_eligible, fill = income_quintile)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_grid(threshold ~ coverage) +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  # scale_fill_manual(values = cbp1, drop = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Strategy", y = "OOP expenditures averted per person eligible for HEF (USD 2019)", fill = "Income quintile")
ggsave("figures/oop_by_income.png")


# OOP by sex --------------------------------------------------------------
incr.oop.male.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.female.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.male.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.female.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.oop.male.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.female.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.male.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Male") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.female.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
  summarise(
    oop_medical = sum(oop_medical) * inf.fct,
    oop_nonmedical = sum(oop_nonmedical) * inf.fct,
    oop_total = sum(oop_total) * inf.fct
  ) %>% 
  filter(sex == "Female") %>% 
  ungroup %>% 
  mutate(
    incr_oop_medical = oop_medical - oop_medical[strategy == 'base'],
    incr_oop_nonmedical = oop_nonmedical - oop_nonmedical[strategy == 'base'], 
    incr_oop_total = oop_total - oop_total[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.oop.sex <- rbind(incr.oop.male.20_80, incr.oop.female.20_80, incr.oop.male.20_100, incr.oop.female.20_100, 
                      incr.oop.male.30_80, incr.oop.female.30_80, incr.oop.male.30_100, incr.oop.female.30_100)

# Plot
incr.oop.sex.plot <- incr.oop.sex %>% gather(key = "oop", value = "value", -c(strategy, sex, threshold, coverage, n_hef_eligible))
incr.oop.sex.plot$strategy <- ordered(incr.oop.sex.plot$strategy,
                                      levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.oop.sex.plot$threshold <- factor(incr.oop.sex.plot$threshold, 
                                      levels = c("20%", "30%"),
                                      labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr.oop.sex.plot$coverage <- ordered(incr.oop.sex.plot$coverage, 
                                      levels = c("80%", "100%"),
                                      labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
incr.oop.sex.plot$value <- incr.oop.sex.plot$value * -1

#### *************************************************************************************************************************
# OOP spending averted per HEF eligible by sex
ggplot(data = incr.oop.sex.plot %>% filter(strategy != "Current standard" & oop == "incr_oop_total"), 
       aes(x = strategy, y = value / n_hef_eligible, fill = sex)) +
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

# Updating 'disposable' results (previous error in 'che' function, 'results_12056712L_2020-01-14_2028.csv')
# results20_80 <- results20_80 %>% mutate(
#   che10_disposable_result = che(oop_total, disposable_income)$che10,
#   che25_disposable_result = che(oop_total, disposable_income)$che25,
#   che40_disposable_result = che(oop_total, disposable_income)$che40
# )

# Summary of FRP ----------------------------------------------------------
incr.frp20_80 <- results20_80 %>% 
  group_by(strategy) %>% 
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
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp20_100 <- results20_100 %>% 
  group_by(strategy) %>% 
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
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp30_80 <- results30_80 %>% 
  group_by(strategy) %>% 
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
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp30_100 <- results30_100 %>% 
  group_by(strategy) %>% 
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
  mutate(
    incr_che10 = che10 - che10[strategy == 'base'],
    incr_che25 = che25 - che25[strategy == 'base'], 
    incr_che40 = che40 - che40[strategy == 'base'], 
    incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
    incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
    incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
    incr_pov = pov - pov[strategy == 'base'],
    incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp <- rbind(incr.frp20_80, incr.frp20_100, incr.frp30_80, incr.frp30_100)
che.frp <- incr.frp %>% select(strategy, threshold, coverage, che40, incr_che40)

# FRP by income quintile --------------------------------------------------
incr.frp.income1.20_80 <- results20_80 %>% group_by(strategy, income_quintile) %>% 
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
  filter(income_quintile == 1) %>% 
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
    threshold = "20%",
    coverage = "80%"
  )

incr.frp.income2.20_80 <- results20_80 %>% group_by(strategy, income_quintile) %>%
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
  filter(income_quintile == 2) %>%
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
    threshold = "20%",
    coverage = "80%"
  )

incr.frp.income1.20_100 <- results20_100 %>% group_by(strategy, income_quintile) %>% 
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
  filter(income_quintile == 1) %>% 
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
    threshold = "20%",
    coverage = "100%"
  )

incr.frp.income2.20_100 <- results20_100 %>% group_by(strategy, income_quintile) %>%
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
  filter(income_quintile == 2) %>%
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
    threshold = "20%",
    coverage = "100%"
  )

incr.frp.income1.30_80 <- results30_80 %>% group_by(strategy, income_quintile) %>% 
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
  filter(income_quintile == 1) %>% 
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
    threshold = "30%",
    coverage = "80%"
  )

incr.frp.income2.30_80 <- results30_80 %>% group_by(strategy, income_quintile) %>%
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
  filter(income_quintile == 2) %>%
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
    threshold = "30%",
    coverage = "80%"
  )

incr.frp.income1.30_100 <- results30_100 %>% group_by(strategy, income_quintile) %>% 
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
  filter(income_quintile == 1) %>% 
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
    threshold = "30%",
    coverage = "100%"
  )

incr.frp.income2.30_100 <- results30_100 %>% group_by(strategy, income_quintile) %>%
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
  filter(income_quintile == 2) %>%
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
    threshold = "30%",
    coverage = "100%"
  )

incr.frp.income <- rbind(incr.frp.income1.20_80, incr.frp.income2.20_80,
                         incr.frp.income1.20_100, incr.frp.income2.20_100,
                         incr.frp.income1.30_80, incr.frp.income2.30_80,
                         incr.frp.income1.30_100, incr.frp.income2.30_100)

# Plot
incr.frp.income.plot <- incr.frp.income %>% gather(key = "frp", value = "value", -c(strategy, income_quintile, threshold, coverage))
incr.frp.income.plot$strategy <- ordered(incr.frp.income.plot$strategy,
                                         levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                         labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.frp.income.plot$income_quintile <- ordered(incr.frp.income.plot$income_quintile,
                                                levels = c(1:5))
incr.frp.income.plot$threshold <- factor(incr.frp.income.plot$threshold, 
                                         levels = c("20%", "30%"),
                                         labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr.frp.income.plot$coverage <- ordered(incr.frp.income.plot$coverage, 
                                         levels = c("80%", "100%"),
                                         labels = c("HEF coverage: 80% of direct medical costs", "HEF coverage: 100% of direct medical costs"))
incr.frp.income.plot$frp <- ordered(incr.frp.income.plot$frp, 
                                    levels = c("incr_che10", "incr_che25", "incr_che40", "incr_che10_disposable", "incr_che25_disposable", "incr_che40_disposable", "incr_pov", "incr_pov_disposable"),
                                    labels = c("CHE threshold: 10%", "CHE threshold: 25%", "CHE threshold: 40%", 
                                               "CHE threshold: 10% \n(disposable income)", "CHE threshold: 25% \n(disposable income)", "CHE threshold: 40% \n(disposable income)",
                                               "Impoverishment", "Impoverishment \n(disposable income)"))
incr.frp.income.plot$value <- incr.frp.income.plot$value * -1

#### *************************************************************************************************************************
# CHE averted
ggplot(data = incr.frp.income.plot %>% filter(strategy != "Current standard" & income_quintile == 1 & frp %in% c("CHE threshold: 40%")), 
       aes(x = strategy, y = value)) +
  geom_bar(stat = "identity") + 
  facet_grid(threshold~coverage) +
  theme_bw() + 
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Coverage strategy", y = "Cases of CHE averted") + 
  scale_y_continuous(labels = scales::comma) 
ggsave("figures/che.png")

# FRP by sex --------------------------------------------------------------
incr.frp.male.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Male") %>% 
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
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp.female.20_80 <- results20_80 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Female") %>% 
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
    threshold = "20%",
    coverage = "80%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp.male.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Male") %>% 
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
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp.female.20_100 <- results20_100 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Female") %>% 
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
    threshold = "20%",
    coverage = "100%",
    n_hef_eligible = population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct
  )

incr.frp.male.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Male") %>% 
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
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.frp.female.30_80 <- results30_80 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Female") %>% 
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
    threshold = "30%",
    coverage = "80%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.frp.male.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Male") %>% 
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
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.frp.female.30_100 <- results30_100 %>% 
  group_by(strategy, sex) %>% 
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
  filter(sex == "Female") %>% 
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
    threshold = "30%",
    coverage = "100%",
    n_hef_eligible = population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct
  )

incr.frp.sex <- rbind(incr.frp.male.20_80, incr.frp.female.20_80, incr.frp.male.20_100, incr.frp.female.20_100,
                      incr.frp.male.30_80, incr.frp.female.30_80, incr.frp.male.30_100, incr.frp.female.30_100)

# Plot
incr.frp.sex.plot <- incr.frp.sex %>% gather(key = "frp", value = "value", -c(strategy, sex, threshold, coverage, n_hef_eligible))
incr.frp.sex.plot$strategy <- ordered(incr.frp.sex.plot$strategy,
                                      levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                      labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
incr.frp.sex.plot$threshold <- factor(incr.frp.sex.plot$threshold, 
                                      levels = c("20%", "30%"),
                                      labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
incr.frp.sex.plot$coverage <- ordered(incr.frp.sex.plot$coverage, 
                                         levels = c("80%", "100%"),
                                         labels = c("HEF coverage: 80%", "HEF coverage: 100%"))
incr.frp.sex.plot$frp <- ordered(incr.frp.sex.plot$frp, 
                                 levels = c("incr_che10", "incr_che25", "incr_che40", "incr_che10_disposable", "incr_che25_disposable", "incr_che40_disposable", "incr_pov", "incr_pov_disposable"),
                                 labels = c("CHE threshold: 10%", "CHE threshold: 25%", "CHE threshold: 40%", 
                                            "CHE threshold: 10% \n(disposable income)", "CHE threshold: 25% \n(disposable income)", "CHE threshold: 40% \n(disposable income)",
                                            "Impoverishment", "Impoverishment \n(disposable income)"))
incr.frp.sex.plot$value <- incr.frp.sex.plot$value * -1

#### *************************************************************************************************************************
# CHE averted by sex
ggplot(data = incr.frp.sex.plot %>% filter(strategy != "Current standard" & (frp %in% c("CHE threshold: 40%"))), 
       aes(x = strategy, y = value, fill = sex)) +
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
icer.dalys <- cbind(total.costs, dalys.health[3:4]) %>% 
  mutate(
    icer.dalys = incr_cost_total / -incr_dm_dalys
  )
xtable::xtable(icer.dalys)

icer.dalys.sex <- merge(incr.costs.sex, incr.health.sex, by = c("strategy", "sex", "threshold")) %>% 
  mutate(
    icer.dalys.sex = incr_cost_total / -incr_dm_dalys
  )


### ICER PLOT 
# DALYs
icer.dalys.plot <- icer.dalys
icer.dalys.plot <- icer.dalys.plot %>% mutate(
  n_hef_eligible = ifelse(threshold == "20%", 
                          population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct,
                          population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct)
  )
icer.dalys.plot$strategy <- ordered(icer.dalys.plot$strategy,
                                    levels = c("screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
icer.dalys.plot$threshold <- factor(icer.dalys.plot$threshold, 
                                      levels = c("20%", "30%"),
                                      labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
ggplot(icer.dalys.plot, aes(x = incr_dm_dalys, y = incr_cost_total / n_hef_eligible)) + 
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
icer.dalys.sex.plot <- dplyr::bind_rows(icer.dalys, icer.dalys.sex) %>% 
  select(strategy, sex, threshold, incr_cost_total, incr_dm_dalys) %>% 
  mutate(
    sex = ifelse(is.na(sex), "Total", sex)
  )
icer.dalys.sex.plot <- icer.dalys.sex.plot %>% mutate(
  n_hef_eligible = ifelse(threshold == "20%", 
                          population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct,
                          population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct)
)
icer.dalys.sex.plot$strategy <- ordered(icer.dalys.sex.plot$strategy,
                                    levels = c("screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
icer.dalys.sex.plot$threshold <- factor(icer.dalys.sex.plot$threshold, 
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

ggplot(icer.dalys.sex.plot %>% filter(!is.na(strategy)), aes(x = incr_dm_dalys * -1, y = incr_cost_total / n_hef_eligible)) + 
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
incr.oop.sex
icer.oop.sex <- merge(incr.costs.sex, incr.oop.sex, by = c("strategy", "sex", "threshold", "coverage")) %>% 
  mutate(
    icer.oop.sex = incr_cost_total / incr_oop_total
  ) %>% 
  select(strategy, sex, threshold, coverage, incr_oop_total, icer.oop.sex)

# icer.oop.plot <- dplyr::bind_rows(icer.oop.income, icer.oop.sex) %>% 
#   mutate(
#     category = ifelse(!is.na(income_quintile), "Poorest income quintile", sex)
#   ) %>% 
#   select(-c(income_quintile, sex, cost_diagnosis, cost_tx, cost_complications, oop_medical, oop_nonmedical))
# icer.oop.plot$n_hef_beneficiaries <- ifelse(icer.oop.plot$category == "Poorest income quintile", n_hef_beneficiaries20_80,
#                                             ifelse(icer.oop.plot$category == "Male", n_hef_beneficiaries20_80.sex[2, 2][[1]], n_hef_beneficiaries20_80.sex[1, 2][[1]]))
# icer.oop.plot$strategy <- ordered(icer.oop.plot$strategy,
#                                   levels = c("screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
#                                   labels = c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
# 
# ggplot(icer.oop.plot, aes(x = oop_total / 1e6, y = cost_total / n_hef_beneficiaries)) + 
#   geom_point(aes(color = strategy, shape = coverage), size = 4) + 
#   facet_grid(category ~ .) + 
#   theme_bw() + 
#   scale_color_npg(drop = FALSE) + 
#   scale_shape_manual(values = c(16, 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14)) + 
#   labs(x = "OOP expenditures averted (millions, USD 2019)", y = "Incremental costs per HEF beneficiary (USD 2019)", color = "Coverage strategy", shape = "HEF coverage") + 
#   scale_y_continuous(labels = scales::comma) + 
#   scale_x_continuous(labels = scales::comma) + 
#   geom_hline(aes(yintercept = 0)) +
#   geom_vline(aes(xintercept = 0))

## Financial risk protection
icer.che <- merge(total.costs, che.frp, by = c("strategy", "threshold")) %>% 
  mutate(
    icer.che = incr_cost_total / incr_che40
  )
icer.che <- icer.che[, c("strategy", "threshold", "coverage", "cost_total", "incr_cost_total", "che40", "incr_che40", "icer.che")]
xtable::xtable(icer.che)

icer.frp.sex <-  merge(incr.costs.sex, incr.frp.sex, by = c("strategy", "sex", "threshold", "coverage")) %>% 
  select(strategy, sex, threshold, coverage, incr_cost_total, incr_che40, n_hef_eligible) %>% 
  mutate(
    icer.che.sex = incr_cost_total / incr_che40
  )

icer.frp.plot <- dplyr::bind_rows(icer.che, icer.frp.sex) %>% 
  select(strategy, sex, threshold, coverage, incr_cost_total, incr_che40, n_hef_eligible) %>% 
  mutate(
    sex = ifelse(is.na(sex), "Total", sex),
    n_hef_eligible =  ifelse(threshold == "20%", population20 %>% filter(income < hef_threshold20) %>% nrow() * inf.fct, 
                             ifelse(threshold == "30%", population30 %>% filter(income < hef_threshold30) %>% nrow() * inf.fct, n_hef_eligible))
  )  
icer.frp.plot$strategy <- ordered(icer.frp.plot$strategy,
                                  levels = c("screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                  labels = c("Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
# icer.frp.plot$frp <- ordered(icer.frp.plot$frp, 
#                              levels = c("incr_che10", "incr_che25", "incr_che40", "incr_pov"),
#                              labels = c("Catastrophic health expenditure averted \n(10% threshold)", "Catastrophic health expenditure averted \n(25% threshold)", "Catastrophic health expenditure averted \n(40% threshold)", "Cases of poverty averted"))
icer.frp.plot$threshold <- factor(icer.frp.plot$threshold, 
                                  levels = c("20%", "30%"),
                                  labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))

ggplot(icer.frp.plot %>%filter(!is.na(strategy)), aes(x = incr_che40 * -1, y = incr_cost_total / n_hef_eligible)) + 
  geom_point(aes(color = strategy, shape = coverage), size = 4) + 
  facet_grid(sex ~ threshold, scales = "free_x") + 
  theme_bw() + 
  scale_color_npg(drop = FALSE) + 
  scale_shape_manual(values = c(16, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  labs(x = "Cases of CHE averted", y = "Incremental costs per person eligible for HEF (USD 2019)", color = "Strategy", shape = "HEF Coverage") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0)) 
ggsave("figures/icer_frp.png")


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


