##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 5.0 Uncertainty analysis
#### Copyright Isabelle Feldhaus
#### 21 January 2020
##################################

# 0. Setup ----------------------------------------------------------------

setwd("/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/diabetes_markov/")

OUTPUT_DIR = 'output'

library(tidyverse)
library(markovchain) 
library(readstata13)
library(parallel)
library(foreach)
library(doParallel) # alternative: `doMC`
library(doMC)
library(scales)
source('helpers.R')
source('1.outcomes.R')


# 1. Uncertainty analysis plot --------------------------------------------

# Load bootstrapped results
bootstrapped_outcomes20 <- load_object_from_rdata('output/bootstrapped_outcomes_total_20L_2020-01-26_0133.Rdata') # 20% threshold
bootstrapped_outcomes30 <- load_object_from_rdata('output/bootstrapped_outcomes_total_30L_2020-01-26_1327.Rdata') # 30% threshold

# Load analysis results
population <- readr::read_csv("output/ichar_200000L_2020-01-25_1855.csv")
person_cycle_x <- load_object_from_rdata('output/person_cycle_x_200000L_2020-01-25_1855.RData')
population20 <- add_hef_to_pop(population, person_cycle_x, hef_quantile = 0.2) %>% mutate(threshold = "20%")
population30 <- add_hef_to_pop(population, person_cycle_x, hef_quantile = 0.3) %>% mutate(threshold = "30%")

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

# Compute outcomes
outcomes20_80 <- compute_outcomes(population20, results20_80, hef_threshold = 0.2)
outcomes20_100 <- compute_outcomes(population20, results20_100, hef_threshold = 0.2)
outcomes30_80 <- compute_outcomes(population30, results30_80, hef_threshold = 0.3) 
outcomes30_100 <- compute_outcomes(population30, results30_100, hef_threshold = 0.3) 

# Combine totals and sex-specific outcomes
outcomes20_80plot <- outcomes20_80 %>% mutate(threshold = "20%", coverage = "80%")
outcomes20_100plot <- outcomes20_100 %>% mutate(threshold = "20%", coverage = "100%")
outcomes30_80plot <- outcomes30_80 %>% mutate(threshold = "30%", coverage = "80%")
outcomes30_100plot <- outcomes30_100 %>% mutate(threshold = "30%", coverage = "100%")
outcomes_plot <- rbind(outcomes20_80plot, outcomes20_100plot, outcomes30_80plot, outcomes30_100plot)

# Plots
bootstrapped_outcomes20$strategy <- ordered(bootstrapped_outcomes20$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
bootstrapped_outcomes30$strategy <- ordered(bootstrapped_outcomes30$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
outcomes_plot$strategy <- ordered(outcomes_plot$strategy,
                                  levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                  labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
outcomes_plot$threshold <- factor(outcomes_plot$threshold, 
                                  levels = c("20%", "30%"),
                                  labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))


# plot_1k_boots_20k_people <- ggplot(bootstrapped_outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), aes(x = -incr_dm_dalys / n_hef_eligible, y = incr_cost_total / n_hef_eligible, colour = strategy)) +
#   geom_point(shape = 1) + 
#   stat_ellipse(color = "black") +
#   geom_point(data = outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), shape = 20, color = "black") +
#   facet_grid(threshold ~ strategy) +
#   theme_bw() 
# # ggsave('figures/uncertainty.png')
# 
# plot_100_boots_50k_people <- ggplot(bootstrapped_outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), aes(x = -incr_dm_dalys / n_hef_eligible, y = incr_cost_total / n_hef_eligible, colour = strategy)) +
#   geom_point(shape = 1) + 
#   stat_ellipse(color = "black") +
#   geom_point(data = outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), shape = 20, color = "black") +
#   facet_grid(threshold ~ strategy) +
#   theme_bw() 
# 
# plot_100_boots_50k_people_wor <- ggplot(bootstrapped_outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), aes(x = -incr_dm_dalys / n_hef_eligible, y = incr_cost_total / n_hef_eligible, colour = strategy)) +
#   geom_point(shape = 1) + 
#   stat_ellipse(color = "black") +
#   geom_point(data = outcomes_plot %>% filter(strategy != "base") %>% filter(sex == "Total"), shape = 20, color = "black") +
#   facet_grid(threshold ~ strategy) +
#   theme_bw() 

# 20% threshold
ggplot(bootstrapped_outcomes20 %>% filter(strategy != "Current standard"), aes(x = -incr_dm_dalys / 1e3, y = incr_cost_total / 1e6, colour = strategy)) +
  geom_point(shape = 1) + 
  stat_ellipse(color = "black") +
  geom_point(data = outcomes_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 20% eligible for HEF"), shape = 20, color = "black") +
  scale_color_npg() +
  facet_grid(sex ~ strategy, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") + 
  labs(x = "Incremental diabetes-related DALYs averted (thousands)", y = "Incremental total costs (millions)") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty20_health.png")

ggplot(bootstrapped_outcomes20 %>% filter(strategy != "Current standard"), aes(x = -incr_che40 / 1e3, y = incr_cost_total / 1e6, colour = strategy)) +
  geom_point(shape = 1) + 
  stat_ellipse(color = "black") +
  geom_point(data = outcomes_plot %>% filter(strategy != "Current standard") %>% 
               filter(threshold == "Poorest 20% eligible for HEF") %>% 
               filter(coverage == "80%"), shape = 20, color = "black") +
  scale_color_npg() +
  facet_grid(sex ~ strategy, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Incremental cases of CHE averted (thousands)", y = "Incremental total costs (millions)") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty20_frp.png")

# 30% threshold
ggplot(bootstrapped_outcomes30 %>% filter(strategy != "Current standard"), aes(x = -incr_dm_dalys / 1e3, y = incr_cost_total / 1e6, colour = strategy)) +
  geom_point(shape = 1) + 
  stat_ellipse(color = "black") +
  geom_point(data = outcomes_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 30% eligible for HEF"), shape = 20, color = "black") +
  scale_color_npg() +
  facet_grid(sex ~ strategy, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") + 
  labs(x = "Incremental diabetes-related DALYs averted (thousands)", y = "Incremental total costs (millions)") + 
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty30_health.png")

ggplot(bootstrapped_outcomes30 %>% filter(strategy != "Current standard"), aes(x = -incr_che40 / 1e3, y = incr_cost_total / 1e6, colour = strategy)) +
  geom_point(shape = 1) + 
  stat_ellipse(color = "black") +
  geom_point(data = outcomes_plot %>% filter(strategy != "Current standard") %>% 
               filter(threshold == "Poorest 30% eligible for HEF") %>% 
               filter(coverage == "80%"), shape = 20, color = "black") +
  scale_color_npg() +
  facet_grid(sex ~ strategy, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Incremental cases of CHE averted (thousands)", y = "Incremental total costs (millions)") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty30_frp.png")


