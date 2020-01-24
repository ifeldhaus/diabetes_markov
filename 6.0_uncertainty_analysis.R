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
bootstrapped_outcomes20 <- load_object_from_rdata('output/bootstrapped_outcomes_total_20L_2020-01-23_0211.Rdata') # 20% threshold
bootstrapped_outcomes30 <- load_object_from_rdata('output/bootstrapped_outcomes_total_30L_2020-01-23_1213.Rdata')

# Load analysis results
results20 <- readr::read_csv('output/results_24034003L_2020-01-22_1646.csv') # complete results for HEF thresholds 20% and 30% (100% coverage)
population20 <- readr::read_csv("output/ichar_200000L_2020-01-19_2041.csv")
results_pop20 <- merge(population20, results20, by = "id")
outcomes20 <- compute_outcomes(population20, results_pop20, hef_threshold = 0.2)

results30 <- readr::read_csv('output/results_24033414L_2020-01-22_1707.csv')
population30 <- readr::read_csv("output/ichar_200000L_2020-01-19_2306.csv")
results_pop30 <- merge(population30, results30, by = "id")
outcomes30 <- compute_outcomes(population30, results_pop30, hef_threshold = 0.3) 

# Combine totals and sex-specific outcomes
outcomes20_plot <- outcomes20
outcomes20_plot$threshold <- "20%"
outcomes30_plot <- outcomes30
outcomes30_plot$threshold <- "30%"
outcomes_plot <- rbind(outcomes20_plot, outcomes30_plot)

# Plots
bootstrapped_outcomes20$strategy <- ordered(bootstrapped_outcomes20$strategy,
                                          levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                          labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
bootstrapped_outcomes30$strategy <- ordered(bootstrapped_outcomes30$strategy,
                                            levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                            labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
outcomes20_plot$strategy <- ordered(outcomes20_plot$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
outcomes20_plot$threshold <- factor(outcomes20_plot$threshold, 
                                    levels = c("20%", "30%"),
                                    labels = c("Poorest 20% eligible for HEF", "Poorest 30% eligible for HEF"))
outcomes30_plot$strategy <- ordered(outcomes30_plot$strategy,
                                    levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                    labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug therapy + complications"))
outcomes30_plot$threshold <- factor(outcomes30_plot$threshold, 
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
  geom_point(data = outcomes20_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 20% eligible for HEF"), shape = 20, color = "black") +
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
# ggsave("figures/uncertainty20_health.png")

ggplot(bootstrapped_outcomes20 %>% filter(strategy != "Current standard"), aes(x = -incr_che40 / 1e3, y = incr_cost_total / 1e6, colour = strategy)) +
  geom_point(shape = 1) + 
  stat_ellipse(color = "black") +
  geom_point(data = outcomes20_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 20% eligible for HEF"), shape = 20, color = "black") +
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
  geom_point(data = outcomes30_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 30% eligible for HEF"), shape = 20, color = "black") +
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
  geom_point(data = outcomes30_plot %>% filter(strategy != "Current standard") %>% filter(threshold == "Poorest 30% eligible for HEF"), shape = 20, color = "black") +
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


