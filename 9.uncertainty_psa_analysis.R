##################################
#### Extended cost-effectiveness analysis of interventions for diabetes under diverse health coverage in Cambodia
#### 7.0 Probabilistic sensitivity analysis -- analysis of results
#### Copyright Isabelle Feldhaus
#### 27 January 2020
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


# 1. Load results ---------------------------------------------------------

psa_results <- readr::read_csv('output/psa_results_1050L_2020-01-30_1356.csv')
bootstrapped_outcomes <- readr::read_csv('output/bootstrapped_outcomes_105000L_2020-01-30_1356.csv')

psa_results$strategy <- ordered(psa_results$strategy,
                                levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug \ntherapy + complications"))

bootstrapped_outcomes$strategy <- ordered(bootstrapped_outcomes$strategy,
                                          levels = c("base", "screen_only", "tx_only", "screen_tx", "comp_only", "tx_comp", "screen_tx_comp"),
                                          labels = c("Current standard", "Diagnostics only", "Drug therapy only", "Diagnostics + drug therapy", "Complications only", "Drug therapy + complications", "Diagnostics + drug \ntherapy + complications"))

combined_uncertainty <- dplyr::bind_rows(psa_results, bootstrapped_outcomes)

ggplot(combined_uncertainty %>% filter(sex == "Total") %>% filter(strategy != "Current standard"), 
       aes(x = -incr_dm_dalys / 1e3, y = incr_cost_total / n_hef_eligible, color = strategy)) + 
  geom_point(alpha = 0.2) + 
  facet_grid(. ~ strategy) + 
  scale_color_manual(values = npg_palette[2:7]) + 
  stat_ellipse(color = "black", alpha = 0.6) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") + 
  labs(x = "Incremental diabetes-related DALYs averted (thousands)", y = "Incremental total costs per person eligible for HEF (2019 USD)") +
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty_dalys.png")

ggplot(combined_uncertainty %>% filter(sex == "Total") %>% filter(strategy != "Current standard"), 
       aes(x = -incr_che40 / 1e3, y = incr_cost_total / n_hef_eligible, color = strategy)) + 
  geom_point(alpha = 0.2) + 
  facet_grid(. ~ strategy) + 
  scale_color_manual(values = npg_palette[2:7]) + 
  stat_ellipse(color = "black", alpha = 0.6) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") + 
  labs(x = "Incremental diabetes-related CHE cases averted (thousands)", y = "Incremental total costs per person eligible for HEF (2019 USD)") +
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty_che.png")

combined_uncertainty_plot <- combined_uncertainty %>% 
  select(strategy, sex, incr_cost_total, incr_dm_dalys, incr_che40, n_hef_eligible) %>% 
  gather(key = "outcome", value = "value", -c(strategy, sex, n_hef_eligible, incr_cost_total))
combined_uncertainty_plot$outcome <- ordered(combined_uncertainty_plot$outcome,
                                             levels = c("incr_dm_dalys", "incr_che40"),
                                             labels = c("DALYs", "CHE cases"))

ggplot(combined_uncertainty_plot %>% filter(sex == "Total") %>% filter(strategy != "Current standard"), 
       aes(x = -value / 1e3, y = incr_cost_total / n_hef_eligible, color = strategy)) + 
  geom_point(alpha = 0.2) + 
  facet_grid(strategy ~ outcome) + 
  scale_color_manual(values = npg_palette[2:7]) + 
  stat_ellipse(color = "black", alpha = 0.6) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 12),
        legend.position = "none") + 
  labs(x = "Incremental effects averted (thousands)", y = "Incremental total costs per person eligible for HEF (2019 USD)") +
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::comma) + 
  geom_hline(aes(yintercept = 0), alpha = 0.4) +
  geom_vline(aes(xintercept = 0), alpha = 0.4)
ggsave("figures/uncertainty2.png", width = 180, height = 375, units = "mm")

