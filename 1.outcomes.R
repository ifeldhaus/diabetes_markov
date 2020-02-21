# Functions to compute outcomes


compute_outcomes <- function(population, results_pop, hef_threshold) {
  
  hef_threshold_amt <- quantile(population$income, hef_threshold)[[1]]
  inf.fct <- 16e6 / nrow(population)
  
  outcomes <- results_pop %>% 
    group_by(strategy) %>% 
    summarize(
      cost_diagnosis = sum(cost_diagnosis) * inf.fct,
      cost_tx = sum(cost_tx) * inf.fct,
      cost_complications = sum(cost_complications) * inf.fct, 
      cost_total = sum(cost_total) * inf.fct,
      dm_deaths = sum(state == "dm_death") * inf.fct,
      dm_dalys = sum(dalys_dm) * inf.fct,
      che10 = sum(che10_result) * inf.fct,
      che25 = sum(che25_result) * inf.fct,
      che40 = sum(che40_result) * inf.fct,
      pov = sum(pov_result, na.rm = T) * inf.fct,
      che10_disposable = sum(che10_disposable_result, na.rm = T) * inf.fct,
      che25_disposable = sum(che25_disposable_result, na.rm = T) * inf.fct,
      che40_disposable = sum(che40_disposable_result, na.rm = T) * inf.fct,
      pov_disposable = sum(pov_disposable_result, na.rm = T) * inf.fct
    ) %>% 
    ungroup %>% 
    mutate(
      incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
      incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
      incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
      incr_cost_total = cost_total - cost_total[strategy == 'base'],
      incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
      incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
      incr_che10 = che10 - che10[strategy == 'base'],
      incr_che25 = che25 - che25[strategy == 'base'], 
      incr_che40 = che40 - che40[strategy == 'base'], 
      incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
      incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
      incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
      incr_pov = pov - pov[strategy == 'base'],
      incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base'],
      n_hef_eligible = sum(population$income < hef_threshold_amt) * inf.fct
    )
  
  outcomes_by_sex <- results_pop %>% 
    group_by(strategy, sex) %>% 
    summarize(
      cost_diagnosis = sum(cost_diagnosis) * inf.fct,
      cost_tx = sum(cost_tx) * inf.fct,
      cost_complications = sum(cost_complications) * inf.fct, 
      cost_total = sum(cost_total) * inf.fct,
      dm_deaths = sum(state == "dm_death") * inf.fct,
      dm_dalys = sum(dalys_dm) * inf.fct,
      che10 = sum(che10_result) * inf.fct,
      che25 = sum(che25_result) * inf.fct,
      che40 = sum(che40_result) * inf.fct,
      pov = sum(pov_result, na.rm = T) * inf.fct,
      che10_disposable = sum(che10_disposable_result, na.rm = T) * inf.fct,
      che25_disposable = sum(che25_disposable_result, na.rm = T) * inf.fct,
      che40_disposable = sum(che40_disposable_result, na.rm = T) * inf.fct,
      pov_disposable = sum(pov_disposable_result, na.rm = T) * inf.fct
    ) %>% 
    group_by(sex) %>% 
    mutate(
      incr_cost_diagnosis = cost_diagnosis - cost_diagnosis[strategy == 'base'],
      incr_cost_tx = cost_tx - cost_tx[strategy == 'base'], 
      incr_cost_complications = cost_complications - cost_complications[strategy == 'base'], 
      incr_cost_total = cost_total - cost_total[strategy == 'base'],
      incr_dm_deaths = dm_deaths - dm_deaths[strategy == 'base'],
      incr_dm_dalys = dm_dalys - dm_dalys[strategy == 'base'],
      incr_che10 = che10 - che10[strategy == 'base'],
      incr_che25 = che25 - che25[strategy == 'base'], 
      incr_che40 = che40 - che40[strategy == 'base'], 
      incr_che10_disposable = che10_disposable - che10_disposable[strategy == 'base'],
      incr_che25_disposable = che25_disposable - che25_disposable[strategy == 'base'],
      incr_che40_disposable = che40_disposable - che40_disposable[strategy == 'base'],
      incr_pov = pov - pov[strategy == 'base'],
      incr_pov_disposable = pov_disposable - pov_disposable[strategy == 'base']
    )
  
  all_outcomes <- merge(
    outcomes_by_sex,
    population %>% group_by(sex) %>% summarize( n_hef_eligible = sum(income < hef_threshold_amt) * inf.fct),
    by = 'sex'
  ) %>% rbind(
    outcomes %>% mutate(sex = 'Total')
  )
  
  return(all_outcomes)
  
}


