# Compute costs & OOP

discount_costs <- function(value, dr, time) {
  final_amt <- value * (1 / (1 + dr)^time)
  return(final_amt)
}

compute_costs <- function(strategy, current_state, next_state, p_provider, outpatient_costs, hospitalization_costs,
                          p_op, p_hosp, hef, hef_utilization, dr, i_cycle) {
  
  if (next_state == "healthy" | next_state == "diabetes" | next_state == "dm_death" | next_state == "other_death") {
    cost_diagnosis <- 0
    cost_tx <- 0
    cost_complications <- 0
    oop_medical <- 0
    oop_nonmedical <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
  }
  
  if (next_state == "notx") { 
    if (current_state == "healthy" | current_state == "diabetes") { # incidence of diabetes, diagnosed
      # Determine where diagnosis occurred depending on careseeking behavior
      # ASSUMPTION: if previously in `healthy` or `undiagnosed` state, likely diagnosed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
      
      if (strategy != "base" & hef == 1 & hef_utilization == 1) {
        oop_medical <- cost_diagnosis * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle)
      } else {
        oop_medical <- cost_diagnosis
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle)
      }
      
    } else { # in all other cases, patient has been previously diagnosed and given diet consultation (no drug treatment)
      cost_diagnosis <- 0
      oop_medical <- 0
      oop_nonmedical <- 0
    }
    cost_tx <- 0
    cost_complications <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
  }
  
  if (next_state == "oad") {
    if (current_state == "healthy" | current_state == "diabetes") { 
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
    } 
    if (current_state == "notx") {
      # ASSUMPTION: change in therapy prescribed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_laboratory + c_outpatient, dr, i_cycle) 
    }
    else { 
      cost_diagnosis <- 0 # Assumes no further outpatient visits nor regular laboratory testing coming from other states
    }
    cost_tx <- discount_costs(c_oad, dr, i_cycle) 
    cost_complications <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage) 
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    }
  }
  
  if (next_state == "ins") {
    if (current_state == "healthy" | current_state == "diabetes") { 
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
    } 
    if (current_state == "notx") {
      # ASSUMPTION: change in therapy prescribed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_laboratory + c_outpatient, dr, i_cycle)  
    }
    else { 
      cost_diagnosis <- 0 # Assumes no further outpatient visits nor regular laboratory testing coming from other states
    }
    cost_tx <- discount_costs(c_ins, dr, i_cycle) 
    cost_complications <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage)
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    }
  }
  
  if (next_state == "nephro") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_nephro, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_nephro, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications 
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage)
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
    }
  }
  
  if (next_state == "retino") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_retino, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_retino, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (op_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "neuro") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_neuro, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_neuro, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (op_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "ang") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_ang, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_ang, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (op_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "pvd") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_pvd, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, i_cycle) 
        cost_complications <- discount_costs(c_pvd, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (op_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "mi") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_mi, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_mi, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications 
    
    if (hosp_visit == 1) {
      if (hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "stroke") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_stroke, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_stroke, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (hosp_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (next_state == "fail") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (current_state == "healthy" | current_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_fail, dr, i_cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, i_cycle) 
        cost_complications <- discount_costs(c_fail, dr, i_cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (hosp_visit == 1) {
      if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
        oop_medical <- cost_total * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, i_cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  oop_total <- oop_medical + oop_nonmedical
  costs <- list(
    "cost_diagnosis" = cost_diagnosis,
    "cost_tx" = cost_tx, 
    "cost_complications" = cost_complications,
    "cost_total"= cost_total,
    "oop_medical" = oop_medical,
    "oop_nonmedical" = oop_nonmedical, 
    "oop_total" = oop_total
  )
  
  return(costs)
}
