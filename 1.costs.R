# Compute costs & OOP

discount_costs <- function(value, dr, time) {
  final_amt <- value * (1 / (1 + dr)^time)
  return(final_amt)
}

compute_costs <- function(init_state, end_state, p_provider, outpatient_costs, hospitalization_costs,
                          p_op, p_hosp, hef, hef_utilization, dr, cycle) {
  
  if (end_state == "healthy" | end_state == "diabetes" | end_state == "dm_death" | end_state == "other_death") {
    cost_diagnosis <- 0
    cost_tx <- 0
    cost_complications <- 0
    oop_medical <- 0
    oop_nonmedical <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
  }
  
  if (end_state == "notx") { 
    if (init_state == "healthy" | init_state == "diabetes") { # incidence of diabetes, diagnosed
      # Determine where diagnosis occurred depending on careseeking behavior
      # ASSUMPTION: if previously in `healthy` or `undiagnosed` state, likely diagnosed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
      
      if (strategy != "base" & hef == 1 & hef_utilization == 1) {
        oop_medical <- cost_diagnosis * (1 - coverage)
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle)
      } else {
        oop_medical <- cost_diagnosis
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle)
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
  
  if (end_state == "oad") {
    if (init_state == "healthy" | init_state == "diabetes") { 
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
    } 
    if (init_state == "notx") {
      # ASSUMPTION: change in therapy prescribed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_laboratory + c_outpatient, dr, cycle) 
    }
    else { 
      cost_diagnosis <- 0 # Assumes no further outpatient visits nor regular laboratory testing coming from other states
    }
    cost_tx <- discount_costs(c_oad, dr, cycle) 
    cost_complications <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage) 
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    }
  }
  
  if (end_state == "ins") {
    if (init_state == "healthy" | init_state == "diabetes") { 
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
    } 
    if (init_state == "notx") {
      # ASSUMPTION: change in therapy prescribed through outpatient care
      provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
      c_outpatient <- outpatient_costs[provider]
      cost_diagnosis <- discount_costs(c_laboratory + c_outpatient, dr, cycle)  
    }
    else { 
      cost_diagnosis <- 0 # Assumes no further outpatient visits nor regular laboratory testing coming from other states
    }
    cost_tx <- discount_costs(c_ins, dr, cycle) 
    cost_complications <- 0
    cost_total <- cost_diagnosis + cost_tx + cost_complications
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage)
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    }
  }
  
  if (end_state == "nephro") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_nephro, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_nephro, dr, cycle) 
      }
    } else {
      cost_diagnosis <- 0
      cost_complications <- 0
    }
    cost_tx <- 0 # OAD or insulin costs to be extended as necessary upon calculation of final results 
    cost_total <- cost_diagnosis + cost_tx + cost_complications 
    
    if (strategy != "base" & hef == 1 & hef_utilization == 1 & cost_diagnosis != 0) {
      oop_medical <- cost_total * (1 - coverage)
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    } else {
      oop_medical <- cost_total
      oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
    }
  }
  
  if (end_state == "retino") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_retino, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_retino, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "neuro") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_neuro, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_neuro, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "ang") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_ang, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_ang, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "pvd") {
    op_visit <- rbinom(n = 1, size = 1, p = p_op)
    if (op_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_pvd, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_outpatient <- outpatient_costs[provider]
        cost_diagnosis <- discount_costs(c_outpatient, dr, cycle) 
        cost_complications <- discount_costs(c_pvd, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_op, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "mi") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_mi, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_mi, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "stroke") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_stroke, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_stroke, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
      }
    } else {
      oop_medical <- 0
      oop_nonmedical <- 0
    }
  }
  
  if (end_state == "fail") {
    hosp_visit <- rbinom(n = 1, size = 1, p = p_hosp)
    if (hosp_visit == 1) {
      if (init_state == "healthy" | init_state == "diabetes") { 
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_screening + c_laboratory + c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_fail, dr, cycle) 
      } else {
        provider <- sample(c("HC", "CPA1", "CPA2", "CPA3"), size = 1, prob = p_provider, replace = T)
        c_hospitalization <- hospitalization_costs[provider]
        cost_diagnosis <- discount_costs(c_hospitalization, dr, cycle) 
        cost_complications <- discount_costs(c_fail, dr, cycle) 
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
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
      } else {
        oop_medical <- cost_total
        oop_nonmedical <- discount_costs(oop_transport_hosp, dr, cycle) 
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
