# Compute DALYs

compute_dalys <- function(end_state, age, sex) {
  if (end_state == "healthy" | end_state == "diabetes") {
    dalys_dm <- 0
    dalys_other <- 0
  }
  
  if (end_state == "notx" | end_state == "oad" | end_state == "ins") {
    dalys_dm <- dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "nephro") {
    dalys_dm <- dw_nephropathy_stage5 + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "retino") {
    dalys_dm <- dw_retinopathy + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "neuro") {
    dalys_dm <- dw_neuropathy + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "ang") {
    dalys_dm <- dw_angina_moderate + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "pvd") {
    dalys_dm <- dw_pvd + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "mi") {
    dalys_dm <- dw_mi + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "stroke") {
    dalys_dm <- dw_stroke_level5 + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "fail") {
    dalys_dm <- dw_failure_severe + dw_uncomplicated
    dalys_other <- 0
  }
  
  if (end_state == "dm_death") {
    dalys_dm <- ifelse(sex == "Female", max(71 - age, 0), max(67 - age, 0)) # Male life expectancy: 67 years; Female life expectancy: 71 years
    dalys_other <- 0
  }
  
  if (end_state == "other_death") {
    dalys_dm <- 0
    dalys_other <- ifelse(sex == "Female", max(71 - age, 0), max(67 - age, 0))
  }
  
  dalys_total <- dalys_dm + dalys_other
  dalys <- list(
    "dalys_dm" = dalys_dm, 
    "dalys_other" = dalys_other,
    "dalys_total" = dalys_total
  )
  
  return(dalys)
}
