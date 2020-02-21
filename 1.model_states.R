
# Model helpers -----------------------------------------------------------

## States
stateNames <- c("healthy", "diabetes", "notx", "oad", "ins", "nephro", "retino", "neuro", "ang", "pvd", "mi", "stroke", "fail", "dm_death", "other_death")

## Distribution of Cambodian population
initialStates_KHR <- function(age, sex) {
  prevalence <- dm_prevalence[dm_prevalence$age == age & dm_prevalence$sex_name == sex, "val"]
  prev_nephro <- ifelse(age > 52, 0.012, 0) # Thomas et al., 2014
  prev_retino <- ifelse(age > 52, 0.012, 0) # Sogbesan and Yutho, 2000
  prev_neuro <- ifelse(age > 52, 0.05, 0) # Lim et al., 2002
  p_oad <- 0.224
  p_ins <- 0.017
  p_diet <- (1 - p_oad - p_ins)
  p_diag_base <- 0.3696
  init_states_KHR <-
    c((1 - (prevalence - (prevalence * p_diag_base * (p_diet + p_oad + p_ins + prev_nephro + prev_retino + prev_neuro))) - (prevalence * p_diag_base * (p_diet + p_oad + p_ins + prev_nephro + prev_retino + prev_neuro))),
    prevalence - (prevalence * p_diag_base * (p_diet + p_oad + p_ins + prev_nephro + prev_retino + prev_neuro)),
    prevalence * p_diag_base * p_diet,
    prevalence * p_diag_base * p_oad,
    prevalence * p_diag_base * p_ins,
    prevalence * p_diag_base * prev_nephro,
    prevalence * p_diag_base * prev_retino,
    prevalence * p_diag_base * prev_neuro,
    0,
    0,
    0,
    0,
    0,
    0,
    0
    ) # CVD-related data (among diabetics) not available
  return(init_states_KHR)
}


draw_state_from_x <- function(state_probs, x) {
  names(state_probs)[x < cumsum(state_probs)][1]
}
