
# Model helpers -----------------------------------------------------------

## States
stateNames <- c("healthy", "diabetes", "notx", "oad", "ins", "nephro", "retino", "neuro", "ang", "pvd", "mi", "stroke", "fail", "dm_death", "other_death")

## Distribution of Cambodian population
initialStates_KHR <- function(age, sex, p_diet, p_oad) {
  prevalence <- ifelse(sex == "Male", 0.057, 0.061) # WHO Fact Sheet, Cambodia
  prev_nephro <- ifelse(age > 52, 0.012, 0) # Thomas et al., 2014
  prev_retino <- ifelse(age > 52, 0.012, 0) # Sogbesan and Yutho, 2000
  prev_neuro <- ifelse(age > 52, 0.05, 0) # Lim et al., 2002
  p_oad <- 0.224
  p_ins <- 0.017
  p_diet <- (1 - p_oad - p_ins)
  init_states_KHR <- c((1 - prevalence * (1 + p_diet + p_oad + p_ins + prev_nephro + prev_retino + prev_neuro)), 
                       prevalence, prevalence * p_diet, prevalence * p_oad, prevalence * p_ins, 
                       prevalence * prev_nephro, prevalence * prev_retino, prevalence * prev_neuro, 
                       0, 0, 0, 0, 0, 0, 0) # CVD-related data (among diabetics) not available
  return(init_states_KHR)
}
