### Financial risk protection
## Catastrophic health expenditure (CHE)
che_threshold = c(0.1, 0.25, 0.4)
che <- function(oop_individual, income) {
  che_10 <- ifelse(income < 0, NA, 
                   ifelse(oop_individual > income * che_threshold[1], 1, 0))
  che_25 <- ifelse(income < 0, NA,
                   ifelse(oop_individual > income * che_threshold[2], 1, 0))
  che_40 <- ifelse(income < 0, NA, 
                   ifelse(oop_individual > income * che_threshold[3], 1, 0))
  
  return(list("che10" = che_10, "che25" = che_25, "che40" = che_40))
}

## Poverty cases
poverty_line = 52.21 * 12 # 7,112 KHR per day (* 30 days) / 4086.45 KHR per USD (2 Oct 2019)
pov <- function(oop_individual, income) {
  if (income > poverty_line) {
    ifelse(income - oop_individual <= poverty_line, 1, 0)
  } else {
    NA
  }
}