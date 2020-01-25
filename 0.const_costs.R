
## Direct medical costs; Source: Flessa & Zembok, 2014; adjusted to 2019 USD
c_screening <- 1.10 # i.e. screening; unit cost FPG-test
c_laboratory <- 1.41 # diagnosis cost
c_oad <- 27.86
c_ins <- 125.80
c_outpatient_cpa3 <- 44.40 # Maximum prices -- Costing of Health Care Services in Three Provinces of Cambodia (2018)
c_outpatient_cpa2 <- 6.28 # Secondary hospital
c_outpatient_cpa1 <- 10.32 # Primary hospital
c_outpatient_hc <- 4.15 
c_hospitalization_cpa3 <- 40.85
c_hospitalization_cpa2 <- 29.52
c_hospitalization_cpa1 <- 59.73
c_hospitalization_hc <- 4.77

outpatient_costs <- c(
  "HC" = c_outpatient_hc,
  "CPA1" = c_outpatient_cpa1,
  "CPA2" = c_outpatient_cpa2,
  "CPA3" = c_outpatient_cpa3
)

n_days <- 5
hospitalization_costs <- c(
  "HC" = c_hospitalization_hc,
  "CPA1" = c_hospitalization_cpa1,
  "CPA2" = c_hospitalization_cpa2,
  "CPA3" = c_hospitalization_cpa3
) * n_days

## Direct medical costs associated with each complication
c_nephro <- 6358 # Hemodialysis, systematic review of LMICs; Int$ 3,424 to Int$ 42,785 overall, but used Sri Lanka numbers: Int$ 5,869–8,804; Mushi et al., 2015 -- adjusted for inflation
c_retino <- 330.1 # Cost of intravitreal injection in Indonesia, Sasongko et al., 2019
c_neuro <- 2.62 * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
c_ang <- 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
c_pvd <- 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
c_mi <- 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation
c_stroke <- 2.62 * 365 # Acetylsalicylic acid - 500 mg cap/tab (aspirin) for a year: 1.93 (1.76, 3.08); WHO/HAI database, MSH 2004 -- adjusted for inflation
c_fail <- 7.07 * 365 # CVD medication (beta blocker) - Atenolol, 50 mg for a year: 5.20 (5.03, 5.55); WHO/HAI database, MSH 2004 -- adjusted for inflation

### Individual expenditures
## Indirect (non-medical) costs
oop_transport_op <- 0.92
oop_transport_hosp <- 11.65

## Discount rate
dr <- 0.03

## Subsistence expenditure
subsistence <- 47.96339 * 12 # 196K KHR (CSES 2017) / 4086.45 KHR per USD (2 Oct 2019); average monthly value per capital, 2017

