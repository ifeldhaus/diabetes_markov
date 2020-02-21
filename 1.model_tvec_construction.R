

## 5. Transition probabilities ---------------------------------------------

#  Complications; Various sources see Flessa & Zembok 2014 -- total: 0.1189
p_nephro <- 0.01
p_retino <- 0.02124
p_nephro_to_retino <- p_retino * 5.01 # Jeng et al., 2016 -- adjusted HRs for NPDR and PDR were 5.01 (95% confidence interval (CI) = 4.68–5.37) and 9.7 (95% CI = 8.15–11.5), respectively
p_neuro <- 0.04658
p_ang <- 0.00668
p_pvd <- 0.00847 
p_mi <- 0.01735 
p_stroke <- 0.00529
p_fail <- 0.00329
e_nephro_to_cvd <- 1.83 # RR, Gerstein et al., 2001
p_pvd_to_mi <- 0.17 # 0.34 2-year incidence; Rossi et al., 2002
p_pvd_to_ang <- 0.05 # 0.10 2-year incidence; Rossi et al., 2002
e_ang_to_stroke <- 1.14 # HR, Eisen et al., 2016 
e_ang_to_fail <- 1.17 # HR, Eisen et al., 2016
# Note: Other effects non-significant (see Eisen et al., 2016, Grenon et al., 2013)

#  Complication-related death
p_mi_death <- 0.7068 # see Flessa & Zembok 2014; Turner et al. 1998
p_stroke_death <- 0.38136 # see Flessa & Zembok 2014; Turner et al. 1998
p_nephro_death <- 0.311 # Afkarian et al. 2013
p_retino_death <- 0 # Frith & Loprinzi 2018; had to change because it increases beyond 1
p_neuro_death <- 0 
p_ang_death <- 0 
p_pvd_death <- 0.016 * 0.099 # combined probability of (1) amputation and (2) subsequent death within one year; Hoffstad et al. 2015
p_fail_death <- 0.327 # Bell et al. 2019; Bertoni et al., 2004

#  Treatment
p_oad <- 0.224 # Taniguchi et al., 2016; differs significantly from King et al., 2005
p_ins <- 0.017 # Van Olmen et al. 2015
p_diet <- (1 - p_oad - p_ins)
p_oad_to_ins <- 0.04 # Ringborg et al. 2010

## Effects of treatments on transition probabilities (relative risks)
#  Diet (no treatment)
e_diet <- 0.90 # Assumption

#  OAD -- focus on metaformin (and glyburide/sulphonylureas) as most likely candidates for prescription in LMIC settings; probability of incidence (or mortality but used as incidence in transition matrix)
e_oad <- 0.68 # Turner et al. 1998 (UKPDS): Any diabetes-related endpoint -- 0.68 (0.53, 0.87)
e_oad_mi <- 0.61 # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS) (39% risk reduction; 0.61 (0.41-0.89))
e_oad_stroke <- 0.59 # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS) -- for metaformin: 0.59 (0.29-1.18)
e_oad_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 74% risk reduction
e_oad_retino <- e_oad # Turner et al. 1998 (UKPDS)
e_oad_neuro <- 0.94 # Juster-Switlyk et al. 2016: Type 1 reduces 60-70%, Type 2 reduces 5-7%
e_oad_ang <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_pvd <- 0.74 # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS) -- 0.74 (0.26-2.09)
e_oad_fail <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_dm_death <- 0.58 # Turner et al. 1998 (UKPDS): Diabetes-related death -- 0.58 (0.37-0.91)
e_oad_all_cause_death <- 0.64 # Turner et al. 1998 (UKPDS): all-cause mortality -- 0.64 (0.45-0.91)

#  Insulin -- same as above
e_ins <- 0.68 # Turner et al. 1998 (UKPDS): Any diabetes-related endpoint -- 0.68 (0.53, 0.87)
e_ins_mi <- 0.39 # Turner et al. 1998 (UKPDS) - did not differ from metaformin group; Herman et al. 2017; ranges to increased risk for CV events - see Table 1 (uses last registry entry for MACE)
e_ins_stroke <- 0.59 # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS): 0.59 (0.29-1.18)
e_ins_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 74% risk reduction
e_ins_retino <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_neuro <- 0.94 # Juster-Switlyk et al. 2016: Type 1 reduces 60-70%, Type 2 reduces 5-7%
e_ins_ang <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_pvd <- 0.74 # Turner et al. 1998 (UKPDS) -- 0.74 (0.26-2.09)
e_ins_fail <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_dm_death <- 0.58 # Turner et al. 1998 (UKPDS): Diabetes-related death -- 0.58 (0.37-0.91)
e_ins_all_cause_death <- 0.64 # Turner et al. 1998 (UKPDS): all-cause mortality -- 0.64 (0.45-0.91)

# General probablities
p_ad <- 0.125  # Normal
p_diag_base  <- 0.3696  # Normal

# Effects of strategies
e_screen_on_diag <- 1.5
e_screen_on_util <- 2 # RR
e_comp_on_util <- 2 # RR
p_ad_if_tx <- 0.40  # Increased adherence


# Main tvec function ------------------------------------------------------

build_tvec <- function(current_state, strategy, hef, age, sex, p_check = 0.0001){

  ## Set relevant probabilities
  is_strategy_screen <- strategy == "screen_only" | strategy == "screen_tx" | strategy == "screen_tx_comp"
  is_strategy_tx <- !(strategy == "base" | strategy == "screen_only" | strategy == "comp_only")
  
  #  Diagnosis; Source: Flessa & Zembok 2014
  p_diag <- p_diag_base * ifelse(is_strategy_screen & hef, e_screen_on_diag, 1)
  
  p_oad_ad <- ifelse(is_strategy_tx & hef, p_ad_if_tx, p_ad)
  p_ins_ad <- ifelse(is_strategy_tx & hef, p_ad_if_tx, p_ad)
  
  # Assign parameters based on individual characteristics
  p_other_death <- all_cause_mortality[ all_cause_mortality$age == age & all_cause_mortality$sex == sex, "mr"]
  p_dm_death <- dm_mortality[ dm_mortality$age == age & dm_mortality$sex_name == sex, "val"]
  p_dm <- dm_incidence[dm_incidence$age == age & dm_incidence$sex_name == sex, "val"]
  
  ## Initialize the transition matrix
  tvec <- rep(NA, 15)
  names(tvec) <- stateNames
  
  
  if (current_state == 'healthy') {
    tvec["diabetes"] <- p_dm * (1 - p_diag)
    tvec["notx"] <- p_dm * p_diag * p_diet
    tvec["oad"] <- p_dm * p_diag * p_oad
    tvec["ins"] <- p_dm * p_diag * p_ins
    tvec["nephro"] <- p_dm * p_nephro
    tvec["retino"] <- p_dm * p_retino
    tvec["neuro"] <- p_dm * p_neuro
    tvec["ang"] <- p_dm * p_ang
    tvec["pvd"] <- p_dm * p_pvd
    tvec["mi"] <- p_dm * p_mi
    tvec["stroke"] <- p_dm * p_stroke
    tvec["fail"] <- p_dm * p_fail
    tvec["dm_death"] <- p_dm * p_dm_death
    tvec["other_death"] <- p_other_death
    tvec["healthy"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'diabetes') {
    tvec["healthy"] <- 0
    tvec["notx"] <- p_diag * p_diet
    tvec["oad"] <- p_diag * p_oad
    tvec["ins"] <- p_diag * p_ins
    tvec["nephro"] <- p_nephro 
    tvec["retino"] <- p_retino 
    tvec["neuro"] <- p_neuro 
    tvec["ang"] <- p_ang 
    tvec["pvd"] <- p_pvd 
    tvec["mi"] <- p_mi 
    tvec["stroke"] <- p_stroke 
    tvec["fail"] <- p_fail 
    tvec["dm_death"] <- p_dm_death
    tvec["other_death"] <- p_other_death
    tvec["diabetes"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'notx') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["oad"] <- p_oad 
    tvec["ins"] <- p_ins
    tvec["nephro"] <- p_nephro * e_diet
    tvec["retino"] <- p_retino * e_diet
    tvec["neuro"] <- p_neuro * e_diet
    tvec["ang"] <- p_ang * e_diet
    tvec["pvd"] <- p_pvd * e_diet
    tvec["mi"] <- p_mi * e_diet
    tvec["stroke"] <- p_stroke * e_diet
    tvec["fail"] <- p_fail * e_diet
    tvec["dm_death"] <- p_dm_death
    tvec["other_death"] <- p_other_death
    tvec["notx"] <- 1 - sum(tvec[], na.rm = T)
  }
  
  if (current_state == 'oad') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["ins"] <- p_oad_to_ins
    tvec["nephro"] <- p_nephro * (1 - p_oad_ad * e_oad_nephro)
    tvec["retino"] <- p_retino * (1 - p_oad_ad * e_oad_retino)
    tvec["neuro"] <- p_neuro * (1 - p_oad_ad * e_oad_neuro)
    tvec["ang"] <- p_ang * (1 - p_oad_ad * e_oad_ang)
    tvec["pvd"] <- p_pvd * (1 - p_oad_ad * e_oad_pvd)
    tvec["mi"] <- p_mi * (1 - p_oad_ad * e_oad_mi)
    tvec["stroke"] <- p_stroke * (1 - p_oad_ad * e_oad_stroke)
    tvec["fail"] <- p_fail * (1 - p_oad_ad * e_oad_fail)
    tvec["dm_death"] <- p_dm_death * (1 - p_oad_ad * e_oad_dm_death)
    tvec["other_death"] <- p_other_death * (1 - p_oad_ad * e_oad_all_cause_death)
    tvec["oad"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'ins') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["nephro"] <- p_nephro * (1 - p_ins_ad * e_ins_nephro)
    tvec["retino"] <- p_retino * (1 - p_ins_ad * e_ins_retino)
    tvec["neuro"] <- p_neuro * (1 - p_ins_ad * e_ins_neuro)
    tvec["ang"] <- p_ang * (1 - p_ins_ad * e_ins_ang)
    tvec["pvd"] <- p_pvd * (1 - p_ins_ad * e_ins_pvd)
    tvec["mi"] <- p_mi * (1 - p_ins_ad * e_ins_mi)
    tvec["stroke"] <- p_stroke * (1 - p_ins_ad * e_ins_stroke)
    tvec["fail"] <- p_fail * (1 - p_ins_ad * e_ins_fail)
    tvec["dm_death"] <- p_dm_death * (1 - p_ins_ad * e_ins_dm_death)
    tvec["other_death"] <- p_other_death * (1 - p_oad_ad * e_ins_all_cause_death)
    tvec["ins"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'nephro') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["retino"] <- p_nephro_to_retino
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_ang * e_nephro_to_cvd
    tvec["pvd"] <- p_pvd * e_nephro_to_cvd
    tvec["mi"] <- p_mi * e_nephro_to_cvd
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail * e_nephro_to_cvd
    tvec["dm_death"] <- p_nephro_death
    tvec["other_death"] <- p_other_death
    tvec["nephro"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'retino') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_ang
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_retino_death
    tvec["other_death"] <- p_other_death
    tvec["retino"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'neuro') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["ang"] <- p_ang
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_neuro_death
    tvec["other_death"] <- p_other_death
    tvec["neuro"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'ang') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["neuro"] <- p_neuro
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke * e_ang_to_stroke
    tvec["fail"] <- p_fail * e_ang_to_fail
    tvec["dm_death"] <- p_ang_death
    tvec["other_death"] <- p_other_death
    tvec["ang"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'pvd') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_pvd_to_ang
    tvec["mi"] <- p_pvd_to_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_pvd_death
    tvec["other_death"] <- p_other_death
    tvec["pvd"] <- 1 - sum(tvec, na.rm = T)
  }
  
  if (current_state == 'mi') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_ang
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_mi_death
    tvec["other_death"] <- p_other_death
    p_mi_diagnosed <- (1 - sum(tvec, na.rm = T))
    tvec["notx"] <- p_mi_diagnosed * p_diet
    tvec["oad"] <- p_mi_diagnosed * p_oad
    tvec["ins"] <- p_mi_diagnosed * p_ins
  }
  
  if (current_state == 'stroke') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_ang
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_mi_death
    tvec["other_death"] <- p_other_death
    p_stroke_diagnosed <- (1 - sum(tvec, na.rm = T))
    tvec["notx"] <- p_stroke_diagnosed * p_diet
    tvec["oad"] <- p_stroke_diagnosed * p_oad
    tvec["ins"] <- p_stroke_diagnosed * p_ins
  }
  
  if (current_state == 'fail') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["nephro"] <- p_nephro
    tvec["retino"] <- p_retino
    tvec["neuro"] <- p_neuro
    tvec["ang"] <- p_ang
    tvec["pvd"] <- p_pvd
    tvec["mi"] <- p_mi
    tvec["stroke"] <- p_stroke
    tvec["fail"] <- p_fail
    tvec["dm_death"] <- p_mi_death
    tvec["other_death"] <- p_other_death
    p_fail_diagnosed <- (1 - sum(tvec, na.rm = T))
    tvec["notx"] <- p_fail_diagnosed * p_diet
    tvec["oad"] <- p_fail_diagnosed * p_oad
    tvec["ins"] <- p_fail_diagnosed * p_ins
  }
  
  if (current_state == 'dm_death') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- 0
    tvec["retino"] <- 0
    tvec["neuro"] <- 0
    tvec["ang"] <- 0
    tvec["pvd"] <- 0
    tvec["mi"] <- 0
    tvec["stroke"] <- 0
    tvec["fail"] <- 0
    tvec["dm_death"] <- 1
    tvec["other_death"] <- 0
  }
  
  if (current_state == 'other_death') {
    tvec["healthy"] <- 0
    tvec["diabetes"] <- 0
    tvec["notx"] <- 0
    tvec["oad"] <- 0
    tvec["ins"] <- 0
    tvec["nephro"] <- 0
    tvec["retino"] <- 0
    tvec["neuro"] <- 0
    tvec["ang"] <- 0
    tvec["pvd"] <- 0
    tvec["mi"] <- 0
    tvec["stroke"] <- 0
    tvec["fail"] <- 0
    tvec["dm_death"] <- 0
    tvec["other_death"] <- 1
  }
  
  return(tvec)
  
}
