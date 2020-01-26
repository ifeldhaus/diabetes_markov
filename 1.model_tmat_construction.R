

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
p_fail_death <- 0.327 # Bell et al. 2003

#  Treatment
p_oad <- 0.224 # Taniguchi et al., 2016; differs significantly from King et al., 2005
p_ins <- 0.017 # Van Olmen et al. 2015
p_diet <- (1 - p_oad - p_ins)
p_oad_to_ins <- 0.04 # Ringborg et al. 2010

## Effects of treatments on transition probabilities (relative risks)
#  Diet (no treatment)
e_diet <- 0.90 # Assumption

#  OAD -- focus on metaformin (and glyburide/sulphonylureas) as most likely candidates for prescription in LMIC settings; probability of incidence (or mortality but used as incidence in transition matrix)
e_oad <- 0.5 # Assumption
e_oad_mi <- 0.39 # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS) [61% risk reduction? -> to check]
e_oad_stroke <- e_oad # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS) -- for metaformin
e_oad_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
e_oad_retino <- e_oad # Turner et al. 1998 (UKPDS)
e_oad_neuro <- 0.95 # Juster-Switlyk et al. 2016
e_oad_ang <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_pvd <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)
e_oad_fail <- e_oad # Chaudhury et al. 2017; Turner et al. 1998 (UKPDS)

#  Insulin -- same as above
e_ins <- 0.5 # Assumption
e_ins_mi <- 0.39 # Turner et al. 1998 (UKPDS) - did not differ from metaformin group; Herman et al. 2017; ranges to increased risk for CV events - see Table 1 (uses last registry entry for MACE)
e_ins_stroke <- e_ins # Avgerinos et al. 2017; Turner et al. 1998 (UKPDS)
e_ins_nephro <- 0.3 # Turner et al. 1998 (UKPDS); extrapolated but may be between 67% and 72% risk reduction
e_ins_retino <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_neuro <- 0.95 # Turner et al. 1998 (UKPDS)
e_ins_ang <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_pvd <- e_ins # Turner et al. 1998 (UKPDS)
e_ins_fail <- e_ins # Turner et al. 1998 (UKPDS)

# General probablities
p_ad <- 0.125  # Normal
p_diag_base  <- 0.3696  # Normal

# Effects of strategies
e_screen_on_diag <- 1.5
e_screen_on_util <- 2 # RR
e_comp_on_util <- 2 # RR
p_ad_if_tx <- 0.40  # Increased adherence


# Main tmat function ------------------------------------------------------

build_tmat <- function(strategy, hef, age, sex, p_check = 0.0001){

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
  tmat <- matrix(NA, nrow = 15, ncol = 15)
  colnames(tmat) <- stateNames
  rownames(tmat) <- stateNames
  
  ## Define all values of tmat
  diag(tmat) <- 0
  
  tmat["healthy", "diabetes"] <- p_dm * (1 - p_diag)
  tmat["healthy", "notx"] <- p_dm * p_diag * p_diet
  tmat["healthy", "oad"] <- p_dm * p_diag * p_oad
  tmat["healthy", "ins"] <- p_dm * p_diag * p_ins
  tmat["healthy", "nephro"] <- p_dm * p_nephro
  tmat["healthy", "retino"] <- p_dm * p_retino
  tmat["healthy", "neuro"] <- p_dm * p_neuro
  tmat["healthy", "ang"] <- p_dm * p_ang
  tmat["healthy", "pvd"] <- p_dm * p_pvd
  tmat["healthy", "mi"] <- p_dm * p_mi
  tmat["healthy", "stroke"] <- p_dm * p_stroke
  tmat["healthy", "fail"] <- p_dm * p_fail
  tmat["healthy", "dm_death"] <- p_dm * p_dm_death
  tmat["healthy", "other_death"] <- p_other_death
  tmat["healthy", "healthy"] <- 1 - sum(tmat["healthy",], na.rm = T)
  
  tmat["diabetes", "healthy"] <- 0
  tmat["diabetes", "notx"] <- p_diag * p_diet
  tmat["diabetes", "oad"] <- p_diag * p_oad
  tmat["diabetes", "ins"] <- p_diag * p_ins
  tmat["diabetes", "nephro"] <- p_nephro 
  tmat["diabetes", "retino"] <- p_retino 
  tmat["diabetes", "neuro"] <- p_neuro 
  tmat["diabetes", "ang"] <- p_ang 
  tmat["diabetes", "pvd"] <- p_pvd 
  tmat["diabetes", "mi"] <- p_mi 
  tmat["diabetes", "stroke"] <- p_stroke 
  tmat["diabetes", "fail"] <- p_fail 
  tmat["diabetes", "dm_death"] <- p_dm_death
  tmat["diabetes", "other_death"] <- p_other_death
  tmat["diabetes", "diabetes"] <- 1 - sum(tmat["diabetes",], na.rm = T)
  
  tmat["notx", "healthy"] <- 0
  tmat["notx", "diabetes"] <- 0
  tmat["notx", "oad"] <- p_oad 
  tmat["notx", "ins"] <- p_ins
  tmat["notx", "nephro"] <- p_nephro * e_diet
  tmat["notx", "retino"] <- p_retino * e_diet
  tmat["notx", "neuro"] <- p_neuro * e_diet
  tmat["notx", "ang"] <- p_ang * e_diet
  tmat["notx", "pvd"] <- p_pvd * e_diet
  tmat["notx", "mi"] <- p_mi * e_diet
  tmat["notx", "stroke"] <- p_stroke * e_diet
  tmat["notx", "fail"] <- p_fail * e_diet
  tmat["notx", "dm_death"] <- p_dm_death
  tmat["notx", "other_death"] <- p_other_death
  tmat["notx", "notx"] <- 1 - sum(tmat["notx",], na.rm = T)
  
  tmat["oad", "healthy"] <- 0
  tmat["oad", "diabetes"] <- 0
  tmat["oad", "notx"] <- 0
  tmat["oad", "ins"] <- p_oad_to_ins
  tmat["oad", "nephro"] <- p_nephro * (1 - p_oad_ad * e_oad_nephro)
  tmat["oad", "retino"] <- p_retino * (1 - p_oad_ad * e_oad_retino)
  tmat["oad", "neuro"] <- p_neuro * (1 - p_oad_ad * e_oad_neuro)
  tmat["oad", "ang"] <- p_ang * (1 - p_oad_ad * e_oad_ang)
  tmat["oad", "pvd"] <- p_pvd * (1 - p_oad_ad * e_oad_pvd)
  tmat["oad", "mi"] <- p_mi * (1 - p_oad_ad * e_oad_mi)
  tmat["oad", "stroke"] <- p_stroke * (1 - p_oad_ad * e_oad_stroke)
  tmat["oad", "fail"] <- p_fail * (1 - p_oad_ad * e_oad_fail)
  tmat["oad", "dm_death"] <- p_dm_death
  tmat["oad", "other_death"] <- p_other_death
  tmat["oad", "oad"] <- 1 - sum(tmat["oad",], na.rm = T)
  
  tmat["ins", "healthy"] <- 0
  tmat["ins", "diabetes"] <- 0
  tmat["ins", "notx"] <- 0
  tmat["ins", "oad"] <- 0
  tmat["ins", "nephro"] <- p_nephro * (1 - p_ins_ad * e_ins_nephro)
  tmat["ins", "retino"] <- p_retino * (1 - p_ins_ad * e_ins_retino)
  tmat["ins", "neuro"] <- p_neuro * (1 - p_ins_ad * e_ins_neuro)
  tmat["ins", "ang"] <- p_ang * (1 - p_ins_ad * e_ins_ang)
  tmat["ins", "pvd"] <- p_pvd * (1 - p_ins_ad * e_ins_pvd)
  tmat["ins", "mi"] <- p_mi * (1 - p_ins_ad * e_ins_mi)
  tmat["ins", "stroke"] <- p_stroke * (1 - p_ins_ad * e_ins_stroke)
  tmat["ins", "fail"] <- p_fail * (1 - p_ins_ad * e_ins_fail)
  tmat["ins", "dm_death"] <- p_dm_death
  tmat["ins", "other_death"] <- p_other_death
  tmat["ins", "ins"] <- 1 - sum(tmat["ins",], na.rm = T)
  
  tmat["nephro", "healthy"] <- 0
  tmat["nephro", "diabetes"] <- 0
  tmat["nephro", "notx"] <- 0
  tmat["nephro", "oad"] <- 0
  tmat["nephro", "ins"] <- 0
  tmat["nephro", "retino"] <- p_nephro_to_retino
  tmat["nephro", "neuro"] <- p_neuro
  tmat["nephro", "ang"] <- p_ang * e_nephro_to_cvd
  tmat["nephro", "pvd"] <- p_pvd * e_nephro_to_cvd
  tmat["nephro", "mi"] <- p_mi * e_nephro_to_cvd
  tmat["nephro", "stroke"] <- p_stroke
  tmat["nephro", "fail"] <- p_fail * e_nephro_to_cvd
  tmat["nephro", "dm_death"] <- p_nephro_death
  tmat["nephro", "other_death"] <- p_other_death
  tmat["nephro", "nephro"] <- 1 - sum(tmat["nephro",], na.rm = T)
  
  tmat["retino", "healthy"] <- 0
  tmat["retino", "diabetes"] <- 0
  tmat["retino", "notx"] <- 0
  tmat["retino", "oad"] <- 0
  tmat["retino", "ins"] <- 0
  tmat["retino", "nephro"] <- p_nephro
  tmat["retino", "neuro"] <- p_neuro
  tmat["retino", "ang"] <- p_ang
  tmat["retino", "pvd"] <- p_pvd
  tmat["retino", "mi"] <- p_mi
  tmat["retino", "stroke"] <- p_stroke
  tmat["retino", "fail"] <- p_fail
  tmat["retino", "dm_death"] <- p_retino_death
  tmat["retino", "other_death"] <- p_other_death
  tmat["retino", "retino"] <- 1 - sum(tmat["retino",], na.rm = T)
  
  tmat["neuro", "healthy"] <- 0
  tmat["neuro", "diabetes"] <- 0
  tmat["neuro", "notx"] <- 0
  tmat["neuro", "oad"] <- 0
  tmat["neuro", "ins"] <- 0
  tmat["neuro", "nephro"] <- p_nephro
  tmat["neuro", "retino"] <- p_retino
  tmat["neuro", "ang"] <- p_ang
  tmat["neuro", "pvd"] <- p_pvd
  tmat["neuro", "mi"] <- p_mi
  tmat["neuro", "stroke"] <- p_stroke
  tmat["neuro", "fail"] <- p_fail
  tmat["neuro", "dm_death"] <- p_neuro_death
  tmat["neuro", "other_death"] <- p_other_death
  tmat["neuro", "neuro"] <- 1 - sum(tmat["neuro",], na.rm = T)
  
  tmat["ang", "healthy"] <- 0
  tmat["ang", "diabetes"] <- 0
  tmat["ang", "notx"] <- 0
  tmat["ang", "oad"] <- 0
  tmat["ang", "ins"] <- 0
  tmat["ang", "nephro"] <- p_nephro
  tmat["ang", "retino"] <- p_retino
  tmat["ang", "neuro"] <- p_neuro
  tmat["ang", "pvd"] <- p_pvd
  tmat["ang", "mi"] <- p_mi
  tmat["ang", "stroke"] <- p_stroke * e_ang_to_stroke
  tmat["ang", "fail"] <- p_fail * e_ang_to_fail
  tmat["ang", "dm_death"] <- p_ang_death
  tmat["ang", "other_death"] <- p_other_death
  tmat["ang", "ang"] <- 1 - sum(tmat["ang",], na.rm = T)
  
  tmat["pvd", "healthy"] <- 0
  tmat["pvd", "diabetes"] <- 0
  tmat["pvd", "notx"] <- 0
  tmat["pvd", "oad"] <- 0
  tmat["pvd", "ins"] <- 0
  tmat["pvd", "nephro"] <- p_nephro
  tmat["pvd", "retino"] <- p_retino
  tmat["pvd", "neuro"] <- p_neuro
  tmat["pvd", "ang"] <- p_pvd_to_ang
  tmat["pvd", "mi"] <- p_pvd_to_mi
  tmat["pvd", "stroke"] <- p_stroke
  tmat["pvd", "fail"] <- p_fail
  tmat["pvd", "dm_death"] <- p_pvd_death
  tmat["pvd", "other_death"] <- p_other_death
  tmat["pvd", "pvd"] <- 1 - sum(tmat["pvd",], na.rm = T)
  
  tmat["mi", "healthy"] <- 0
  tmat["mi", "diabetes"] <- 0
  tmat["mi", "nephro"] <- p_nephro
  tmat["mi", "retino"] <- p_retino
  tmat["mi", "neuro"] <- p_neuro
  tmat["mi", "ang"] <- p_ang
  tmat["mi", "pvd"] <- p_pvd
  tmat["mi", "mi"] <- p_mi
  tmat["mi", "stroke"] <- p_stroke
  tmat["mi", "fail"] <- p_fail
  tmat["mi", "dm_death"] <- p_mi_death
  tmat["mi", "other_death"] <- p_other_death
  p_mi_diagnosed <- (1 - sum(tmat["mi", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
  tmat["mi", "notx"] <- p_mi_diagnosed * p_diet
  tmat["mi", "oad"] <- p_mi_diagnosed * p_oad
  tmat["mi", "ins"] <- p_mi_diagnosed * p_ins
  
  tmat["stroke", "healthy"] <- 0
  tmat["stroke", "diabetes"] <- 0
  tmat["stroke", "nephro"] <- p_nephro
  tmat["stroke", "retino"] <- p_retino
  tmat["stroke", "neuro"] <- p_neuro
  tmat["stroke", "ang"] <- p_ang
  tmat["stroke", "pvd"] <- p_pvd
  tmat["stroke", "mi"] <- p_mi
  tmat["stroke", "stroke"] <- p_stroke
  tmat["stroke", "fail"] <- p_fail
  tmat["stroke", "dm_death"] <- p_mi_death
  tmat["stroke", "other_death"] <- p_other_death
  p_stroke_diagnosed <- (1 - sum(tmat["stroke", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
  tmat["stroke", "notx"] <- p_stroke_diagnosed * p_diet
  tmat["stroke", "oad"] <- p_stroke_diagnosed * p_oad
  tmat["stroke", "ins"] <- p_stroke_diagnosed * p_ins
  
  tmat["fail", "healthy"] <- 0
  tmat["fail", "diabetes"] <- 0
  tmat["fail", "nephro"] <- p_nephro
  tmat["fail", "retino"] <- p_retino
  tmat["fail", "neuro"] <- p_neuro
  tmat["fail", "ang"] <- p_ang
  tmat["fail", "pvd"] <- p_pvd
  tmat["fail", "mi"] <- p_mi
  tmat["fail", "stroke"] <- p_stroke
  tmat["fail", "fail"] <- p_fail
  tmat["fail", "dm_death"] <- p_mi_death
  tmat["fail", "other_death"] <- p_other_death
  p_fail_diagnosed <- (1 - sum(tmat["fail", !(colnames(tmat) %in% c("notx", "oad", "ins"))], na.rm = T))
  tmat["fail", "notx"] <- p_fail_diagnosed * p_diet
  tmat["fail", "oad"] <- p_fail_diagnosed * p_oad
  tmat["fail", "ins"] <- p_fail_diagnosed * p_ins
  
  tmat["dm_death", "healthy"] <- 0
  tmat["dm_death", "diabetes"] <- 0
  tmat["dm_death", "notx"] <- 0
  tmat["dm_death", "oad"] <- 0
  tmat["dm_death", "ins"] <- 0
  tmat["dm_death", "nephro"] <- 0
  tmat["dm_death", "retino"] <- 0
  tmat["dm_death", "neuro"] <- 0
  tmat["dm_death", "ang"] <- 0
  tmat["dm_death", "pvd"] <- 0
  tmat["dm_death", "mi"] <- 0
  tmat["dm_death", "stroke"] <- 0
  tmat["dm_death", "fail"] <- 0
  tmat["dm_death", "dm_death"] <- 1
  tmat["dm_death", "other_death"] <- 0
  
  tmat["other_death", "healthy"] <- 0
  tmat["other_death", "diabetes"] <- 0
  tmat["other_death", "notx"] <- 0
  tmat["other_death", "oad"] <- 0
  tmat["other_death", "ins"] <- 0
  tmat["other_death", "nephro"] <- 0
  tmat["other_death", "retino"] <- 0
  tmat["other_death", "neuro"] <- 0
  tmat["other_death", "ang"] <- 0
  tmat["other_death", "pvd"] <- 0
  tmat["other_death", "mi"] <- 0
  tmat["other_death", "stroke"] <- 0
  tmat["other_death", "fail"] <- 0
  tmat["other_death", "dm_death"] <- 0
  tmat["other_death", "other_death"] <- 1
  
  # Check the tmat sometimes
  if (runif(1) < p_check) {
    mc <- new("markovchain",
            transitionMatrix = tmat)
  }
  
  
  return(tmat)
  
}
