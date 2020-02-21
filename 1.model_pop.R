add_hef_to_pop <- function(pop_df, pop_x, hef_quantile){
  stopifnot(hef_quantile <= BASE_HEF_QUANTILE)
  pop_df %>% mutate(
    hef_eligible = if (hef_quantile == BASE_HEF_QUANTILE) 1 else income < quantile(income, hef_quantile/BASE_HEF_QUANTILE)[[1]],
    hef = hef_eligible & ( pop_x[, dim(pop_x)[2]] < p_hef_enrollment ),
    hef_utilization = ifelse(hef, pop_x[, dim(pop_x)[2] - 1] < p_hef_utilization, NA)
  )
}
