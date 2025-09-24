<<<<<<< HEAD
=======
#' Extract transition trend categories for plotting
#' @param patient_data Patient data
#' @return Data frame with transition trend categories
add_transition_trends <- function(patient_data) {
  
  
  
  patient_data <- patient_data %>%
    left_join(trend_types)
  
  return(patient_data)
}

>>>>>>> c6e20c8b285831bacedd41ae31ed081e6d34a949
#' Format confidence intervals
format_ci <- function(estimate, lower, upper, digits = 2) {
  paste0(
    round(estimate, digits), 
    " (95% CI: ", 
    round(lower, digits), 
    "—", 
    round(upper, digits), 
    ")"
  )
}

add_trend_type <- function(df, from_name = "from_state", to_name = "to_state") {
  df %>%
    left_join(trend_types, by = join_by(from_name == "from", to_name == "to"))
}

add_covariate_labels <- function(df, cov_col_name = "covariate") {
  cov_sub <- unique(df[[cov_col_name]])
  df %>%
    mutate(covariate = factor())
}