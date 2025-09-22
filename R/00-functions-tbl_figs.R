#' Extract transition trend categories for plotting
#' @param patient_data Patient data
#' @return Data frame with transition trend categories
add_transition_trends <- function(patient_data) {
  
  
  
  patient_data <- patient_data %>%
    left_join(trend_types)
  
  return(patient_data)
}

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
