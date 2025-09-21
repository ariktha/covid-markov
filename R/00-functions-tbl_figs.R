#' Extract transition trend categories for plotting
#' @param patient_data Patient data
#' @return Data frame with transition trend categories
add_transition_trends <- function(patient_data) {
  
  trend_types <- cross_join(
    tibble(from = c("M", "M1", "M2", "M3", "MS", "S", "S1", "S2")),
    tibble(to = c("M", "M1", "M2", "M3", "MS", "S", "S1", "S2", "R", "D"))
  ) %>%
    mutate(
      from_type = substring(from, 1, 1),
      to_type = substring(to, 1, 1)
    ) %>%
    mutate(
      trend = case_when(
        from == to ~ "Self-transition",
        to == "R" ~ "Recovery",
        to == "D" ~ "Death",
        from_type == "M" & to_type == "M" & 
          as.numeric(sub("M", "", to)) > as.numeric(sub("M", "", from)) ~ "Worse",
        from_type == "M" & to_type == "M" & 
          as.numeric(sub("M", "", to)) < as.numeric(sub("M", "", from)) ~ "Better",
        from_type == "M" & to_type == "S" ~ "Worse",
        from_type == "S" & to_type == "M" ~ "Better",
        from_type == "S" & to_type == "S" & 
          as.numeric(sub("S", "", to)) > as.numeric(sub("S", "", from)) ~ "Worse",
        from_type == "S" & to_type == "S" & 
          as.numeric(sub("S", "", to)) < as.numeric(sub("S", "", from)) ~ "Better",
        from_type == "M" & to == "MS" ~ "Other",
        from == "MS" & to_type == "M" ~ "Other",
        TRUE ~ "Other"
      )
    ) %>%
    filter(!is.na(trend)) %>%
    dplyr::select(from, to, trend)
  
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
