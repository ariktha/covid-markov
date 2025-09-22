rm(list = ls())

# Load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(broom)
library(data.table)
library(msm)
library(parallel)
library(naniar)
library(purrr)
library(glue)
library(splines)
library(stringr)
library(future)
library(future.apply)
library(R.utils)

# Set options and parameters -----------------------------------------------

options(readr.show_col_types = FALSE)

n.cores = 10
# setOption("future.globals.maxSize", 891289600) # 850MB

# Covariate specifications

continuous_covariates <- c("age", "cci_score", "BMI")

key_covariates <- c("age", "sex", "race", "ethnicity", "language", 
                    "insurance_type", "smoking", "BMI", "bmi_cat", 
                    "COVID_vax", "cci_score", "chf", "cci_cat", 
                    "copd", "dnr_on_admit")

key_covariates_labels <- c("Age", "Sex", "Race", "Ethnicity", "Language", 
                           "Insurance type", "Smoking status", "BMI", "BMI category", 
                           "COVID-19 vaccination", "Charlson comorbidity index", "CHF", "CCI category", 
                           "COPD", "DNR on admission")

# Sections to run ----------------------------------------------------------

run_data_setup <- FALSE

fit_no_cov_models <- FALSE
comp_no_cov_models <- FALSE

## Check Markov assumption ----

fit_markov_models <- FALSE
comp_markov_models <- FALSE

## Covariate effects ----------

fit_univar_models <- FALSE
comp_univar_models <- FALSE

fit_spline_models <- TRUE
comp_spline_models <- TRUE
hr_spline_models <- TRUE

fit_multivar_models <- TRUE
comp_multivar_models <- TRUE

fit_trans_models <- FALSE
comp_trans_models <- FALSE

## Time-homogeneity ----------

fit_time_vary_models <- FALSE
comp_time_vary_models <- FALSE

run_long_stay_analysis <- FALSE

# Config for results compilation --------------------------------

config_core <- list(
  prevalence = list(skip = TRUE),
  residuals = list(skip = TRUE)
)

config_core_cov <- list(
  prevalence = list(skip = TRUE),
  residuals = list(skip = TRUE),
  comparison = list(include_across = FALSE)
)

# Transition types and colors ------------

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

# Define colors for consistency
trend_colors <- c(
  "Worse" = "indianred2",
  "Better" = "darkseagreen4",
  "Death" = "indianred4",
  "Recovery" = "darkseagreen2",
  "Self-transition" = "goldenrod",
  "Other" = "gray50"
)
