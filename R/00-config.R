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
library(patchwork)

# Set options and parameters -----------------------------------------------

options(readr.show_col_types = FALSE)

n.cores = 6
# setOption("future.globals.maxSize", 891289600) # 850MB

time_vec <- c(1, 7, 15, 30)

# Covariate specifications

continuous_covariates <- c("age", "cci_score", "BMI")

key_covariates <- c("age", "age_cat", "sex", "race", "ethnicity", "language", 
                    "insurance_type", "smoking", "BMI", "bmi_cat", 
                    "vax", "cci_score", "chf", "cci_cat", 
                    "copd", "dnr_on_admit", "COVID_tx")

key_covariates_labels <- c("Age", "Age category", "Sex", "Race", "Ethnicity", "Language", 
                           "Insurance type", "Smoking status", "BMI", "BMI category", 
                           "COVID-19 vaccination", "Charlson comorbidity index", "CHF", "CCI category", 
                           "COPD", "DNR on admission", "COVID-19 treatment")

# Sections to run ----------------------------------------------------------

run_data_setup <- TRUE



do_no_cov <- FALSE
fit_no_cov_models <- TRUE
comp_no_cov_models <- TRUE

## Check Markov assumption ----

do_markov <- FALSE
fit_markov_models <- TRUE
comp_markov_models <- TRUE

## Covariate effects ----------

do_univar <- TRUE
univar_timeout_times <- c(1, 5, 10)
fit_univar_models <- TRUE
refit_failed_univar_models <- TRUE
comp_univar_models <- TRUE

# Not implemented
do_multivar <- FALSE
fit_multivar_models <- TRUE
comp_multivar_models <- TRUE

do_constrained_trans <- TRUE
fit_trans_models <- TRUE
comp_trans_models <- TRUE

## Time-homogeneity ----------

do_time_inhomogeneity <- TRUE
fit_time_vary_models <- TRUE
comp_time_vary_models <- TRUE

## Long stay sensitivity analysis ----------

run_long_stay_analysis <- TRUE

# Config for results compilation --------------------------------

config_core <- list(
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


stage_order <- c("4", "5", "6", "7", "8", "9", "10", "11")
state_order <- c("M", "M1", "M2", "M3", "MS", "S", "S1", "S2", "D", "R")

models <- c(
  "base_model" = "Base model",
  "reduced_trans" = "Simplified transitions",
  "hx_sev" = "History of severe",
  "mod_2" = "2 moderate states",
  "mod_3" = "3 moderate states",
  "sev_2" = "2 severe states"
)

# model_colors <- c(
#   "black",
#   
# )



models_df <- tibble(
  model = names(models),
  model_name = models
)

