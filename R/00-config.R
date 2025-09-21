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


rm(list = ls())

## Set options
options(readr.show_col_types = FALSE)

n.cores = 10

# Sections to run

fit_no_cov_models <- FALSE
comp_no_cov_models <- FALSE

fit_cov_models <- FALSE
comp_cov_models <- FALSE

fit_spline_models <- FALSE
comp_spline_models <- FALSE
hr_spline_models <- FALSE

# Config for results compilation

config_core <- list(
  prevalence = list(skip = TRUE),
  residuals = list(skip = TRUE)
)

config_core_cov <- list(
  prevalence = list(skip = TRUE),
  residuals = list(skip = TRUE),
  comparison = list(include_across = FALSE)
)

# Defining transition types and corresponding colors

# Define colors for consistency
trend_colors <- c(
  "Worse" = "indianred2",
  "Better" = "darkseagreen4",
  "Death" = "indianred4",
  "Recovery" = "darkseagreen2",
  "Self-transition" = "goldenrod",
  "Other" = "gray50"
)
