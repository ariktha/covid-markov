library(tidyverse)
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

rm(list = ls())

## Set options
options(readr.show_col_types = FALSE)

n.cores = 6

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
