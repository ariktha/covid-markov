library(tidyverse)
library(broom)
library(data.table)
library(msm)
library(parallel)
library(naniar)
library(purrr)
library(glue)

rm(list = ls())

## Set options
options(readr.show_col_types = FALSE)

n.cores = 6

# Defining transition types and corresponding colors
