# COVID-19 Clinical Progression Analysis

## Project Overview
This project analyzes the clinical progression of COVID-19 patients using multi-state Markov models. The analysis examines different state configurations (4-state, 5-state, 6-state, and 7-state models) to understand disease progression patterns.

## Project Structure

### Data Processing Scripts
- `scripts/00-config.R`: Configuration file with library imports and global options
- `scripts/00-functions.R`: Helper functions for data processing and analysis
- `scripts/01a-dm-filter_enc.R`: Filters and processes encounter data
- `scripts/01b-dm-process-datafiles.R`: Processes raw data files
- `scripts/01c-dm-process-rds-files.R`: Processes RDS format data files

### Analysis Scripts
- `scripts/02a-base-models.R`: Implements base Markov models
- `scripts/02b-check-markov.R`: Validates Markov assumptions
- `scripts/02c-covariate-selection.R`: Performs covariate selection for models

### State Definitions
The models use the following state definitions based on disease severity:
- M = Moderate
- S = Severe  
- R = Recovered
- D = Dead
- MS = Moderate with history of Severe

Model configurations:
- 4-state model: M, S, R, D
- 5-state model: M, MS, S, R, D  
- 6-state model: M1, M2, MS, S, R, D
- 7-state model: M1, M2, M3, MS, S, R, D

## Data Requirements
The analysis requires the following input data:
- Demographics and event data
- Medication administration data
- Flowsheet data
- COVID test results
- Diagnosis codes
- Vaccination records

## Key Features
- Filters encounters based on COVID-19 test positivity
- Processes longitudinal patient data into daily states
- Implements multi-state Markov models
- Tests Markov assumptions
- Performs covariate selection
- Compares models using DRAIC and likelihood ratio tests

## Dependencies
- R packages:
  - tidyverse
  - msm
  - here
  - data.table
  - broom

## Usage
1. Place raw data files in `data/raw/` directory
2. Run scripts in numerical order (01a -> 01b -> 01c -> 02a -> 02b -> 02c)
3. Results are saved in `data/temp/` directory

## Output Files
- `data/temp/base_model_results.RData`: Contains fitted base models and results
- `data/temp/pt_enc_mod_staging.rds`: Patient encounter staging data
- `data/temp/crude_base_rates.rds`: Crude transition rates
- `data/temp/fitted_base_models.rds`: Fitted model objects
- `data/temp/tidied_base_results.rds`: Tidied model results 