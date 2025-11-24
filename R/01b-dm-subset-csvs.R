library(here)

## Setup

source(here("R", "00-config.R"))

## Read encounter ids and patient ids

enc_ids <- readRDS(here("data", "preprocess", "encounter_ids.rds"))
pat_ids <- readRDS(here("data", "preprocess", "patient_ids.rds"))

# Reading and saving data from the files that are specific to patients or encounters

# Define the file paths and filtering parameters
file_paths <- list.files(here("data", "raw"), pattern = "*.csv", full.names = TRUE) 
primary_filter_column <- "deid_pat_id"  
primary_filter_values <- pat_ids

secondary_filter_column <- "deid_enc_id"  
secondary_filter_values <- enc_ids

# Function to filter a file
filter_file <- function(file_path, primary_filter_column, primary_filter_values, 
                        secondary_filter_column, secondary_filter_values) {
  message("Processing file: ", file_path)  # Log progress
  
  # Read the file with pipe as the delimiter
  full_data <- tryCatch({
    fread(
      file_path,
      sep = "|",
      fill = TRUE,            # Fill missing fields with NA
      quote = "",             # Disable quote handling
      showProgress = TRUE     # Could suppress progress bar for cleaner output
    )
  }, error = function(e) {
    message("Error reading file: ", file_path, "\n", e)
    return(NULL)
  })
  
  if (is.null(full_data)) return(NULL)  # Skip file if read failed
  
  # Determine which filter column to use
  if (primary_filter_column %in% colnames(full_data)) {
    filter_column <- primary_filter_column
    filter_values <- primary_filter_values
  } else if (secondary_filter_column %in% colnames(full_data)) {
    filter_column <- secondary_filter_column
    filter_values <- secondary_filter_values
  } else {
    message("No relevant filter column found in file: ", file_path)
    return(NULL)  # Skip file if neither column is present
  }
  
  # Filter the data based on the chosen filter column and values
  filtered_data <- full_data[full_data[[filter_column]] %in% filter_values]
  
  return(filtered_data)
}

# Process all files and store each filtered dataset in a named list
filtered_files <- lapply(file_paths, function(file) {
  tryCatch(
    filter_file(
      file, 
      primary_filter_column, primary_filter_values, 
      secondary_filter_column, secondary_filter_values
    ),
    error = function(e) {
      message("Error processing file: ", file, "\n", e)
      return(NULL)
    }
  )
})

# Assign names to the list based on the file names
names(filtered_files) <- basename(file_paths)

# Save each filtered dataset to its own RDS file
output_dir <- here("data", "preprocess")
dir.create(output_dir, showWarnings = FALSE)  # Ensure output directory exists

for (file_name in names(filtered_files)) {
  filtered_data <- filtered_files[[file_name]]
  if (!is.null(filtered_data)) {
    # Construct the output RDS file path
    var_name <- tools::file_path_sans_ext(basename(file_name))
    var_name <- substr(var_name, 1, nchar(var_name) - 8)
    
    output_file <- file.path(output_dir, paste0(var_name, ".rds"))
    
    # Save the filtered dataset as an RDS file
    saveRDS(filtered_data, output_file)
    
    message("Saved filtered data to: ", output_file)  # Log saved file
  }
}

# Reading and saving data from the files that aren't specific to patients or encounters

o2_dev_names <- read_csv(here("data", "raw", "others", "o2_dev_names.csv"))
diagnosis_classification <- read_csv(here("data", "raw", "others", "primary_diagnoses_SB.csv"))

saveRDS(o2_dev_names, here("data", "preprocess", "o2_dev_names.rds"))
saveRDS(diagnosis_classification, here("data", "preprocess", "primary_diagnoses_SB.rds"))

