# Load Data Module - Using DROMA packages
# This module replaces the old Rda file loading with DROMA_Set database connection

# Source required modules and functions
source("Modules/DataAdapter.R")
source("Functions/FeatureLoading.R")

# Initialize DROMA data connection
# The initialization happens once and data is available globally
init_success <- initializeDROMAData()

# Check if initialization was successful
if (!init_success) {
  stop("Failed to initialize DROMA data. Please check database configuration.")
}

# Make data available globally
# Create reactive values for data access
loadDataValues <- reactiveValues(
  initialized = TRUE
)

# Function to load molecular profiles (legacy compatibility)
loadMolecularProfiles <- function(type, name = NULL, dataType = NULL, tumorType = NULL) {
  return(getMolecularData(type, name, dataType, tumorType))
}

# Function to load drug data (legacy compatibility)
loadDrugData <- function(drugName, dataType = NULL, tumorType = NULL) {
  return(getDrugData(drugName, dataType, tumorType))
}

# Get available data types
getAvailableDataTypes <- function() {
  return(names(getSampleTypes()))
}

# Get available tumor types
getAvailableTumorTypes <- function() {
  if (.dromadata$initialized) {
    meta <- getSampleMetadata()
    if ("tumor_type" %in% colnames(meta)) {
      return(unique(meta$tumor_type[!is.na(meta$tumor_type)]))
    }
  }
  return(character(0))
}

# Get available drugs
getAvailableDrugs <- function() {
  if (.dromadata$initialized) {
    meta <- getTreatmentMetadata()
    if (nrow(meta) > 0) {
      # Try different possible column names
      drug_col <- NULL
      if ("treatment_name" %in% colnames(meta)) {
        drug_col <- "treatment_name"
      } else if ("DrugName" %in% colnames(meta)) {
        drug_col <- "DrugName"
      } else if ("drug_name" %in% colnames(meta)) {
        drug_col <- "drug_name"
      }
      
      if (!is.null(drug_col)) {
        drug_names <- unique(meta[[drug_col]])
        drug_names <- drug_names[drug_names != "" & !is.na(drug_names)]
        if (length(drug_names) > 0) {
          return(sort(drug_names))
        }
      }
    }
  }
  return(character(0))
}

# Legacy compatibility - these were the old variable names
# Now they are dynamically loaded from DROMA database
getOmicsData <- function() {
  warning("getOmicsData() is deprecated. Please use getMolecularData() instead.")
  return(NULL)
}

getDrugSensitivityData <- function() {
  warning("getDrugSensitivityData() is deprecated. Please use getDrugData() instead.")
  return(NULL)
}

# For backward compatibility with existing code
# These functions provide the same interface as before but use DROMA backend
loadOldData <- function() {
  # This function is kept for compatibility
  # The actual data loading is now done via DROMA packages
  return(invisible(TRUE))
}