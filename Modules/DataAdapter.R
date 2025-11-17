# Data Adapter for DROMA Web Integration
# This file provides compatibility functions to bridge between old and new data structures

# Global variables for DROMA objects
.dromadata <- reactiveValues(
  connection = NULL,
  multiDS = NULL,
  initialized = FALSE
)

# Initialize DROMA data connection
initializeDROMAData <- function() {
  if (!.dromadata$initialized) {
    tryCatch({
      # Load database path from config.yml
      if (requireNamespace("config", quietly = TRUE)) {
        config_data <- config::get()
        db_path <- config_data$database_path
      } else {
        # Fallback to default path if config package not available
        db_path <- "data/droma.sqlite"
      }
      
      # Normalize path (handle relative paths)
      if (!file.exists(db_path)) {
        # Try relative to current working directory
        db_path <- normalizePath(file.path(getwd(), db_path), mustWork = FALSE)
      }
      
      # Check if database exists
      if (!file.exists(db_path)) {
        stop("DROMA database not found at: ", db_path,
             "\nPlease check config.yml or ensure droma.sqlite is in the data/ directory")
      }

      # Connect to database
      .dromadata$connection <- DROMA.Set::connectDROMADatabase(db_path)

      # Create MultiDromaSet object
      .dromadata$multiDS <- DROMA.Set::createMultiDromaSetFromAllProjects(db_path)

      .dromadata$initialized <- TRUE
      return(TRUE)

    }, error = function(e) {
      message("Error initializing DROMA data: ", e$message)
      return(FALSE)
    })
  }
  return(.dromadata$initialized)
}

# Get available projects
getAvailableProjects <- function() {
  if (.dromadata$initialized) {
    return(.dromadata$multiDS@name)
  }
  return(character(0))
}

# Get available molecular profiles
getAvailableMolecularProfiles <- function() {
  if (.dromadata$initialized) {
    # Get from first project as example
    if (length(.dromadata$multiDS@DromaSets) > 0) {
      firstDS <- .dromadata$multiDS@DromaSets[[1]]
      return(names(firstDS@molecularProfiles))
    }
  }
  return(character(0))
}

# Get available treatment responses
getAvailableTreatmentResponses <- function() {
  if (.dromadata$initialized) {
    # Get from first project as example
    if (length(.dromadata$multiDS@DromaSets) > 0) {
      firstDS <- .dromadata$multiDS@DromaSets[[1]]
      return(names(firstDS@treatmentResponse))
    }
  }
  return(character(0))
}

# Get drug data using DROMA_R
getDrugData <- function(drug_name, data_type_filter = NULL, tumor_type_filter = NULL) {
  if (.dromadata$initialized) {
    tryCatch({
      result <- DROMA.R::getDrugSensitivityData(
        dromaset_object = .dromadata$multiDS,
        select_drugs = drug_name,
        data_type = if(is.null(data_type_filter)) "all" else data_type_filter,
        tumor_type = if(is.null(tumor_type_filter)) "all" else tumor_type_filter
      )
      return(result)
    }, error = function(e) {
      message("Error getting drug data: ", e$message)
      return(NULL)
    })
  }
  return(NULL)
}

# Get molecular data
getMolecularData <- function(feature_type, feature_name, data_type_filter = NULL, tumor_type_filter = NULL) {
  if (.dromadata$initialized) {
    tryCatch({
      # Check object type and use appropriate function
      if (inherits(.dromadata$multiDS, "MultiDromaSet")) {
        # Use multi-project function for MultiDromaSet
        data <- DROMA.Set::loadMultiProjectMolecularProfiles(
          .dromadata$multiDS,
          feature_type = feature_type,
          select_features = feature_name,
          data_type = if(is.null(data_type_filter)) "all" else data_type_filter,
          tumor_type = if(is.null(tumor_type_filter)) "all" else tumor_type_filter
        )
      } else if (inherits(.dromadata$multiDS, "DromaSet")) {
        # Use single-project function for DromaSet
        data <- DROMA.Set::loadMolecularProfiles(
          .dromadata$multiDS,
          feature_type = feature_type,
          select_features = feature_name,
          return_data = TRUE,
          data_type = if(is.null(data_type_filter)) "all" else data_type_filter,
          tumor_type = if(is.null(tumor_type_filter)) "all" else tumor_type_filter
        )
        # Convert to list format for consistency
        if (!is.list(data)) {
          data <- list(data)
          names(data) <- .dromadata$multiDS@name
        }
      } else {
        stop("Unknown object type: ", class(.dromadata$multiDS))
      }
      return(data)
    }, error = function(e) {
      message("Error getting molecular data: ", e$message)
      return(NULL)
    })
  }
  return(NULL)
}

# Analyze drug-omic pair
analyzeDrugOmicPairWrapper <- function(drug_name, omic_type, omic_name,
                                      data_type_filter = NULL, tumor_type_filter = NULL,
                                      merge_studies = FALSE, output_option = "ggplot") {
  if (.dromadata$initialized) {
    tryCatch({
      # Get z-score setting from global state
      zscore_enabled <- TRUE
      if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
        zscore_enabled <- isTRUE(base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$enabled)
      }
      
      result <- DROMA.R::analyzeDrugOmicPair(
        dromaset_object = .dromadata$multiDS,
        feature_type = omic_type,
        select_features = omic_name,
        select_drugs = drug_name,
        data_type = if(is.null(data_type_filter)) "all" else data_type_filter,
        tumor_type = if(is.null(tumor_type_filter)) "all" else tumor_type_filter,
        merged_enabled = merge_studies,
        zscore = zscore_enabled
      )
      return(result)
    }, error = function(e) {
      message("Error in drug-omic analysis: ", e$message)
      return(NULL)
    })
  }
  return(NULL)
}

# Batch feature analysis
batchFindSignificantFeaturesWrapper <- function(feature1_type, feature2_type,
                                               feature1_names = NULL, feature2_names = NULL,
                                               data_type_filter = NULL, tumor_type_filter = NULL,
                                               cores = NULL) {
  if (.dromadata$initialized) {
    tryCatch({
      if (is.null(cores)) {
        cores <- parallel::detectCores() - 1
      }

      result <- DROMA.R::batchFindSignificantFeatures(
        dromaset_object = .dromadata$multiDS,
        feature1_type = feature1_type,
        feature1_name = feature1_names,
        feature2_type = feature2_type,
        feature2_name = feature2_names,
        data_type = if(is.null(data_type_filter)) "all" else data_type_filter,
        tumor_type = if(is.null(tumor_type_filter)) "all" else tumor_type_filter,
        cores = cores
      )
      return(result)
    }, error = function(e) {
      message("Error in batch feature analysis: ", e$message)
      return(NULL)
    })
  }
  return(NULL)
}

# Convert old data format to new format (for compatibility)
convertToLegacyFormat <- function(data) {
  if (is.null(data)) {
    return(NULL)
  }

  # If it's already a data frame, return as is
  if (is.data.frame(data)) {
    return(data)
  }

  # Convert S4 objects to data frame if needed
  if (is(data, "DromaSet") || is(data, "MultiDromaSet")) {
    # Extract relevant information and convert to list/data frame
    # This is a simplified conversion - adjust based on actual needs
    return(list(data = data))
  }

  return(data)
}

# Get sample metadata
getSampleMetadata <- function() {
  if (.dromadata$initialized) {
    return(.dromadata$multiDS@sampleMetadata)
  }
  return(data.frame())
}

# Get treatment metadata
getTreatmentMetadata <- function() {
  if (.dromadata$initialized) {
    return(.dromadata$multiDS@treatmentMetadata)
  }
  return(data.frame())
}

# Z-score normalization wrapper
zscoreNormalizeWrapper <- function(data, by = "row") {
  # Use DROMA_Set's zscore normalization or custom implementation
  if (is.data.frame(data)) {
    if (by == "row") {
      return(t(scale(t(data))))
    } else {
      return(scale(data))
    }
  }
  return(data)
}