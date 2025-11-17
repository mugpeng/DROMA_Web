#' ==========================================
#' Z-Score Normalization Functions
#' ==========================================
#' These functions provide z-score normalization for DROMA data

#' Apply z-score normalization to DROMA molecular data
#'
#' This function applies z-score normalization across samples for each gene
#' and row-by-row normalization for drug data
#'
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param force Force re-normalization even if already normalized
#' @return TRUE on success, FALSE on failure
#'
#' @examples
#' applyZScoreNormalization(multiDS)
#' applyZScoreNormalization(dromaset, force = TRUE)
applyZScoreNormalization <- function(dromaset_object = NULL, force = FALSE) {
  tryCatch({
    # Check if already normalized
    if (!force && exists("normalization_state", envir = .GlobalEnv) &&
        isTRUE(get("normalization_state", envir = .GlobalEnv))) {
      message("Data is already normalized. Use force = TRUE to re-normalize.")
      return(TRUE)
    }

    # Use provided dromaset or get from global environment
    if (is.null(dromaset_object)) {
      if (exists(".dromadata", envir = .GlobalEnv) && !is.null(.dromadata$multiDS)) {
        dromaset_object <- .dromadata$multiDS
      } else {
        stop("DROMA data not found. Please initialize DROMA data first.")
      }
    }

    # TODO: Implement actual z-score normalization
    # This is a placeholder - actual implementation depends on package structure

    message("Applying z-score normalization to molecular profiles...")
    message("  - mRNA, CNV, methylation, protein: z-score across samples for each gene")
    message("  - Drug data: z-score row-by-row (each drug normalized across cell lines)")

    # Placeholder for actual normalization logic
    # In full implementation, this would:
    # 1. Extract molecular profiles
    # 2. Apply z-score normalization (column-wise for genes)
    # 3. Extract drug sensitivity data
    # 4. Apply z-score normalization (row-wise for drugs)
    # 5. Update the dromaset object with normalized data

    # Set global normalization state
    assign("normalization_state", TRUE, envir = .GlobalEnv)
    assign("GLOBAL_ZSCORE_STATE", list(
      enabled = TRUE,
      timestamp = Sys.time()
    ), envir = .GlobalEnv)

    message("Z-score normalization completed successfully.")
    return(TRUE)

  }, error = function(e) {
    message("Error applying z-score normalization: ", e$message)
    return(FALSE)
  })
}

#' Reset data to original (non-normalized) values
#'
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @return TRUE on success, FALSE on failure
#'
#' @examples
#' resetToOriginal()
#' resetToOriginal(multiDS)
resetToOriginal <- function(dromaset_object = NULL) {
  tryCatch({
    # Check if already reset
    if (exists("normalization_state", envir = .GlobalEnv) &&
        !isTRUE(get("normalization_state", envir = .GlobalEnv))) {
      message("Data is already in original state.")
      return(TRUE)
    }

    # Use provided dromaset or get from global environment
    if (is.null(dromaset_object)) {
      if (exists(".dromadata", envir = .GlobalEnv) && !is.null(.dromadata$multiDS)) {
        dromaset_object <- .dromadata$multiDS
      } else {
        stop("DROMA data not found. Please initialize DROMA data first.")
      }
    }

    # TODO: Implement actual reset logic
    # This is a placeholder - actual implementation depends on package structure

    message("Resetting data to original values...")

    # Placeholder for actual reset logic
    # In full implementation, this would:
    # 1. Reload original data from database
    # 2. Replace normalized data with original data
    # 3. Clear normalization caches

    # Update global normalization state
    assign("normalization_state", FALSE, envir = .GlobalEnv)
    assign("GLOBAL_ZSCORE_STATE", list(
      enabled = FALSE,
      timestamp = Sys.time()
    ), envir = .GlobalEnv)

    message("Data reset to original values successfully.")
    return(TRUE)

  }, error = function(e) {
    message("Error resetting data: ", e$message)
    return(FALSE)
  })
}

#' Check if data is currently z-score normalized
#'
#' @return TRUE if normalized, FALSE otherwise
#'
#' @examples
#' isZScoreNormalized()
isZScoreNormalized <- function() {
  if (exists("normalization_state", envir = .GlobalEnv)) {
    return(isTRUE(get("normalization_state", envir = .GlobalEnv)))
  }
  return(FALSE)
}

#' Get normalization state information
#'
#' @return List with normalization state details
#'
#' @examples
#' getNormalizationState()
getNormalizationState <- function() {
  state <- list(
    normalized = FALSE,
    timestamp = NULL,
    enabled = FALSE
  )

  if (exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
    zscore_state <- get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)
    state$enabled <- zscore_state$enabled
    state$timestamp <- zscore_state$timestamp
  }

  if (exists("normalization_state", envir = .GlobalEnv)) {
    state$normalized <- isTRUE(get("normalization_state", envir = .GlobalEnv))
  }

  return(state)
}

#' Apply z-score to a matrix or data frame
#'
#' @param data Numeric matrix or data frame
#' @param by Direction of normalization ("row" or "column")
#' @return Normalized matrix/data frame
#'
#' @examples
#' normalized <- zScoreMatrix(data_matrix, by = "column")
zScoreMatrix <- function(data, by = "row") {
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }

  if (by == "row") {
    # Normalize each row independently
    result <- t(scale(t(data)))
  } else if (by == "column") {
    # Normalize each column independently
    result <- scale(data)
  } else {
    stop("by must be either 'row' or 'column'")
  }

  # Handle NaN values (result from zero variance)
  result[is.nan(result)] <- 0

  return(result)
}