#' ==========================================
#' Dynamic Feature Loading Functions
#' ==========================================
#' These functions provide dynamic loading of features from DROMA database
#' with fallback options for error handling

#' Helper function to provide fallback features when database query fails
#'
#' @param feature_type The type of molecular feature (mRNA, meth, cnv, etc.)
#' @return A character vector of feature names
#'
#' @examples
#' getFallbackFeatures("mRNA")
#' getFallbackFeatures("mutation_gene")
getFallbackFeatures <- function(feature_type) {
  switch(feature_type,
         "mRNA" = c("ABCC3", "EGFR", "TP53", "KRAS", "MYC", "BRCA1", "BRCA2"),
         "meth" = c("MGMT", "CDKN2A", "MLH1", "RASSF1", "GSTP1"),
         "proteinrppa" = c("EGFR", "AKT", "ERK", "HER2", "BCL2"),
         "proteinms" = c("EGFR", "KRAS", "TP53", "MYC", "PIK3CA"),
         "cnv" = c("EGFR", "MYC", "TP53", "CDKN2A", "HER2", "BRCA1"),
         "drug" = c("Paclitaxel", "Cisplatin", "Doxorubicin", "5-FU", "Erlotinib",
                    "Gefitinib", "Trastuzumab", "Tamoxifen", "Imatinib"),
         "mutation_gene" = c("TP53", "KRAS", "EGFR", "BRAF", "PIK3CA", "PTEN", "NRAS",
                             "APC", "VHL", "IDH1", "IDH2"),
         "mutation_site" = c("TP53_R175H", "KRAS_G12D", "EGFR_L858R", "BRAF_V600E",
                            "PIK3CA_H1047R", "NRAS_Q61R", "IDH1_R132H"),
         "fusion" = c("EML4-ALK", "BCR-ABL", "TMPRSS2-ERG", "PML-RARA", "NPM-ALK"),
         character(0)
  )
}

#' Get available features dynamically from database with fallback
#'
#' @param feature_type The type of molecular feature
#' @param connection Optional database connection
#' @return A character vector of available feature names
#'
#' @examples
#' getAvailableFeatures("mRNA")
#' getAvailableFeatures("drug")
getAvailableFeatures <- function(feature_type, connection = NULL) {
  features <- tryCatch({
    # Initialize DROMA data if not already done
    if (!exists("initializeDROMAData", envir = .GlobalEnv)) {
      source("Modules/DataAdapter.R", local = FALSE)
    }

    if (!initializeDROMAData()) {
      return(getFallbackFeatures(feature_type))
    }

    # Special handling for drug features
    if (feature_type == "drug") {
      # Get available drugs from treatment metadata
      if (!is.null(.dromadata$multiDS@treatmentMetadata)) {
        drug_names <- unique(.dromadata$multiDS@treatmentMetadata$treatment_name)
        # Filter out empty values
        drug_names <- drug_names[drug_names != "" & !is.na(drug_names)]
        if (length(drug_names) > 0) {
          return(sort(drug_names))
        }
      }
      return(getFallbackFeatures("drug"))
    } else {
      # Get molecular features from database
      feature_data <- DROMA.Set::getFeatureFromDatabase(
        feature_type = feature_type,
        select_features = "all",
        connection = if(is.null(connection)) .dromadata$connection else connection
      )

      # Extract unique feature names
      if (!is.null(feature_data) && nrow(feature_data) > 0) {
        if (feature_type %in% c("mutation_gene", "mutation_site", "fusion")) {
          sort(unique(feature_data$feature_id))
        } else {
          sort(unique(feature_data$feature_id))
        }
      } else {
        getFallbackFeatures(feature_type)
      }
    }
  }, error = function(e) {
    message("Error loading features for ", feature_type, ": ", e$message)
    getFallbackFeatures(feature_type)
  })

  return(features)
}

#' Get available sample types with display names
#'
#' @return A named character vector of sample types
#'
#' @examples
#' getSampleTypes()
getSampleTypes <- function() {
  c("CellLine" = "Cell Lines",
    "PDC" = "Patient-Derived Cells",
    "PDO" = "Patient-Derived Organoids",
    "PDX" = "Patient-Derived Xenografts")
}

#' Filter features by frequency or availability
#'
#' @param features Character vector of feature names
#' @param min_samples Minimum number of samples required (optional)
#' @param connection Database connection
#' @return Filtered character vector
#'
#' @examples
#' filterFeaturesByFrequency(c("EGFR", "TP53", "MYC"), min_samples = 10)
filterFeaturesByFrequency <- function(features, min_samples = NULL, connection = NULL) {
  if (is.null(min_samples) || length(features) == 0) {
    return(features)
  }

  # TODO: Implement frequency filtering logic
  # This would query the database to count samples per feature
  # and filter out features with less than min_samples

  return(features)
}

#' Cache features for performance optimization
#'
#' @param feature_type The type of feature
#' @param features Character vector of features to cache
#' @return Invisible TRUE on success
#'
#' @examples
#' cacheFeatures("mRNA", c("EGFR", "TP53"))
cacheFeatures <- function(feature_type, features) {
  if (!exists(".feature_cache", envir = .GlobalEnv)) {
    assign(".feature_cache", list(), envir = .GlobalEnv)
  }

  cache <- get(".feature_cache", envir = .GlobalEnv)
  cache[[feature_type]] <- features
  assign(".feature_cache", cache, envir = .GlobalEnv)

  invisible(TRUE)
}

#' Get cached features if available
#'
#' @param feature_type The type of feature
#' @return Cached features or NULL if not available
#'
#' @examples
#' getCachedFeatures("mRNA")
getCachedFeatures <- function(feature_type) {
  if (exists(".feature_cache", envir = .GlobalEnv)) {
    cache <- get(".feature_cache", envir = .GlobalEnv)
    if (feature_type %in% names(cache)) {
      return(cache[[feature_type]])
    }
  }
  return(NULL)
}