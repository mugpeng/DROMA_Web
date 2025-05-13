bright_palette_26 <- c(
  "#FF5733", # Bright red-orange
  "#33C4FF", # Bright sky blue
  "#FF33E9", # Bright pink
  "#33FF57", # Bright green
  "#FFD133", # Bright yellow
  "#8B33FF", # Bright purple
  "#FF8B33", # Bright orange
  "#33FFC4", # Bright teal
  "#FF3355", # Bright red
  "#33FFED", # Bright cyan
  "#C433FF", # Bright violet
  "#BBFF33", # Bright lime
  "#FF33BB", # Bright magenta
  "#33FFB2", # Bright mint
  "#3357FF", # Bright blue
  "#FF9933", # Bright amber
  "#33FF78", # Bright seafoam
  "#ED33FF", # Bright fuchsia
  "#57FF33", # Bright chartreuse
  "#FF33C4", # Bright hot pink
  "#33A5FF", # Bright azure
  "#FFB233", # Bright gold
  "#3378FF", # Bright royal blue
  "#FF5733", # Bright vermilion
  "#33FF9E", # Bright turquoise
  "#D6FF33"  # Bright yellow-green
)

#' Process Drug Sensitivity Data
#'
#' Creates a combined dataframe with raw and normalized drug sensitivity values
#'
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @return A dataframe with combined raw and normalized drug sensitivity values
#' @export
process_drug_data <- function(drug_name, data_type = "all", tumor_type = "all") {
  if (drug_name == "") {
    stop("Please select a drug.")
  }
  
  # Get z-score normalized data
  normalized_data_list <- selFeatures(
    select_feas_type = "drug", 
    select_feas = drug_name,
    data_type = data_type, 
    tumor_type = tumor_type
  )
  
  # Get raw data
  raw_data_list <- selFeatures(
    select_feas_type = "drug_raw", 
    select_feas = drug_name,
    data_type = data_type, 
    tumor_type = tumor_type
  )
  
  # Get all study names (using union of keys from both lists)
  all_studies <- union(names(normalized_data_list), names(raw_data_list))
  
  # Create an empty list to hold combined data for each study
  combined_data <- list()
  
  # For each study, process and combine the data
  for (study in all_studies) {
    if (study %in% names(normalized_data_list) && study %in% names(raw_data_list)) {
      # Get values for this study
      normalized_values <- normalized_data_list[[study]]
      raw_values <- raw_data_list[[study]]
      
      # Get common sample IDs
      common_samples <- intersect(names(normalized_values), names(raw_values))
      
      if (length(common_samples) > 0) {
        # Create dataframe with both values
        study_df <- data.frame(
          sampleid = common_samples,
          zscore_value = as.numeric(normalized_values[common_samples]),
          raw_value = as.numeric(raw_values[common_samples]),
          study = study,
          stringsAsFactors = FALSE
        )
        
        combined_data[[study]] <- study_df
      }
    }
  }
  
  # Combine all studies into a single dataframe
  if (length(combined_data) > 0) {
    result_df <- bind_rows(combined_data)
    result_df <- na.omit(result_df)
    
    return(result_df)
  } else {
    return(NULL)
  }
}

#' Annotate Drug Sensitivity Data
#'
#' Adds sample annotations to drug sensitivity data
#'
#' @param drug_data Dataframe containing drug sensitivity data from process_drug_data
#' @param sample_annotations Dataframe containing sample annotations (default: uses global sample_anno)
#' @return A dataframe with drug sensitivity data and sample annotations
#' @export
annotate_drug_data <- function(drug_data, sample_annotations = NULL) {
  if (is.null(drug_data)) {
    return(NULL)
  }
  
  # Use provided sample annotations or global sample_anno
  annotations <- sample_annotations
  if (is.null(annotations)) {
    if (exists("sample_anno", envir = .GlobalEnv)) {
      annotations <- base::get("sample_anno", envir = .GlobalEnv)
    } else {
      warning("sample_anno object not found. Returning non-annotated data.")
      return(drug_data)
    }
  }
  
  # Merge drug data with annotations
  merged_df <- left_join(drug_data, 
                         annotations %>% select(-ProjectID), # Remove ProjectID from the join
                         by = c("sampleid" = "SampleID"))
  
  # Clean up
  if ("ProjectRawName" %in% colnames(merged_df)) {
    merged_df$ProjectRawName <- NULL
  }
  
  merged_df <- unique(merged_df)
  return(merged_df)
}

#' Format Drug Data Table
#'
#' Creates a formatted datatable for drug sensitivity data
#'
#' @param drug_data Dataframe containing drug sensitivity data
#' @param caption Caption text for the table
#' @return A formatted DT::datatable object
#' @export
format_drug_table <- function(drug_data, caption = "Drug sensitivity data - showing both raw and Z-score normalized values") {
  if (is.null(drug_data) || nrow(drug_data) == 0) {
    return(NULL)
  }
  
  DT::datatable(
    drug_data,
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
      htmltools::strong(caption)
    ),
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons',
    rownames = FALSE,
    filter = 'top'
  ) %>%
  # Format numeric columns to 3 decimal places
  DT::formatRound(columns = c('zscore_value', 'raw_value'), digits = 3)
}

#' Get Drug Sensitivity Data
#'
#' Wrapper function that processes and returns drug sensitivity data
#'
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param include_annotations Logical indicating whether to include sample annotations
#' @param sample_annotations Optional dataframe containing sample annotations
#' @return A dataframe with drug sensitivity data
#' @export
get_drug_sensitivity_data <- function(drug_name, 
                                     data_type = "all", 
                                     tumor_type = "all",
                                     include_annotations = TRUE,
                                     sample_annotations = NULL) {
  # Process drug data
  drug_data <- process_drug_data(drug_name, data_type, tumor_type)
  
  # Add annotations if requested
  if (include_annotations) {
    drug_data <- annotate_drug_data(drug_data, sample_annotations)
  }
  
  return(drug_data)
}
