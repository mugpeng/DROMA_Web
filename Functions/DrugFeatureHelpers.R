#' Drug Feature Visualization Helper Functions
#' 
#' This file contains helper functions for drug feature analysis visualization.
#' These functions wrap DROMA_R package functions for use in the Shiny app.

#' Create Drug Comparison Plot
#'
#' @description Creates comparison plots for drug sensitivity data based on different comparison variables.
#' This function wraps DROMA_R package visualization functions.
#'
#' @param data Data frame containing drug sensitivity data with comparison variable
#' @param comparison_var Character string specifying the variable to compare by
#' @param value_column Character string specifying the value column name (default: "value")
#' @param value_label Character string for the value label (default: "Drug Sensitivity")
#' @param num_bins Numeric, number of bins for continuous variables (default: 4)
#' @param show_groups_boxplot Logical, whether to show grouped boxplot for continuous variables (default: TRUE)
#' @return A ggplot2 object or grid of plots
#'
#' @examples
#' \dontrun{
#' # For categorical variable
#' plot <- create_drug_comparison_plot(data, comparison_var = "TumorType")
#' 
#' # For continuous variable
#' plot <- create_drug_comparison_plot(data, comparison_var = "Age", 
#'                                     num_bins = 4, show_groups_boxplot = TRUE)
#' }
create_drug_comparison_plot <- function(data, comparison_var, value_column = "value", 
                                        value_label = "Drug Sensitivity", 
                                        num_bins = 4, show_groups_boxplot = TRUE) {
  # Handle missing values in the comparison variable
  data <- data[!is.na(data[[comparison_var]]), ]
  
  if (nrow(data) == 0) {
    return(ggplot2::ggplot() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                             label = "No data available for this comparison") + 
             ggplot2::theme_void())
  }
  
  # Check if the comparison variable is numeric/continuous
  if (is.numeric(data[[comparison_var]])) {
    # For continuous variables - use DROMA_R functions
    p1 <- DROMA.R::plotContinuousComparison(
      data = data, 
      cont_column = comparison_var, 
      value_column = value_column, 
      value_label = value_label
    )
    
    # Also create a boxplot with bins if requested
    if (show_groups_boxplot) {
      # Create grouped boxplot using DROMA_R function
      p2 <- DROMA.R::plotContinuousGroups(
        data = data, 
        cont_column = comparison_var, 
        value_column = value_column, 
        value_label = value_label,
        num_bins = num_bins
      )
      
      # Return grid of both plots
      return(gridExtra::grid.arrange(p1, p2, ncol = 2))
    }
    
    return(p1)
  } else {
    # For categorical variables - use DROMA_R function
    return(DROMA.R::plotCategoryComparison(
      data = data, 
      category_column = comparison_var, 
      value_column = value_column, 
      value_label = value_label
    ))
  }
}

