# Process mutation data variables
process_mutation_data <- function(variable_names) {
  for(var_name in variable_names) {
    data <- base::get(var_name)
    
    # Get prefix by removing "mut"
    prefix <- gsub("mut", "", var_name)
    
    # Handle two-column datasets (gCSI_mut and Xeva_mut)
    if(ncol(data) == 2) {
      # Rename directly to xx_mutation_gene
      new_name <- paste0(prefix, "mutation_gene")
      assign(new_name, data, envir = .GlobalEnv)
    } 
    # Handle four-column datasets
    else if(ncol(data) == 4) {
      # Create gene dataset (columns 1,2)
      gene_data <- unique(data[, c(1, 2)])
      gene_name <- paste0(prefix, "mutation_gene")
      
      # Create site dataset (columns 4,2)
      site_data <- unique(data[, c(4, 2)])
      site_name <- paste0(prefix, "mutation_site")
      
      # Assign both new datasets to global environment
      assign(gene_name, gene_data, envir = .GlobalEnv)
      assign(site_name, site_data, envir = .GlobalEnv)
    }
    else {
      warning(paste("Unexpected number of columns in", var_name, "- skipping"))
    }
    
    # Remove original variable
    if(exists(var_name, envir = .GlobalEnv)) {
      rm(list = var_name, envir = .GlobalEnv)
    }
  }
}