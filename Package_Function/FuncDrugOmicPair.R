# Con ----
# Function to pair omics and drug data
pairDrugOmic <- function(myOmics, myDrugs, merged = FALSE){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      omic_sel <- na.omit(omic_sel); drug_sel <- na.omit(drug_sel)
      if(length(na.omit(omic_sel)) == 0 | length(na.omit(drug_sel)) == 0){ return(NULL) }
      intersected_cells <- intersect(names(omic_sel), names(drug_sel))
      omic_sel <- omic_sel[match(intersected_cells, names(omic_sel))]
      drug_sel <- drug_sel[match(intersected_cells, names(drug_sel))]
      if(length(na.omit(omic_sel)) == 0 | length(na.omit(drug_sel)) == 0){ return(NULL) }
      list("omic" = omic_sel,
           "drug" = drug_sel)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  if(length(pair_list2) < 1) {stop("Please try to another drug-omic pair. This pair do not have result.")}
  # If merged is TRUE and z-score normalization is enabled, create a merged dataset
  if(merged) {
    # Create a merged dataset by combining all pairs using the more efficient approach
    combined_list <- list(
      # Combine all omic vectors
      unlist(lapply(pair_list2, function(x) x$omic)),
      
      # Combine all drug vectors
      unlist(lapply(pair_list2, function(x) x$drug))
    )
    pair_list2[["merged_dataset"]] <- list(
      "omic" = combined_list[[1]],
      "drug" = combined_list[[2]]
    )
  }
  pair_list2
}

# Perform meta-analysis on continuous drug-omic pairs
analyze_continuous_drugomic <- function(myPairs) {
  # Initialize list to store correlation results
  test_list <- list()
  valid_indices <- c()
  
  # Analyze each pair
  for (x in seq_along(myPairs)) {
    # Skip merged dataset for meta-analysis
    if (names(myPairs)[x] == "merged_dataset") next
    
    # Try to perform correlation test
    tryCatch({
      omic_sel <- myPairs[[x]]$omic
      drug_sel <- myPairs[[x]]$drug
      
      # Check for valid data
      if (length(omic_sel) < 3 || length(drug_sel) < 3) next
      
      # Perform correlation test
      cor_re <- cor.test(omic_sel, drug_sel, method = "pearson")
      test_list[[x]] <- data.frame(
        p = cor_re$p.value,
        effect = cor_re$estimate,
        N = length(omic_sel)
      )
      # Track valid indices for proper study name mapping
      valid_indices <- c(valid_indices, x)
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Return NULL if no correlations could be calculated
  if (length(test_list) < 1) return(NULL)
  
  # Prepare data for meta-analysis
  meta_df <- do.call(rbind, test_list)
  # Use valid_indices to correctly map study names
  meta_df$study <- names(myPairs)[valid_indices]
  meta_df$se <- sqrt((1 - meta_df$effect^2) / (meta_df$N - 2))
  meta_df$z <- 0.5 * log((1 + meta_df$effect) / (1 - meta_df$effect))  # Fisher's z 
  meta_df$se_z <- 1 / sqrt(meta_df$N - 3)         
  
  # Perform meta-analysis
  tryCatch({
    # Only perform meta-analysis if we have at least 2 studies
    if (nrow(meta_df) >= 2) {
      cal_meta_re <- metagen(TE = z, 
                             seTE = se_z, 
                             data = meta_df, 
                             sm = "Z",
                             studlab = study)
      return(cal_meta_re)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

# Discrete ----
# Function to pair discrete omics and drug data
pairDrugOmic2 <- function(myOmics, myDrugs, merged = FALSE){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      yes_drugs <- na.omit(drug_sel[names(drug_sel) %in% omic_sel])
      no_drugs <- na.omit(drug_sel[!names(drug_sel) %in% omic_sel])
      if(length(yes_drugs) == 0 | length(no_drugs) == 0){ 
        return(NULL)
      }
      list(yes = yes_drugs,
           no = no_drugs)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  if(length(pair_list2) < 1){stop("Please try to another drug-omic pair. This pair do not have result.")}
  # If merged is TRUE and z-score normalization is enabled, create a merged dataset
  if(merged) {
    # Create a merged dataset by combining all pairs using the more efficient approach
    combined_list <- list(
      # Combine all yes vectors
      unlist(lapply(pair_list2, function(x) x$yes)),
      
      # Combine all no vectors
      unlist(lapply(pair_list2, function(x) x$no))
    )
    
    # Only add merged dataset if we have enough data points
    pair_list2[["merged_dataset"]] <- list(
      "yes" = combined_list[[1]],
      "no" = combined_list[[2]]
    )
  }
  pair_list2
}

# Perform meta-analysis on discrete drug-omic pairs
analyze_discrete_drugomic <- function(myPairs) {
  # Initialize list to store test results
  test_list <- list()
  valid_indices <- c()

  # Analyze each pair
  for (x in seq_along(myPairs)) {
    # Skip merged dataset for meta-analysis
    if (names(myPairs)[x] == "merged_dataset") next
    
    # Try to perform Wilcoxon test and effect size calculation
    tryCatch({
      yes_drugs <- myPairs[[x]]$yes
      no_drugs <- myPairs[[x]]$no
      
      # Check for valid data
      if (length(yes_drugs) < 3 || length(no_drugs) < 3) next
      
      # Perform statistical test
      wilcox_test <- wilcox.test(no_drugs, yes_drugs)
      cliff_delta <- cliff.delta(no_drugs, yes_drugs)
      
      # Store results
      test_list[[x]] <- data.frame(
        p = wilcox_test$p.value,
        effect = cliff_delta$estimate,
        N = length(yes_drugs) + length(no_drugs),
        n1 = length(yes_drugs),
        n2 = length(no_drugs)
      )
      valid_indices <- c(valid_indices, x)
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Return NULL if no tests could be performed
  if (length(test_list) < 1) return(NULL)
  
  # Prepare data for meta-analysis
  meta_df <- do.call(rbind, test_list)
  meta_df$study <- names(myPairs)[valid_indices]
  
  # Calculate standard error for Cliff's Delta
  meta_df$se <- sqrt((1 - meta_df$effect^2) * (meta_df$n1 + meta_df$n2 + 1) / 
                       (12 * meta_df$n1 * meta_df$n2))
  
  # Perform meta-analysis
  tryCatch({
    # Only perform meta-analysis if we have at least 2 studies
    if (nrow(meta_df) >= 2) {
      meta_result <- metagen(TE = effect, 
                             seTE = se, 
                             data = meta_df,
                             sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
                             studlab = study)
      return(meta_result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

# All ----
# Function to handle both continuous and discrete omics data in 
# one function
oneDrugOmicPair <- function(select_omics_type, select_omics,
                            select_drugs,
                            data_type = "all", tumor_type = "all",
                            merged_enabled = TRUE,
                            meta_enabled = TRUE){
  # Get drug data
  myDrugs <- selFeatures("drug", select_drugs, 
                        data_type = data_type, 
                        tumor_type = tumor_type)
  myOmics <- selFeatures(select_omics_type, select_omics,
                         data_type = data_type, 
                         tumor_type = tumor_type)
  # Initialize result list
  result <- list()
  
  # Handle continuous omics data
  if(select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")){
    # Pair data
    myPairs <- pairDrugOmic(myOmics, myDrugs, merged = merged_enabled)
    
    # Create plots
    result$plot <- plot_all_continuous_drugomic(myPairs)
    
    # Perform meta-analysis
    if(meta_enabled){
      meta_result <- analyze_continuous_drugomic(myPairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    # Store data
    result$data <- myPairs
  } else {
    # Pair data
    myPairs <- pairDrugOmic2(myOmics, myDrugs, merged = merged_enabled)
    
    # Create plots
    result$plot <- plot_all_discrete_drugomic(myPairs)
    
    # Perform meta-analysis
    if(meta_enabled){
      meta_result <- analyze_discrete_drugomic(myPairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    # Store data
    result$data <- myPairs
  }
  
  # Return results if there's a plot
  if (is.null(result$plot)) {
    return(list())
  } else {
    return(result)
  }
}

