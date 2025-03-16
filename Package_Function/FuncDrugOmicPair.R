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
  
  if(length(pair_list2) < 1) {stop("Please try to another drug-omic pair. This pair do not have result.")}
  pair_list2
}

# Function to plot paired omics and drug data (continuous)
plotDrugOmicPair_con <- function(myPairs, meta = T){
  p_list <- list()
  cor_list <- list()
  for(x in 1:length(myPairs)){
    # Try to perform correlation test and plotting
    tryCatch({
      omic_sel <- myPairs[[x]][[1]]
      drug_sel <- myPairs[[x]][[2]]
      
      # Check for valid data
      if(length(omic_sel) < 3 || length(drug_sel) < 3 ) next
      
      cor_df <- data.frame(
        genes = omic_sel,
        drugs = drug_sel
      )
      
      # Create plot
      current_plot <- ggscatter(cor_df, x = "genes", y = "drugs",
                                alpha = 0.2) +
        stat_cor(size = 6, method = "spearman") + 
        stat_smooth(formula = y ~ x, method = "lm") + 
        theme_bw() +
        theme(
          axis.title = element_blank(),
          title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12)
        ) + 
        ggtitle(names(myPairs)[x])
      
      p_list[[x]] <- current_plot
      
      # Perform correlation test if meta-analysis is requested
      # Skip meta-analysis for merged dataset
      if(meta == T && names(myPairs)[x] != "merged_dataset"){
        cor_re <- cor.test(omic_sel, drug_sel, method = "pearson")
        cor_list[[x]] <- data.frame(
          p = cor_re$p.value,
          effect = cor_re$estimate,
          N = length(omic_sel)
        )
      }
    }, error = function(e) {
    })
  }
  # Remove NULL entries from lists
  p_list <- p_list[!sapply(p_list, is.null)]
  # Prepare return list
  re <- list()
  
  # Only create plot if we have valid plots
  if(length(p_list) > 0) {
    re[[1]] <- wrap_plots(p_list, ncol = 3)
  } else {
    re[[1]] <- NULL
  }
  
  # Perform meta-analysis only if we have correlation results
  if(meta == T && length(cor_list) > 1) {
    tryCatch({
      cal_re <- do.call(rbind, cor_list)
      cal_re$study <- names(myPairs)[!names(myPairs) %in% "merged_dataset"][!sapply(cor_list, is.null)]
      cal_re$se <- sqrt((1 - cal_re$effect^2) / (cal_re$N - 2))
      cal_re$z <- 0.5 * log((1 + cal_re$effect) / (1 - cal_re$effect))  # Fisher's z 
      cal_re$se_z <- 1 / sqrt(cal_re$N - 3)         
      
      cal_meta_re <- metagen(TE = z, 
                             seTE = se_z, 
                             data = cal_re, 
                             sm = "Z",
                             studlab = study)
      re[[2]] <- cal_meta_re
    }, error = function(e) {
      # message("Meta-analysis failed: ", as.character(e))
      re[[2]] <- NULL
    })
  } else {
    re[[2]] <- NULL
  }
  if(length(re) < 1) {stop("Please try to another drug-omic pair. This pair do not have result.")}
  return(re)
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
  
  if(length(pair_list2) < 1){stop("Please try to another drug-omic pair. This pair do not have result.")}
  pair_list2
}

# Function to plot paired omics and drug data (discrete)
plotDrugOmicPair_dis <- function(myPairs, meta = T){
  p_list <- list()
  test_list <- list()
  
  for(x in 1:length(myPairs)){
    # Try to perform wilcoxon test and plotting
    tryCatch({
      yes_drugs <- myPairs[[x]][[1]]
      no_drugs <- myPairs[[x]][[2]]
      
      # Check for valid data
      if(length(yes_drugs) < 3 || length(no_drugs) < 3) next
      
      box_df <- data.frame(
        drugs = c(no_drugs, yes_drugs),
        events = rep(c("no","yes"), times = c(length(no_drugs), length(yes_drugs))))
      
      # Create plot
      current_plot <- ggboxplot(data = box_df, x = "events", y = "drugs",
                                fill = "events", palette = c("#BEBADAFF", "#FB8072FF"),
                                add = "jitter", add.params = list(alpha = 0.15)) + 
        stat_compare_means(size = 6, label.x = 0.8,
                           label.y = (max(box_df$drugs) - max(box_df$drugs)/8),
                           label = "p.format") + 
        theme_bw() + 
        theme(
          axis.title = element_blank(),
          title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12),
          legend.position = "none"
        ) + 
        coord_cartesian(ylim = c(0, max(box_df$drugs) + 
                                   max(box_df$drugs)/20)) + 
        ggtitle(names(myPairs)[x])
      
      p_list[[x]] <- current_plot
      
      # Calculate effect size and other statistics for meta-analysis
      # Skip meta-analysis for merged dataset
      if(meta == T && names(myPairs)[x] != "merged_dataset"){
        wilcox_test <- wilcox.test(no_drugs, yes_drugs)
        cliff_delta <- cliff.delta(no_drugs, yes_drugs)
        
        test_list[[x]] <- data.frame(
          p = wilcox_test$p.value,
          effect = cliff_delta$estimate,
          N = length(yes_drugs) + length(no_drugs),
          n1 = length(yes_drugs),
          n2 = length(no_drugs)
        )
      }
    }, error = function(e) {
    })
  }
  
  # Remove NULL entries from lists
  p_list <- p_list[!sapply(p_list, is.null)]
  
  # Prepare return list
  re <- list()
  
  # Only create plot if we have valid plots
  if(length(p_list) > 0) {
    re[[1]] <- wrap_plots(p_list, ncol = 3)
  } else {
    re[[1]] <- NULL
  }
  
  # Perform meta-analysis only if we have test results
  if(meta == T && length(test_list) > 1) {
    tryCatch({
      # Combine all test results
      meta_df <- do.call(rbind, test_list)
      meta_df$study <- names(myPairs)[!names(myPairs) %in% "merged_dataset"][!sapply(test_list, is.null)]
      
      # Calculate standard error for Cliff's Delta
      meta_df$se <- sqrt((1 - meta_df$effect^2) * (meta_df$n1 + meta_df$n2 + 1) / 
                           (12 * meta_df$n1 * meta_df$n2))
      
      # Perform meta-analysis
      meta_result <- metagen(TE = effect, 
                             seTE = se, 
                             data = meta_df,
                             sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
                             studlab = study)
      
      re[[2]] <- meta_result
    }, error = function(e) {
      re[[2]] <- NULL
    })
  } else {
    re[[2]] <- NULL
  }
  if(length(re) < 1) {stop("Please try to another drug-omic pair. This pair do not have result.")}
  return(re)
}

# All ----
# Function to handle both continuous and discrete omics data in 
# one function
oneDrugOmicPair <- function(select_omics_type, select_omics,
                            select_drugs){
  myDrugs <- selDrugs(select_drugs)
  
  # Check if z-score normalization is enabled
  merged_enabled <- FALSE
  if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
    merged_enabled <- isTRUE(get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$enabled)
  }
  
  if(select_omics_type %in% c("exp", "meth", "proteinrppa", "cnv",
                              "proteinms")){
    myOmics <- selOmics(select_omics_type, select_omics)
    myPairs <- pairDrugOmic(myOmics, myDrugs, merged = merged_enabled)
    p <- plotDrugOmicPair_con(myPairs)
  } else{
    myOmics <- selOmics2(select_omics_type, select_omics)
    myPairs <- pairDrugOmic2(myOmics, myDrugs, merged = merged_enabled)
    p <- plotDrugOmicPair_dis(myPairs)
  }
  return(p)
}

