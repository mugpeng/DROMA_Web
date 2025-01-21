
# meta calculation ----
metaCalConCon <- function(selected_pair){
  if(length(selected_pair) < 2) return(NULL)
  # test pairs one by one
  cal_list <- lapply(1:length(selected_pair), function(y){
    fea1_sel <- selected_pair[[y]][[1]]
    fea2_sel <- selected_pair[[y]][[2]]
    # Check for minimum length
    if(length(fea1_sel) < 3 || length(fea2_sel) < 3) return(NULL)
    cor_re <- tryCatch(
      cor.test(fea1_sel, fea2_sel,
               method = "pearson"),
      error = function(x){return(NULL)}
    )
    if(is.null(cor_re)) return(NULL)
    data.frame(
      p = cor_re$p.value,
      effect = cor_re$estimate,
      N = length(fea2_sel)
    )
  })
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 2) return(NULL)
  cal_re <- do.call(rbind, cal_list)
  cal_re$se <- sqrt((1 - cal_re$effect^2) / (cal_re$N - 2))
  cal_re$z <- 0.5 * log((1 + cal_re$effect) / (1 - cal_re$effect))  # Fisher's z 
  cal_re$se_z <- 1 / sqrt(cal_re$N - 3)         
  # cal_meta_re2 <- rma(yi = z, sei = se_z, data = cal_re, method = "REML")
  cal_meta_re <- tryCatch(
    suppressWarnings({metagen(TE = z, seTE = se_z, data = cal_re, sm = "Z",
                              control = list(maxiter = 2000,
                                             stepadj = 0.1,
                                             threshold = 0.000001)
    )}),
    error = function(x){return(NULL)}
  )
  cal_meta_re
}

metaCalConDis <- function(selected_pair){
  if(length(selected_pair) < 2) return(NULL)
  cal_list <- lapply(1:length(selected_pair), function(y){
    yes_drugs <- selected_pair[[y]][[1]]
    no_drugs <- selected_pair[[y]][[2]]
    
    # Check for minimum length
    if(length(yes_drugs) < 3 || length(no_drugs) < 3) return(NULL)
    
    wilcox_re <- tryCatch(
      wilcox.test(no_drugs, yes_drugs),
      error = function(x){return(NULL)}
    )
    if(is.null(wilcox_re)) return(NULL)
    
    cliff_delta <- tryCatch(
      cliff.delta(no_drugs, yes_drugs),
      error = function(x){return(NULL)}
    )
    if(is.null(cliff_delta)) return(NULL)
    
    data.frame(
      p = wilcox_re$p.value,
      effect = cliff_delta$estimate,
      N = length(yes_drugs) + length(no_drugs),
      n1 = length(yes_drugs),
      n2 = length(no_drugs)
    )
  })
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 2) return(NULL)
  cal_re <- do.call(rbind, cal_list)
  # Calculate standard error for Cliff's Delta
  cal_re$se <- sqrt((1 - cal_re$effect^2) * (cal_re$n1 + cal_re$n2 + 1) / 
                      (12 * cal_re$n1 * cal_re$n2))
  cal_meta_re <- tryCatch(
    suppressWarnings({meta_result <- metagen(TE = effect, 
                                             seTE = se, 
                                             data = cal_re,
                                             control = list(maxiter = 2000,
                                                            stepadj = 0.1,
                                                            threshold = 0.000001),
                                             sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
    )
    }),
    error = function(x){return(NULL)}
  )
  cal_meta_re
}

metaCalDisDis <- function(selected_pair) {
  # Check if we have enough pairs for meta-analysis
  if(length(selected_pair) < 2) return(NULL)
  
  # Calculate statistics for each pair
  cal_list <- lapply(1:length(selected_pair), function(y) {
    cont_table <- selected_pair[[y]]$cont_table
    
    # Skip if any cell has too few observations (e.g., < 3)
    if(any(cont_table < 3)) return(NULL)
    
    # Calculate odds ratio and its standard error
    tryCatch({
      # Extract values from contingency table
      a <- cont_table[1,1] # yes-yes
      b <- cont_table[1,2] # yes-no
      c <- cont_table[2,1] # no-yes
      d <- cont_table[2,2] # no-no
      
      # Calculate log odds ratio and its standard error
      log_or <- log((a * d)/(b * c))
      se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
      
      # Calculate Fisher's exact test p-value
      fisher_test <- fisher.test(cont_table)
      
      data.frame(
        log_or = log_or,
        se = se_log_or,
        p = fisher_test$p.value,
        N = sum(cont_table)
      )
    }, error = function(x) NULL)
  })
  
  # Remove NULL results and combine
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 2) return(NULL)
  
  cal_re <- do.call(rbind, cal_list)
  
  # Perform meta-analysis using random effects model
  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = log_or,
              seTE = se,
              data = cal_re,
              sm = "OR", # Specify odds ratio as summary measure
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )
    }),
    error = function(x) NULL
  )
  
  cal_meta_re
}

# Plot ----
# meta_df$p_adjust2 <- p.adjust(meta_df$p_value, method = "BH") 
# meta_df$p_adjust <- p.adjust(meta_df$p_value, method = "bonferroni") 
es_t = .4
P_t = .001
label = T
top_label_each = 5
title = NULL

# Plot volcano
plotMetaVolcano <- function(meta_df, 
                            es_t = .4, P_t = .001,
                            label = T,
                            top_label_each = 5,
                            title = NULL){
  meta_df$group <- dplyr::case_when(
    meta_df$effect_size > es_t & meta_df$p_value < P_t ~ "Up",
    meta_df$effect_size < -es_t & meta_df$p_value < P_t ~ "Down",
    TRUE ~ "NS"
  )  
  p <- ggplot(data = meta_df, 
              aes(x = effect_size, 
                  y = -log10(p_value))) +
    geom_point(size=2.5, alpha = 0.6, 
               aes(color = group)) +
    theme_bw() + theme(
      # text = element_text(family="serif"),
      legend.position = "none",
      title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, colour = "black"), 
      axis.text = element_text(size = 15, color = "black"), 
      legend.title = element_text(size = 15, colour = "black"),
      legend.text = element_text(size = 15),
      text =element_text(colour = "black"),
      axis.title.x = element_text(colour = "black")
    ) + 
    ylab("-log10(Pvalue)") + scale_color_manual(values = c("Down" = "#44bce4", 
                                                           "NS" = "grey", 
                                                           "Up" = "#fc7474")) + 
    geom_vline(xintercept=c(-es_t, es_t),lty=4,col="black",lwd=0.5) + 
    geom_hline(yintercept = -log10(P_t),lty=4,col="black",lwd=0.5) + ggtitle(title)
  meta_df2 <- meta_df[!meta_df$group %in% "NS",]
  forlabel_names <- c(
    meta_df2[head(order(meta_df2$effect_size), top_label_each),]$name,
    meta_df2[tail(order(meta_df2$effect_size), top_label_each),]$name
  )
  forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]
  if(dim(forlabel_df)[1] > 0){
    p <- p + 
      geom_point(size = 3, shape = 1, data = forlabel_df) +
      ggrepel::geom_text_repel(
        data = forlabel_df,
        aes(label = name),
        size = 5,
        color="black",
        nudge_y      = -0.2,
        max.overlaps = 1000,
      )
  }
  p
}

# Others ----
# Helper function to format time
format_time <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%d seconds", round(seconds)))
  } else if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    remaining_seconds <- round(seconds %% 60)
    return(sprintf("%d minutes %d seconds", minutes, remaining_seconds))
  } else {
    hours <- floor(seconds / 3600)
    remaining_minutes <- floor((seconds %% 3600) / 60)
    return(sprintf("%d hours %d minutes", hours, remaining_minutes))
  }
}

# Main function----
#' Main function to analyze batch features
#' @export
BatchFindSigFeaturesPlus <- function(feature1_type, 
                                     feature1_name, 
                                     feature2_type,
                                     cores = 1
) {  # Add cores parameter
  if(!all(c(feature1_type, feature2_type) %in% c("mRNA","cnv",
                                                 "meth", "proteinrppa", "proteinms", # continuous
                                                 "drug",
                                                 "mutation_gene", "mutation_site", "fusion" # discrete
  )
  )){
    stop("The select feature type doesn't exist. Please choose from drug, mRNA, meth, cnv, proteinms, proteinrppa, \nmutation_gene, mutation_site, or fusion.")
  }
  
  # Determine feature types
  continuous_types <- c("drug", "cnv", "proteinrppa", 
                        "proteinms", "meth", "mRNA")
  is_continuous1 <- feature1_type %in% continuous_types
  is_continuous2 <- feature2_type %in% continuous_types
  
  # Get selected specific feature1 data
  selected_feas1 <- selFeatures(feature1_type, feature1_name) 
  
  # Get compared feature2 vector
  feas_search_sel <- feas_search[feas_search$type %in% feature2_type,]
  # feas_search_sel <- head(feas_search_sel, 2000)
  
  # Define the worker function
  worker_function <- function(x) {
    results <- tryCatch({
      feature2_name <- feas_search_sel[x,1]
      selected_feas2 <- selFeatures(feature2_type, feature2_name) 
      
      if (is.null(selected_feas2)) return(NULL)
      
      # do statistics test based on four circumstances 
      # con vs con ----
      if (is_continuous1 && is_continuous2) {
        selected_pair <- pairDrugOmic_batch1(selected_feas1, selected_feas2)
        cal_meta_re <- metaCalConCon(selected_pair)
        # dis vs con ----
      } else if ((is_continuous1 && !is_continuous2) || (!is_continuous1 && is_continuous2)) {
        if (is_continuous1 && !is_continuous2){
          selected_pair <- pairDrugOmic_batch2(selected_feas2, selected_feas1)
        } else {
          selected_pair <- pairDrugOmic_batch2(selected_feas1, selected_feas2)
        } 
        cal_meta_re <- metaCalConDis(selected_pair)  
        # dis vs dis ----
      } else {
        selected_pair <- pairDrugOmic_batch3(selected_feas1, selected_feas2, 
                                             feature1_type = feature1_type,
                                             feature2_type = feature2_type,
                                             cells_search = cells_search)
        cal_meta_re <- metaCalDisDis(selected_pair)  
      }
      
      if(is.null(cal_meta_re)) return(NULL)
      results <- data.frame(
        p_value = cal_meta_re[["pval.random"]],
        effect_size = cal_meta_re[["TE.random"]],
        N = length(cal_meta_re[["studlab"]])
      )
    }, error = function(e) {
      # Log error but continue processing
      message(sprintf("Error processing feature %d: %s", x, e$message))
      NULL
    })
    # Update progress even if there was an error
    # if(!is.null(progress_callback)) {
    #   progress_callback(x)
    # }
    return(results)
  }
  
  # Use parallel processing if cores > 1
  if (cores > 1) {
    # Initialize snowfall
    sfInit(parallel = TRUE, cpus = cores)
    
    # Export required data and functions
    sfExport("selected_feas1", "feas_search_sel", 
             "is_continuous1", "is_continuous2",
             "feature1_type", "feature2_type", "cells_search", "fea_list")
    
    # Export functions
    sfExport("selFeatures", "pairDrugOmic_batch1", "pairDrugOmic_batch2", 
             "pairDrugOmic_batch3", "metaCalConCon", "metaCalConDis", "metaCalDisDis")
    
    # Export specific data
    feas_sel_vec <- fea_list[[feature2_type]]
    feas_sel_vec2 <- paste0(feas_sel_vec, "_", feature2_type)
    sfExport(list = feas_sel_vec2)
    
    # Load required packages on worker nodes
    sfLibrary(meta)
    sfLibrary(metafor)
    sfLibrary(effsize)
    
    # # Process in chunks
    # n_chunks <- ceiling(nrow(feas_search_sel) / chunk_size)
    # chunks <- split(1:nrow(feas_search_sel), 
    #                 cut(1:nrow(feas_search_sel), n_chunks, labels = FALSE))
    # 
    # # Process each chunk
    # cal_re_list <- unlist(lapply(chunks, function(chunk_indices) {
    #   chunk_results <- sfLapply(chunk_indices, worker_function)
    #   chunk_results[!sapply(chunk_results, is.null)]
    # }), recursive = FALSE)
    # cal_re_list <- sfClusterApplyLB(1:nrow(feas_search_sel), worker_function)
    
    # Run parallel computation
    cal_re_list <- sfLapply(1:nrow(feas_search_sel), worker_function)
    
    # Clean up snowfall
    sfStop()
  } else {
    # Run sequential computation
    cal_re_list <- lapply(1:nrow(feas_search_sel), worker_function)
  }
  
  # Process results
  valid_results <- !sapply(cal_re_list, is.null)
  fea_names <- feas_search_sel$name[valid_results]
  cal_re_list <- cal_re_list[valid_results]
  
  if (length(cal_re_list) == 0) {
    return(NULL)
  }
  
  cal_re_df <- do.call(rbind, cal_re_list)
  cal_re_df$name <- fea_names
  
  return(cal_re_df)
}

# Get pair data ----
# con with con
pairDrugOmic_batch1 <- function(myOmics, myDrugs){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      omic_sel <- na.omit(omic_sel); drug_sel <- na.omit(drug_sel)
      if(length(na.omit(omic_sel)) == 0 | length(na.omit(drug_sel)) == 0){ return(NULL) }
      intersected_cells <- intersect(names(omic_sel), names(drug_sel))
      omic_sel <- omic_sel[match(intersected_cells, names(omic_sel))]
      drug_sel <- drug_sel[match(intersected_cells, names(drug_sel))]
      list("omic" = omic_sel,
           "drug" = drug_sel)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  pair_list2
}

# Discrete with discrete
pairDrugOmic_batch2 <- function(myOmics, myDrugs){
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
  return(pair_list2)
}

# Discrete with discrete
pairDrugOmic_batch3 <- function(my_feas1, my_feas2, 
                                feature1_type, feature2_type,
                                cells_search) {
  # Create pairs list across all features in dataset 1
  pair_list3 <- lapply(1:length(my_feas1), function(x) {
    fea1_sel <- my_feas1[[x]]
    # For each feature in dataset 1, compare with all features in dataset 2
    pair_list <- lapply(1:length(my_feas2), function(y) {
      fea2_sel <- my_feas2[[y]]
      # Get all unique cells/samples
      # all_cells <- unique(c(fea1_sel, fea2_sel))
      # Skip if either feature has no data
      if (length(fea1_sel) == 0 || length(fea2_sel) == 0) {
        return(NULL)
      }
      cells_search_sel <- cells_search[cells_search$type %in% c(feature1_type, feature2_type) &
                                         cells_search$datasets %in% c(names(my_feas1)[x], names(my_feas2)[y]),]
      all_cells <- unique(cells_search_sel$cells)
      # Create 2x2 contingency table
      yes_yes <- length(intersect(fea1_sel, fea2_sel))
      yes_no <- length(fea1_sel) - yes_yes
      no_yes <- length(fea2_sel) - yes_yes
      no_no <- length(all_cells) - (yes_yes + yes_no + no_yes)
      
      # Skip if any cell count is too low
      if (any(c(yes_yes, yes_no, no_yes, no_no) < 0)) {
        return(NULL)
      }
      
      # Create contingency table
      cont_table <- matrix(
        c(yes_yes, yes_no, 
          no_yes, no_no),
        nrow = 2,
        dimnames = list(
          Feature1 = c("Yes", "No"),
          Feature2 = c("Yes", "No")
        )
      )
      
      list(
        "cont_table" = cont_table
      )
    })
    
    names(pair_list) <- paste0(names(my_feas1)[x], "_",
                               names(my_feas2))
    pair_list
  })
  
  # Flatten list and remove NULL entries
  pair_list3 <- unlist(pair_list3, recursive = FALSE)
  pair_list3 <- pair_list3[!sapply(pair_list3, is.null)]
  
  return(pair_list3)
}
