# Function to select feature data include omic data and drug data
selFeatures <- function(select_feas_type, select_feas) {
  if(!select_feas_type %in% c("mRNA","cnv",
                               "meth", "proteinrppa", "proteinms", # continuous
                               "mutation_gene", "mutation_site", "fusion", # discrete
                               "drug"
  )){
    stop("The select feature type doesn't exsit. Please choose from drug, mRNA, meth, cnv, proteinms, proteinrppa, \nmutation_gene, mutation_site, or fusion.")
  }
  feas_sel_vec <- fea_list[[select_feas_type]]
  if(select_feas_type %in% c("mRNA","cnv",
                              "meth", "proteinrppa", "proteinms",
                              "drug")){
    fea_sel_list <- lapply(feas_sel_vec, function(x){
      # x = feas_sel_vec[1]
      fea <- as.data.frame(base::get(paste0(x, "_", select_feas_type)))
      fea_sel <- as.numeric(fea[select_feas,])
      if(all(is.na(fea_sel))){
        return(NULL)
      }
      names(fea_sel) <- colnames(fea)
      fea_sel
    })
  } else{
    fea_sel_list <- lapply(feas_sel_vec, function(x){
      # x = feas_sel_vec[1]
      fea <- as.data.frame(base::get(paste0(x, "_", select_feas_type)))
      fea_sel <- fea$cells[fea[[1]] %in% select_feas]
      if(all(is.na(fea_sel))){
        return(NULL)
      }
      fea_sel
    })
  }
  names(fea_sel_list) <- feas_sel_vec
  fea_sel_list <- fea_sel_list[!sapply(fea_sel_list, is.null)]
  if(length(fea_sel_list) < 1){stop("The select feature doesn't exsit, check if it is in your selected feature types. \nsuch as, 'YM-155' for 'drug'.")}
  fea_sel_list
}

