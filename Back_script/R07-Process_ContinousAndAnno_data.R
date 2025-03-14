load("Input/01/protein.Rda")
load("Input/01/exp.Rda")
load("Input/01/cnv.Rda")
load("Input/01/meth.Rda")
load("Input/02/drug_0309.Rda")
load("Input/05/pdo_deng_0309.Rda")
load("Input/06/anno_0309.Rda")

# Filter ----
ccle_mRNA <- ccle_exp
gdsc_mRNA <- gdsc_exp

# Filter genes with zero standard deviation
filter_zero_sd <- function(mat) {
  # Calculate standard deviation for each gene (row), handling NA values
  sd_genes <- apply(mat, 1, function(x) sd(x, na.rm = TRUE))
  
  # Identify genes with non-zero standard deviation
  non_zero_sd <- which(sd_genes > 0 & !is.na(sd_genes))
  
  # Return the filtered matrix with only non-zero SD genes
  return(mat[non_zero_sd, ])
}

# Apply filtering to each dataset
ccle_proteinms <- filter_zero_sd(ccle_proteinms)
ccle_proteinrppa <- filter_zero_sd(ccle_proteinrppa)
ccle_meth <- filter_zero_sd(ccle_meth)
ccle_cnv <- filter_zero_sd(ccle_cnv)
gdsc_cnv <- filter_zero_sd(gdsc_cnv)
gCSI_cnv <- filter_zero_sd(gCSI_cnv)
ccle_mRNA <- filter_zero_sd(ccle_mRNA)
gdsc_mRNA <- filter_zero_sd(gdsc_mRNA)
deng1_mRNA <- filter_zero_sd(deng1_mRNA)
deng2_mRNA <- filter_zero_sd(deng2_mRNA)
deng3_mRNA <- filter_zero_sd(deng3_mRNA)

# Print number of genes before and after filtering
print("Genes before and after filtering (removing zero standard deviation genes):")
print(paste("CCLE protein MS:", nrow(ccle_proteinms), "genes retained"))
print(paste("CCLE protein RPPA:", nrow(ccle_proteinrppa), "genes retained"))
print(paste("CCLE methylation:", nrow(ccle_meth), "genes retained"))
print(paste("CCLE CNV:", nrow(ccle_cnv), "genes retained"))
print(paste("GDSC CNV:", nrow(gdsc_cnv), "genes retained"))
print(paste("gCSI CNV:", nrow(gCSI_cnv), "genes retained"))
print(paste("CCLE mRNA:", nrow(ccle_mRNA), "genes retained"))
print(paste("GDSC mRNA:", nrow(gdsc_mRNA), "genes retained"))
print(paste("Deng1 mRNA:", nrow(deng1_mRNA), "genes retained"))
print(paste("Deng2 mRNA:", nrow(deng2_mRNA), "genes retained"))
print(paste("Deng3 mRNA:", nrow(deng3_mRNA), "genes retained"))

# Process anno ----
drug_anno2 <- apply(drug_anno, 2, function(x){
  x[x %in% c(0, "")] <- NA
  x
})
drug_anno <- as.data.frame(drug_anno2)

# Save ----
save(
  ccle_proteinms,
  ccle_proteinrppa,
  file = "Input/01/protein.Rda"
)

save(
  ccle_meth,
  file = "Input/01/meth.Rda"
)

save(
  ccle_cnv,
  gdsc_cnv,
  gCSI_cnv,
  file = "Input/01/cnv.Rda"
)

save(
  ccle_mRNA,
  gdsc_mRNA,
  deng1_mRNA,
  deng2_mRNA,
  deng3_mRNA,
  file = "Input/01/exp.Rda"
)

save(
  prism_drug,
  ctrp1_drug,
  ctrp2_drug,
  gCSI_drug,
  gdsc1_drug,
  gdsc2_drug,
  deng1_drug,
  deng2_drug,
  deng3_drug,
  file = "Input/02/drug_0314.Rda"
)

save(
  drug_anno,
  sample_anno,
  file = "Input/06/anno_0314.Rda"
)
