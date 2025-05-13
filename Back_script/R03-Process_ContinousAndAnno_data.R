load("Input/01/mRNA.Rda")
load("Input/01/cnv.Rda")
load("Input/01/meth.Rda")
load("Input/01/protein.Rda")
load("Input/02/drug.Rda")

# Filter ----
# Filter rows with too many zeros or NAs
filter_missing_data <- function(mat, min_valid = 3) {
  # For omics data (typically contains zeros)
  mat <- as.matrix(mat)
  if(all(is.numeric(mat))) {
    # Count valid (non-zero, non-NA) values in each row
    valid_counts <- apply(mat, 1, function(x) sum(x != 0 & !is.na(x)))
    # Keep rows with at least min_valid valid values
    keep_rows <- which(valid_counts >= min_valid)
    return(mat[keep_rows, ])
  } 
  # For drug data (typically contains NAs)
  else {
    # Count non-NA values in each row
    valid_counts <- apply(mat, 1, function(x) sum(!is.na(x)))
    # Keep rows with at least min_valid valid values
    keep_rows <- which(valid_counts >= min_valid)
    return(mat[keep_rows, ])
  }
}

# Apply missing data filtering to omics datasets
print("Filtering rows to keep only those with at least 3 valid (non-zero, non-NA) values:")

# Process omics data
if(exists("ccle_proteinms")) {
  rows_before <- nrow(ccle_proteinms)
  ccle_proteinms <- filter_missing_data(ccle_proteinms)
  print(paste("CCLE protein MS: removed", rows_before - nrow(ccle_proteinms), "rows, retained", nrow(ccle_proteinms)))
}

if(exists("ccle_proteinrppa")) {
  rows_before <- nrow(ccle_proteinrppa)
  ccle_proteinrppa <- filter_missing_data(ccle_proteinrppa)
  print(paste("CCLE protein RPPA: removed", rows_before - nrow(ccle_proteinrppa), "rows, retained", nrow(ccle_proteinrppa)))
}

if(exists("ccle_meth")) {
  rows_before <- nrow(ccle_meth)
  ccle_meth <- filter_missing_data(ccle_meth)
  print(paste("CCLE methylation: removed", rows_before - nrow(ccle_meth), "rows, retained", nrow(ccle_meth)))
}

if(exists("ccle_cnv")) {
  rows_before <- nrow(ccle_cnv)
  ccle_cnv <- filter_missing_data(ccle_cnv)
  print(paste("CCLE CNV: removed", rows_before - nrow(ccle_cnv), "rows, retained", nrow(ccle_cnv)))
}

if(exists("gdsc_cnv")) {
  rows_before <- nrow(gdsc_cnv)
  gdsc_cnv <- filter_missing_data(gdsc_cnv)
  print(paste("GDSC CNV: removed", rows_before - nrow(gdsc_cnv), "rows, retained", nrow(gdsc_cnv)))
}

if(exists("gCSI_cnv")) {
  rows_before <- nrow(gCSI_cnv)
  gCSI_cnv <- filter_missing_data(gCSI_cnv)
  print(paste("gCSI CNV: removed", rows_before - nrow(gCSI_cnv), "rows, retained", nrow(gCSI_cnv)))
}

if(exists("Xeva_cnv")) {
  rows_before <- nrow(Xeva_cnv)
  Xeva_cnv <- filter_missing_data(Xeva_cnv)
  print(paste("Xeva CNV: removed", rows_before - nrow(Xeva_cnv), "rows, retained", nrow(Xeva_cnv)))
}

if(exists("ccle_mRNA")) {
  rows_before <- nrow(ccle_mRNA)
  ccle_mRNA <- filter_missing_data(ccle_mRNA)
  print(paste("CCLE mRNA: removed", rows_before - nrow(ccle_mRNA), "rows, retained", nrow(ccle_mRNA)))
}

if(exists("gdsc_mRNA")) {
  rows_before <- nrow(gdsc_mRNA)
  gdsc_mRNA <- filter_missing_data(gdsc_mRNA)
  print(paste("GDSC mRNA: removed", rows_before - nrow(gdsc_mRNA), "rows, retained", nrow(gdsc_mRNA)))
}

if(exists("NCI60_mRNA")) {
  rows_before <- nrow(NCI60_mRNA)
  NCI60_mRNA <- filter_missing_data(NCI60_mRNA)
  print(paste("NCI60 mRNA: removed", rows_before - nrow(NCI60_mRNA), "rows, retained", nrow(NCI60_mRNA)))
}

if(exists("tavor_mRNA")) {
  rows_before <- nrow(tavor_mRNA)
  tavor_mRNA <- filter_missing_data(tavor_mRNA)
  print(paste("Tavor mRNA: removed", rows_before - nrow(tavor_mRNA), "rows, retained", nrow(tavor_mRNA)))
}

if(exists("UMPDO1_mRNA")) {
  rows_before <- nrow(UMPDO1_mRNA)
  UMPDO1_mRNA <- filter_missing_data(UMPDO1_mRNA)
  print(paste("UMPDO1 mRNA: removed", rows_before - nrow(UMPDO1_mRNA), "rows, retained", nrow(UMPDO1_mRNA)))
}

if(exists("UMPDO2_mRNA")) {
  rows_before <- nrow(UMPDO2_mRNA)
  UMPDO2_mRNA <- filter_missing_data(UMPDO2_mRNA)
  print(paste("UMPDO2 mRNA: removed", rows_before - nrow(UMPDO2_mRNA), "rows, retained", nrow(UMPDO2_mRNA)))
}

if(exists("UMPDO3_mRNA")) {
  rows_before <- nrow(UMPDO3_mRNA)
  UMPDO3_mRNA <- filter_missing_data(UMPDO3_mRNA)
  print(paste("UMPDO3 mRNA: removed", rows_before - nrow(UMPDO3_mRNA), "rows, retained", nrow(UMPDO3_mRNA)))
}

if(exists("Xeva_mRNA")) {
  rows_before <- nrow(Xeva_mRNA)
  Xeva_mRNA <- filter_missing_data(Xeva_mRNA)
  print(paste("Xeva mRNA: removed", rows_before - nrow(Xeva_mRNA), "rows, retained", nrow(Xeva_mRNA)))
}

# Apply missing data filtering to drug datasets
if(exists("prism_drug")) {
  rows_before <- nrow(prism_drug)
  prism_drug <- filter_missing_data(prism_drug)
  print(paste("PRISM drug: removed", rows_before - nrow(prism_drug), "rows, retained", nrow(prism_drug)))
}

if(exists("ctrp1_drug")) {
  rows_before <- nrow(ctrp1_drug)
  ctrp1_drug <- filter_missing_data(ctrp1_drug)
  print(paste("CTRP1 drug: removed", rows_before - nrow(ctrp1_drug), "rows, retained", nrow(ctrp1_drug)))
}

if(exists("ctrp2_drug")) {
  rows_before <- nrow(ctrp2_drug)
  ctrp2_drug <- filter_missing_data(ctrp2_drug)
  print(paste("CTRP2 drug: removed", rows_before - nrow(ctrp2_drug), "rows, retained", nrow(ctrp2_drug)))
}

if(exists("ccle_drug")) {
  rows_before <- nrow(ccle_drug)
  ccle_drug <- filter_missing_data(ccle_drug)
  print(paste("CCLE drug: removed", rows_before - nrow(ccle_drug), "rows, retained", nrow(ccle_drug)))
}

if(exists("FIMM_drug")) {
  rows_before <- nrow(FIMM_drug)
  FIMM_drug <- filter_missing_data(FIMM_drug)
  print(paste("FIMM drug: removed", rows_before - nrow(FIMM_drug), "rows, retained", nrow(FIMM_drug)))
}

if(exists("gCSI_drug")) {
  rows_before <- nrow(gCSI_drug)
  gCSI_drug <- filter_missing_data(gCSI_drug)
  print(paste("gCSI drug: removed", rows_before - nrow(gCSI_drug), "rows, retained", nrow(gCSI_drug)))
}

if(exists("gdsc1_drug")) {
  rows_before <- nrow(gdsc1_drug)
  gdsc1_drug <- filter_missing_data(gdsc1_drug)
  print(paste("GDSC1 drug: removed", rows_before - nrow(gdsc1_drug), "rows, retained", nrow(gdsc1_drug)))
}

if(exists("gdsc2_drug")) {
  rows_before <- nrow(gdsc2_drug)
  gdsc2_drug <- filter_missing_data(gdsc2_drug)
  print(paste("GDSC2 drug: removed", rows_before - nrow(gdsc2_drug), "rows, retained", nrow(gdsc2_drug)))
}

if(exists("GRAY_drug")) {
  rows_before <- nrow(GRAY_drug)
  GRAY_drug <- filter_missing_data(GRAY_drug)
  print(paste("GRAY drug: removed", rows_before - nrow(GRAY_drug), "rows, retained", nrow(GRAY_drug)))
}

if(exists("NCI60_drug")) {
  rows_before <- nrow(NCI60_drug)
  NCI60_drug <- filter_missing_data(NCI60_drug)
  print(paste("NCI60 drug: removed", rows_before - nrow(NCI60_drug), "rows, retained", nrow(NCI60_drug)))
}

if(exists("PDTXBreast_drug")) {
  rows_before <- nrow(PDTXBreast_drug)
  PDTXBreast_drug <- filter_missing_data(PDTXBreast_drug)
  print(paste("PDTXBreast drug: removed", rows_before - nrow(PDTXBreast_drug), "rows, retained", nrow(PDTXBreast_drug)))
}

if(exists("tavor_drug")) {
  rows_before <- nrow(tavor_drug)
  tavor_drug <- filter_missing_data(tavor_drug)
  print(paste("Tavor drug: removed", rows_before - nrow(tavor_drug), "rows, retained", nrow(tavor_drug)))
}

if(exists("UHNBreast_drug")) {
  rows_before <- nrow(UHNBreast_drug)
  UHNBreast_drug <- filter_missing_data(UHNBreast_drug)
  print(paste("UHNBreast drug: removed", rows_before - nrow(UHNBreast_drug), "rows, retained", nrow(UHNBreast_drug)))
}

if(exists("UMPDO1_drug")) {
  rows_before <- nrow(UMPDO1_drug)
  UMPDO1_drug <- filter_missing_data(UMPDO1_drug)
  print(paste("UMPDO1 drug: removed", rows_before - nrow(UMPDO1_drug), "rows, retained", nrow(UMPDO1_drug)))
}

if(exists("UMPDO2_drug")) {
  rows_before <- nrow(UMPDO2_drug)
  UMPDO2_drug <- filter_missing_data(UMPDO2_drug)
  print(paste("UMPDO2 drug: removed", rows_before - nrow(UMPDO2_drug), "rows, retained", nrow(UMPDO2_drug)))
}

if(exists("UMPDO3_drug")) {
  rows_before <- nrow(UMPDO3_drug)
  UMPDO3_drug <- filter_missing_data(UMPDO3_drug)
  print(paste("UMPDO3 drug: removed", rows_before - nrow(UMPDO3_drug), "rows, retained", nrow(UMPDO3_drug)))
}

if(exists("Xeva_drug")) {
  rows_before <- nrow(Xeva_drug)
  Xeva_drug <- filter_missing_data(Xeva_drug)
  print(paste("Xeva drug: removed", rows_before - nrow(Xeva_drug), "rows, retained", nrow(Xeva_drug)))
}

# Add index for samples and drugs ----
unique_drug_names <- unique(drug_anno$DrugName)
drug_id_map <- setNames(paste0("UM_DRUG_", seq_along(unique_drug_names)), unique_drug_names)
drug_anno$IndexID <- drug_id_map[drug_anno$DrugName]

unique_sample_names <- unique(sample_anno$SampleID)
sample_id_map <- setNames(paste0("UM_SAMPLE_", seq_along(unique_sample_names)), unique_sample_names)
sample_anno$IndexID <- sample_id_map[sample_anno$SampleID]

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
  Xeva_cnv,
  file = "Input/01/cnv.Rda"
)

save(
  ccle_mRNA,
  gdsc_mRNA,
  NCI60_mRNA,
  tavor_mRNA,
  UMPDO1_mRNA,
  UMPDO2_mRNA,
  UMPDO3_mRNA,
  Xeva_mRNA,
  file = "Input/01/mRNA.Rda"
)

save(
  # PDC
  tavor_drug,
  PDTXBreast_drug,
  # CellLine
  ctrp1_drug,
  ctrp2_drug,
  gdsc1_drug,
  gdsc2_drug,
  gCSI_drug,
  prism_drug,
  FIMM_drug,
  UHNBreast_drug,
  GRAY_drug,
  NCI60_drug,
  ccle_drug,
  # PDO
  UMPDO1_drug,
  UMPDO2_drug,
  UMPDO3_drug,
  # PDX
  Xeva_drug,
  file = "Input/02/drug.Rda"
)

save(
  drug_anno,
  sample_anno,
  file = "Input/03/anno.Rda"
)
