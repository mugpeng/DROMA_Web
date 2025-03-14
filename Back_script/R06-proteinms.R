ccle_proteinms_raw <- fread("Input/protein_quant_current_normalized.csv.gz")
ccle_proteinms_raw <- as.data.frame(ccle_proteinms_raw)

ccle_proteinms <- ccle_proteinms_raw[,!grepl("Peptides", colnames(ccle_proteinms_raw))]
ccle_proteinms <- ccle_proteinms[,c(2,7:384)]

colnames(ccle_proteinms) <- gsub("_.*$", "", colnames(ccle_proteinms))
colnames(ccle_proteinms) <- toupper(colnames(ccle_proteinms))
ccle_proteinms_gene <- ccle_proteinms$GENE; ccle_proteinms$GENE <- NULL
ccle_proteinms <- as.matrix(ccle_proteinms)
ccle_proteinms[is.na(ccle_proteinms)] <- 0
ccle_proteinms1 <- ccle_proteinms[,colnames(ccle_proteinms) %in% colnames(ccle_proteinms)[duplicated(colnames(ccle_proteinms))]]
ccle_proteinms2 <- ccle_proteinms[,!colnames(ccle_proteinms) %in% colnames(ccle_proteinms)[duplicated(colnames(ccle_proteinms))]]
# ccle_proteinms has duplicated colnames, merge them by mean
ccle_proteinms1 <- sapply(unique(colnames(ccle_proteinms1)), function(col) {
  ms_df <- ccle_proteinms1[,colnames(ccle_proteinms1) %in% col]
  ms_df <- rowMeans(ms_df)
})
ccle_proteinms2 <- cbind(ccle_proteinms1, ccle_proteinms2)
ccle_proteinms2 <- as.data.frame(ccle_proteinms2)
ccle_proteinms2$GENE <- ccle_proteinms_gene
ccle_proteinms2 <- aggregate(.~GENE, max, data = ccle_proteinms2)
rownames(ccle_proteinms2) <- ccle_proteinms2$GENE
ccle_proteinms <- ccle_proteinms2
ccle_proteinms <- ccle_proteinms[!ccle_proteinms$GENE %in% "",]
ccle_proteinms$GENE <- NULL

# remain na in ccle_proteinms
ccle_proteinms <- as.matrix(ccle_proteinms)
ccle_proteinms[ccle_proteinms %in% 0] <- NA
ccle_proteinms <- as.data.frame(ccle_proteinms)

# Rename
load(file = "Input/01/protein.Rda")
ccle_proteinrppa <- ccle_protein
save(
  ccle_proteinms,
  ccle_proteinrppa,
  file = "Input/01/protein.Rda"
)
