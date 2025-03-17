tmp <- list()

# Drugs-omics pairs analysis----
## Search ----
### Omics ----
tmp$omics_search_CNV <- data.frame(
  omics = c(rownames(ccle_cnv),
            rownames(gdsc_cnv),
            rownames(gCSI_cnv)
  ),
  type = "cnv"
) %>% unique()


tmp$omics_search_mRNA <- data.frame(
  omics = c(rownames(ccle_mRNA),
            rownames(gdsc_mRNA),
            rownames(deng1_mRNA),
            rownames(deng2_mRNA),
            rownames(deng3_mRNA)
  ),
  type = "mRNA"
) %>% unique()

tmp$omics_search_meth <- data.frame(
  omics = c(rownames(ccle_meth)),
  type = "meth"
) %>% unique()

tmp$omics_search_proteinrppa <- data.frame(
  omics = c(rownames(ccle_proteinms)),
  type = "proteinms"
) %>% unique()

tmp$omics_search_proteinms <- data.frame(
  omics = c(rownames(ccle_proteinrppa)),
  type = "proteinrppa"
) %>% unique()

tmp$omics_search_mutgenes <- data.frame(
  omics = c(ccle_mut$genes,
            gdsc_mut$genes,
            gCSI_mut$genes
  ),
  type = "mutation_gene"
) %>% unique()

tmp$omics_search_mutsites <- data.frame(
  omics = c(ccle_mut$genes_muts,
            gdsc_mut$genes_muts
  ),
  type = "mutation_site"
) %>% unique()
tmp$omics_search_mutsites <- tmp$omics_search_mutsites[!grepl("noinfo",tmp$omics_search_mutsites$omics),]

tmp$omics_search_fusion <- data.frame(
  omics = c(ccle_fusion$fusion
  ),
  type = "fusion"
) %>% unique()

omics_search <- rbind(
  tmp$omics_search_CNV,
  tmp$omics_search_mRNA,
  tmp$omics_search_meth,
  tmp$omics_search_proteinrppa,
  tmp$omics_search_proteinms,
  tmp$omics_search_fusion,
  tmp$omics_search_mutgenes,
  tmp$omics_search_mutsites
)
omics_search <- unique(omics_search)

### Drugs----
drugs_search <- data.frame(
  drugs = c(rownames(ctrp1_drug),
            rownames(ctrp2_drug),
            rownames(prism_drug),
            rownames(gdsc1_drug),
            rownames(gdsc2_drug),
            rownames(gCSI_drug),
            rownames(deng1_drug),
            rownames(deng2_drug),
            rownames(deng3_drug)
  ),
  type = "drug"
) %>% unique()

# With each database 
drugs_search2 <- data.frame(
  drugs = c(rownames(ctrp1_drug),
            rownames(ctrp2_drug),
            rownames(prism_drug),
            rownames(gdsc1_drug),
            rownames(gdsc2_drug),
            rownames(gCSI_drug),
            rownames(deng1_drug),
            rownames(deng2_drug),
            rownames(deng3_drug)
  ),
  type = c(
    rep("CTRP1", nrow(ctrp1_drug)),
    rep("CTRP2", nrow(ctrp2_drug)),
    rep("Prism", nrow(prism_drug)),
    rep("GDSC1", nrow(gdsc1_drug)),
    rep("GDSC2", nrow(gdsc2_drug)),
    rep("gCSI", nrow(gCSI_drug)),
    rep("DENG1", nrow(deng1_drug)),
    rep("DENG2", nrow(deng2_drug)),
    rep("DENG3", nrow(deng3_drug))
  )
) %>% unique()

# Combined search 
feas_search <- data.frame(
  name = c(omics_search$omics,
            drugs_search$drugs),
  type = c(omics_search$type,
           drugs_search$type)
)

## Omics data ----
### Continuous ----
tmp$omic_sel <- c("mRNA", "meth", "protein", "cnv")
tmp$tmp1 <- ls()[grepl("_drug$", ls())]
tmp$drug_vec <- gsub("_drug", "", tmp$tmp1[!grepl("^p_", tmp$tmp1)])
omics_search_list1 <- list()

### Discrete ----
# Remake mutation_gene and mutation_site
tmp$tmp2 <- ls()[grepl("_mut$", ls())]
gCSI_mut$mutation <- NA
gCSI_mut$genes_muts <- NA
for(i in tmp$tmp2){
  omic <- base::get(i)
  omic_gene <- omic[,c(1,2)] %>% unique()
  omic_site <- omic[,c(4,2)] %>% unique()
  i2 <- paste0(gsub("mut", "", i), "mutation_gene")
  i3 <- paste0(gsub("mut", "", i), "mutation_site")
  assign(i2, omic_gene)
  assign(i3, omic_site)
}
rm(gCSI_mutation_site)

# prepare search vector ----
fea_list <- list()
fea_vec <- c("mRNA", "meth", 
             "proteinrppa", 
             "proteinms",
             "cnv", # continuous 
             "drug", # drug
             "mutation_gene", "mutation_site", "fusion" # discrete
)
fea_list <- sapply(fea_vec, function(x){
  # Get all objects with the pattern "_x" (e.g., "_drug")
  i2 <- paste0("_", x)
  i2 <- ls(globalenv())[grepl(i2, ls(globalenv()))]
  
  # Filter out objects by name patterns
  i2 <- i2[!grepl("p_", i2)]          # Exclude p_ prefixed objects
  i2 <- i2[!grepl("drugs", i2)]       # Exclude objects with "drugs" in name
  i2 <- i2[!grepl("normalize", i2)]   # Exclude zscore_normalize_drug
  i2 <- i2[!grepl("apply", i2)]       # Exclude apply_zscore functions
  i2 <- i2[!grepl("reset", i2)]       # Exclude reset functions
  i2 <- i2[!grepl("function", i2)]    # Exclude any function-related names
  
  # Only keep known dataset prefixes to ensure we only get actual datasets
  # This is a whitelist approach which is safer
  known_prefixes <- c("ccle", "gdsc", "gCSI", "ctrp1", "ctrp2", "prism", 
                      "gdsc1", "gdsc2",
                     "deng1", "deng2", "deng3")
  
  # Extract the prefix (everything before "_x")
  prefixes <- gsub(paste0("_", x), "", i2)
  
  # Only keep objects with known dataset prefixes
  i2 <- i2[prefixes %in% known_prefixes]
  
  # Extract the dataset name by removing the suffix
  i <- gsub(paste0("_", x), "", i2)
  return(i)
})

# Drugs-omics pairs analysis (New)----
# Create cells search for each omics type for dis compare to dis
tmp$cells_search_CNV <- data.frame(
  cells = c(colnames(ccle_cnv),
            colnames(gdsc_cnv),
            colnames(gCSI_cnv)),
  datasets = c(rep("ccle", ncol(ccle_cnv)),
               rep("gdsc", ncol(gdsc_cnv)),
               rep("gCSI", ncol(gCSI_cnv))),
  type = "cnv"
) %>% unique()

tmp$cells_search_mRNA <- data.frame(
  cells = c(colnames(ccle_mRNA),
            colnames(gdsc_mRNA),
            colnames(deng1_mRNA),
            colnames(deng2_mRNA),
            colnames(deng3_mRNA)),
  datasets = c(rep("ccle", ncol(ccle_mRNA)),
               rep("gdsc", ncol(gdsc_mRNA)),
               rep("deng1", ncol(deng1_mRNA)),
               rep("deng2", ncol(deng2_mRNA)),
               rep("deng3", ncol(deng3_mRNA))),
  type = "mRNA"
) %>% unique()

tmp$cells_search_meth <- data.frame(
  cells = colnames(ccle_meth),
  datasets = rep("ccle", ncol(ccle_meth)),
  type = "meth"
) %>% unique()

tmp$cells_search_proteinrppa <- data.frame(
  cells = colnames(ccle_proteinrppa),
  datasets = rep("ccle", ncol(ccle_proteinrppa)),
  type = "proteinrppa"
) %>% unique()

tmp$cells_search_proteinms <- data.frame(
  cells = colnames(ccle_proteinms),
  datasets = rep("ccle", ncol(ccle_proteinms)),
  type = "proteinms"
) %>% unique()

# For mutation data, assuming the cell lines are stored in a column named 'cells' or similar
tmp$cells_search_mutgenes <- data.frame(
  cells = c(ccle_mutation_gene$cells,
            gdsc_mutation_gene$cells,
            gCSI_mutation_gene$cells),
  datasets = c(rep("ccle", length(ccle_mutation_gene$cells)),
               rep("gdsc", length(gdsc_mutation_gene$cells)),
               rep("gCSI", length(gCSI_mutation_gene$cells))),
  type = "mutation_gene"
) %>% unique()

tmp$cells_search_mutsites <- data.frame(
  cells = c(ccle_mutation_site$cells,
            gdsc_mutation_site$cells),
  datasets = c(rep("ccle", length(ccle_mutation_site$cells)),
               rep("gdsc", length(gdsc_mutation_site$cells))),
  type = "mutation_site"
) %>% unique()

tmp$cells_search_fusion <- data.frame(
  cells = ccle_fusion$cells,
  datasets = rep("ccle", length(ccle_fusion$cells)),
  type = "fusion"
) %>% unique()

# Combine all cells search data frames
cells_search <- rbind(
  tmp$cells_search_CNV,
  tmp$cells_search_mRNA,
  tmp$cells_search_meth,
  tmp$cells_search_proteinrppa,
  tmp$cells_search_proteinms,
  tmp$cells_search_fusion,
  tmp$cells_search_mutgenes,
  tmp$cells_search_mutsites
)
cells_search <- unique(cells_search)


# Plot drug counts and cell subtypes ----
all_stat <- data.frame(
  counts = c(dim(ctrp1_drug),
             dim(ctrp2_drug),
             dim(prism_drug),
             dim(gdsc1_drug),
             dim(gdsc2_drug),
             dim(gCSI_drug),
             dim(deng1_drug),
             dim(deng2_drug),
             dim(deng3_drug)),
  source = rep(c("CTRP1", "CTRP2", "PRISM",
                 "GDSC1", "GDSC2", "gCSI",
                 "DENG1", "DENG2", "DENG3"),
               each = 2),
  type = rep(c("Drugs", "Samples"), times = 9)
)
all_stat$source <- factor(all_stat$source, levels = rev(c("CTRP1", "CTRP2", "PRISM",
                                                          "GDSC1", "GDSC2", "gCSI",
                                                          "DENG3", "DENG2", "DENG1")))
p_count_drugandsample <- ggplot(all_stat, aes(x = source, 
                                            y = counts,
                                            fill = type)) + geom_col(position = "dodge") + 
  geom_text(aes(label = counts), size = 5,
            position = position_dodge(0.9), 
            vjust = -0.8) + theme_bw() + 
  theme(
                                              axis.title.x = element_blank(),
                                              title = element_text(size = 17, color = "black"),
                                              axis.title.y = element_text(size = 17, color = "black"),
                                              axis.text = element_text(size = 17, color = "black"),
                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                              axis.title = element_text(size = 17, color = "black"),
                                              legend.title = element_text(size = 17, colour = "black"),
                                              legend.text = element_text(size = 17),
                                              legend.position = "top"
                                            ) + coord_cartesian(ylim = c(0, 1550))

# Plot counts by data type (PDO vs cell lines) ----
sample_counts_by_type <- table(sample_anno$DataType)

# Create a data frame for plotting
datasource_stat <- data.frame(
  counts = c(
    length(unique(drug_anno[drug_anno$type_source %in% "PDO",]$`Harmonized Name`)),
    length(unique(sample_anno[sample_anno$DataType %in% "PDO",]$SampleID)),
    length(unique(drug_anno[drug_anno$type_source %in% "cell",]$`Harmonized Name`)),
    length(unique(sample_anno[sample_anno$DataType %in% "cell",]$SampleID))
  ),
  source = rep(c("PDO", "Cell"), each = 2),
  type = rep(c("Drugs", "Samples"), times = 2)
)

# Create the plot
p_count_datasource <- ggplot(datasource_stat, aes(x = source, 
                                                  y = counts,
                                                  fill = type)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = counts), size = 5,
            position = position_dodge(0.9), 
            vjust = -0.8) + 
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    title = element_text(size = 17, color = "black"),
    axis.title.y = element_text(size = 17, color = "black"),
    axis.text = element_text(size = 17, color = "black"),
    legend.position = "none"
  ) + 
  ylab("Count") + 
  ggtitle("") + coord_cartesian(ylim = c(0, 2300))

p_count_combined <- patchwork::wrap_plots(p_count_datasource, 
                                          p_count_drugandsample + theme(axis.title.y = element_blank()), 
                                          widths = c(4,10)) + 
  plot_layout(guides = 'keep')

# overlap cell and drug ----
drug_list <- list(
  gdsc1 = rownames(gdsc1_drug),
  gdsc2 = rownames(gdsc2_drug),
  ctrp1 = rownames(ctrp1_drug),
  ctrp2 = rownames(ctrp2_drug),
  prism = rownames(prism_drug),
  gCSI = rownames(gCSI_drug),
  deng1 = rownames(deng1_drug),
  deng2 = rownames(deng2_drug),
  deng3 = rownames(deng3_drug)
)
sample_list <- list(
  gdsc1 = colnames(gdsc1_drug),
  gdsc2 = colnames(gdsc2_drug),
  ctrp1 = colnames(ctrp1_drug),
  ctrp2 = colnames(ctrp2_drug),
  prism = colnames(prism_drug),
  gCSI = colnames(gCSI_drug),
  deng1 = colnames(deng1_drug),
  deng2 = colnames(deng2_drug),
  deng3 = colnames(deng3_drug)
)
p_overlap_sample <- upset(fromList(sample_list), mainbar.y.label = "Sample Counts", text.scale = 2,
                        nsets = length(drug_list))
p_overlap_drug <- upset(fromList(drug_list), mainbar.y.label = "Drug counts", text.scale = 2,
                        nsets = length(drug_list))

# Plot molecular characteristics ----
# First, let's create a data frame from your list
tmp$fea_list2 <- fea_list
tmp$fea_list2$drug <- NULL
tmp$data_types <- names(tmp$fea_list2)
tmp$databases <- unique(unlist(tmp$fea_list2))

# Calculate overlap count for each feature
tmp$overlap_counts <- sapply(tmp$fea_list2, length)

# Create a new order based on overlap count (descending) and desired priority
tmp$feature_order <- names(sort(tmp$overlap_counts, decreasing = TRUE))

# Create a matrix of presence/absence
tmp$matrix_data <- expand.grid(
  Feature = tmp$data_types,
  Database = tmp$databases
) %>%
  mutate(Present = mapply(function(feat, db) {
    db %in% tmp$fea_list2[[feat]]
  }, Feature, Database)) %>%
  # Reorder features to match your desired order
  mutate(Feature = factor(Feature, 
                          levels = tmp$feature_order))

# Create mapping between your names and display names
tmp$name_mapping <- c(
  "mRNA" = "Gene Expression",
  "cnv" = "Copy Number",
  "mutation_gene" = "Gene Mutation",
  "mutation_site" = "Gene Site Mutation",
  "fusion" = "Gene Fusion",
  # "drug" = "Chromatin Profile",
  "meth" = "DNA Methylation",
  "proteinrppa" = "Proteome RPPA Expression",
  "proteinms" = "Proteome MS Expression"
)

# Rest of the code remains the same, but update Feature2 ordering
tmp$matrix_data$Feature2 <- unname(tmp$name_mapping[match(tmp$matrix_data$Feature, names(tmp$name_mapping))])
tmp$matrix_data$Feature2 <- factor(tmp$matrix_data$Feature2, 
                                   levels = unique(tmp$name_mapping[tmp$feature_order]))
tmp$matrix_data$Database <- factor(tmp$matrix_data$Database,
                                   levels = rev(c( "ccle", "gdsc", "gCSI", 
                                                  "deng3", "deng2", "deng1")))
# Create the plot
p_mol_character <- ggplot(tmp$matrix_data, aes(x = Feature2, y = Database)) +
  geom_point(aes(size = Present, color = Present)) +
  scale_size_manual(values = c(1, 4)) +
  scale_color_manual(values = c("grey80", "black")) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 17, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )  + labs(x = "", y = "") + # coord_flip() +  
  coord_fixed(ratio = 1) 

# Plot drug MOA Freq ----
drug_MOA_freq <- data.frame(table(drug_anno$MOA))
drug_moa <- drug_MOA_freq[!drug_MOA_freq$Var1 %in% c("0",""),]
all_targets <- c()
# Loop through each row in the dataframe
for (i in 1:nrow(drug_moa)) {
  # Split the compound target by "|"
  targets <- unlist(strsplit(as.character(drug_moa$Var1[i]), "\\|"))
  
  # Repeat each target by its frequency in the original data
  repeated_targets <- rep(targets, drug_moa$Freq[i])
  
  # Add to the all_targets vector
  all_targets <- c(all_targets, repeated_targets)
}

# Create a frequency table of individual targets
target_freq <- table(all_targets)

# Convert to a dataframe
target_freq_df <- data.frame(
  Target = names(target_freq),
  Frequency = as.numeric(target_freq)
)

# Sort by frequency in descending order
drug_MOA_freq <- target_freq_df[order(-target_freq_df$Frequency), ]
top_targets <- drug_MOA_freq %>%
  filter(Frequency >= 10)

# Create a category for targets with lower frequencies
# other_targets <- drug_MOA_freq %>%
#   filter(Frequency < 5) %>%
#   summarize(Target = "Other targets", Frequency = sum(Frequency))

# Combine the datasets
# plot_data <- rbind(top_targets, other_targets)

# Create the treemap
p_drug_moa <- ggplot(top_targets, aes(area = Frequency, fill = Frequency, label = Target)) +
  geom_treemap() +
  geom_treemap_text(
    aes(label = paste0(Target, "\n(", Frequency, ")")),
    colour = "white",
    place = "centre",
    size = 12,
    grow = TRUE
  ) +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  theme_minimal() +
  theme(
    legend.position = "none",
    title = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 17, colour = "black"),
    legend.text = element_text(size = 17)
  ) +
  labs(
    title = "Drug MOA Frequency Distribution"
  )

# Plot sample tumor type ----
tumor_data <- as.data.frame(
  table(sample_anno$TumorType)
)
colnames(tumor_data) <- c("TumorType", "Frequency")


# Group tumor types by system
tumor_data$System <- case_when(
  grepl("lung|aerodigestive|nasopharyngeal", tumor_data$TumorType) ~ "Respiratory",
  grepl("gastrointestinal|stomach|liver|pancreatic", tumor_data$TumorType) ~ "Digestive",
  grepl("breast|ovarian|cervical|endometrial|uterine|vulvar", tumor_data$TumorType) ~ "Reproductive - Female",
  grepl("prostate|testicular", tumor_data$TumorType) ~ "Reproductive - Male",
  grepl("haematopoietic|lymphoid", tumor_data$TumorType) ~ "Blood & Lymphatic",
  grepl("nervous", tumor_data$TumorType) ~ "Nervous System",
  grepl("skin", tumor_data$TumorType) ~ "Integumentary",
  grepl("kidney|bladder", tumor_data$TumorType) ~ "Urinary",
  grepl("sarcoma", tumor_data$TumorType) ~ "Connective Tissue",
  TRUE ~ "Other"
)

# Shorten names for better display
tumor_data$TumorType <- gsub(" cancer", "", tumor_data$TumorType)

# Create the bubble chart
p_tumor_bubble <- ggplot(tumor_data, aes(x = System, y = Frequency, size = Frequency, color = System)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    aes(label = paste0(TumorType, " (", Frequency, ")")
    ),
    color = "black",
    size = 3.2,
    box.padding = 0.7,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  scale_size(range = c(3, 15)) +
  scale_color_brewer(palette = "Set3") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text = element_text(size = 17, color = "black"),
    # panel.grid = element_line(color = "black"),
    legend.position = "none",
    plot.title = element_text(size = 17, color = "black", hjust = 0.5)
  ) +
  labs(
    title = "Tumor Type Distribution by System",
    x = "",
    y = "Number of Samples"
  )



# Others ----
