tmp <- list()

# Make copy ----
# PDX
tavor_drug_raw = tavor_drug
PDTXBreast_drug_raw = PDTXBreast_drug
# CellLine
ccle_drug_raw = ccle_drug
ctrp1_drug_raw = ctrp1_drug
ctrp2_drug_raw = ctrp2_drug
gdsc1_drug_raw = gdsc1_drug
gdsc2_drug_raw = gdsc2_drug
gCSI_drug_raw = gCSI_drug
prism_drug_raw = prism_drug
FIMM_drug_raw = FIMM_drug
UHNBreast_drug_raw = UHNBreast_drug
GRAY_drug_raw = GRAY_drug
NCI60_drug_raw = NCI60_drug
# PDO
UMPDO1_drug_raw = UMPDO1_drug
UMPDO2_drug_raw = UMPDO2_drug
UMPDO3_drug_raw = UMPDO3_drug
# PDX
Xeva_drug_raw = Xeva_drug

# Preprocess ----
### Discrete ----
# Remake mutation_gene and mutation_site
# Use the function with your variables
# mutation_variables <- c("ccle_mut", "gCSI_mut", "gdsc_mut", "Xeva_mut")
tmp$tmp2 <- ls()[grepl("_mut$", ls())]
process_mutation_data(tmp$tmp2)

## prepare search vector ----
fea_list <- list()
fea_vec <- c("mRNA", "meth", 
             "proteinrppa", 
             "proteinms",
             "cnv", # continuous 
             "drug", # drug
             "drug_raw", # drug_raw
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
                      "FIMM", "GRAY", "NCI60", "UHNBreast", 
                      "tavor", "PDTXBreast",
                      "UMPDO1", "UMPDO2", "UMPDO3", "Xeva")
  
  # Extract the prefix (everything before "_x")
  prefixes <- gsub(paste0("_", x), "", i2)
  
  # Only keep objects with known dataset prefixes
  i2 <- i2[prefixes %in% known_prefixes]
  
  # Extract the dataset name by removing the suffix
  i <- gsub(paste0("_", x), "", i2)
  return(i)
})


## Search ----
### Omics ----
tmp$omics_search_CNV <- data.frame(
  omics = c(rownames(ccle_cnv),
            rownames(gdsc_cnv),
            rownames(gCSI_cnv),
            rownames(Xeva_cnv)
  ),
  type = "cnv"
) %>% unique()


tmp$omics_search_mRNA <- data.frame(
  omics = c(rownames(ccle_mRNA),
            rownames(gdsc_mRNA),
            rownames(NCI60_mRNA),
            rownames(tavor_mRNA),
            rownames(UMPDO1_mRNA),
            rownames(UMPDO2_mRNA),
            rownames(UMPDO3_mRNA),
            rownames(Xeva_mRNA)
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
  omics = c(ccle_mutation_gene$genes,
            gdsc_mutation_gene$genes,
            gCSI_mutation_gene$genes,
            Xeva_mutation_gene$genes
  ),
  type = "mutation_gene"
) %>% unique()

tmp$omics_search_mutsites <- data.frame(
  omics = c(ccle_mutation_site$genes_muts,
            gdsc_mutation_site$genes_muts
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
# drugs_search2 <- unique(drug_anno[,c(1,2)])
# colnames(drugs_search2) <- c("drugs", "type")

drugs_search <- unique(
  data.frame(
    drugs = unique(drug_anno$DrugName),
    type = "drug"
  )
)

# Combined search 
feas_search <- data.frame(
  name = c(omics_search$omics,
            drugs_search$drugs),
  type = c(omics_search$type,
           drugs_search$type)
)
### Sample ----
# Create data frames for CNV data
tmp$cells_search_CNV <- data.frame(
  cells = c(colnames(ccle_cnv),
            colnames(gdsc_cnv),
            colnames(gCSI_cnv),
            colnames(Xeva_cnv)),  # Added Xeva_cnv
  datasets = c(rep("ccle", ncol(ccle_cnv)),
               rep("gdsc", ncol(gdsc_cnv)),
               rep("gCSI", ncol(gCSI_cnv)),
               rep("Xeva", ncol(Xeva_cnv))),  # Added Xeva_cnv
  type = "cnv"
) %>% unique()

# Create data frames for mRNA data
tmp$cells_search_mRNA <- data.frame(
  cells = c(colnames(ccle_mRNA),
            colnames(gdsc_mRNA),
            colnames(NCI60_mRNA),    # Added NCI60_mRNA
            colnames(tavor_mRNA),    # Added tavor_mRNA
            colnames(UMPDO1_mRNA),
            colnames(UMPDO2_mRNA),
            colnames(UMPDO3_mRNA),
            colnames(Xeva_mRNA)),    # Added Xeva_mRNA
  datasets = c(rep("ccle", ncol(ccle_mRNA)),
               rep("gdsc", ncol(gdsc_mRNA)),
               rep("NCI60", ncol(NCI60_mRNA)),    # Added NCI60_mRNA
               rep("tavor", ncol(tavor_mRNA)),    # Added tavor_mRNA
               rep("deng1", ncol(UMPDO1_mRNA)),
               rep("deng2", ncol(UMPDO2_mRNA)),
               rep("deng3", ncol(UMPDO3_mRNA)),
               rep("Xeva", ncol(Xeva_mRNA))),    # Added Xeva_mRNA
  type = "mRNA"
) %>% unique()

# Create data frames for methylation data
tmp$cells_search_meth <- data.frame(
  cells = colnames(ccle_meth),
  datasets = rep("ccle", ncol(ccle_meth)),
  type = "meth"
) %>% unique()

# Create data frames for protein RPPA data
tmp$cells_search_proteinrppa <- data.frame(
  cells = colnames(ccle_proteinrppa),
  datasets = rep("ccle", ncol(ccle_proteinrppa)),
  type = "proteinrppa"
) %>% unique()

# Create data frames for protein MS data
tmp$cells_search_proteinms <- data.frame(
  cells = colnames(ccle_proteinms),
  datasets = rep("ccle", ncol(ccle_proteinms)),
  type = "proteinms"
) %>% unique()

# For mutation gene data
tmp$cells_search_mutgenes <- data.frame(
  cells = c(ccle_mutation_gene$cells,
            gdsc_mutation_gene$cells,
            gCSI_mutation_gene$cells,
            Xeva_mutation_gene$cells),  # Added Xeva_mutation_gene
  datasets = c(rep("ccle", length(ccle_mutation_gene$cells)),
               rep("gdsc", length(gdsc_mutation_gene$cells)),
               rep("gCSI", length(gCSI_mutation_gene$cells)),
               rep("Xeva", length(Xeva_mutation_gene$cells))),  # Added Xeva_mutation_gene
  type = "mutation_gene"
) %>% unique()

# For mutation site data
tmp$cells_search_mutsites <- data.frame(
  cells = c(ccle_mutation_site$cells,
            gdsc_mutation_site$cells),
  datasets = c(rep("ccle", length(ccle_mutation_site$cells)),
               rep("gdsc", length(gdsc_mutation_site$cells))),
  type = "mutation_site"
) %>% unique()

# For fusion data
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
# Get all variables ending with "_drug"
drug_vars <- ls(pattern = "_drug$")
drug_vars <- drug_vars[!grepl("normalize", drug_vars)]   # Exclude zscore_normalize_drug

# Define project groups
project_groups <- list(
  "PDC" = c("tavor_drug", "PDTXBreast_drug"),
  "CellLine" = c("ctrp1_drug", "ctrp2_drug", "gdsc1_drug", "gdsc2_drug", "gCSI_drug", 
                 "prism_drug", "FIMM_drug", "UHNBreast_drug", "GRAY_drug", "NCI60_drug", "ccle_drug"),
  "PDO" = c("UMPDO1_drug", "UMPDO2_drug", "UMPDO3_drug"),
  "PDX" = c("Xeva_drug")
)

# Initialize empty data frame
all_stat <- data.frame(
  counts = numeric(),
  source = character(),
  type = character(),
  group = character()
)

# Populate data frame dynamically
for (var in drug_vars) {
  # Get the object
  drug_obj <- base::get(var)
  
  # Extract source name (remove "_drug" suffix)
  source_name <- sub("_drug$", "", var)
  
  # Determine group
  group_name <- "Other"
  for (g in names(project_groups)) {
    if (var %in% project_groups[[g]]) {
      group_name <- g
      break
    }
  }
  
  # Add dimensions to data frame
  all_stat <- rbind(all_stat, data.frame(
    counts = dim(drug_obj),
    source = rep(source_name, 2),
    type = c("Drugs", "Samples"),
    group = rep(group_name, 2)
  ))
}
# Calculate total counts (sum of Drugs and Samples) for each source
source_totals <- aggregate(counts ~ source + group, data = all_stat, FUN = sum)

# Sort by group and then by total counts within each group (ascending)
all_stat$group <- factor(all_stat$group, levels = c("CellLine", "PDC", "PDO", "PDX"))

# Create a lookup table for the total counts by source
source_order <- source_totals[order(source_totals$group, source_totals$counts),] # Removed the minus sign
source_order$rank <- ave(source_order$counts, source_order$group, 
                         FUN = function(x) rank(x, ties.method = "first")) # Removed the minus sign

# Join this rank back to the original data
all_stat <- base::merge(all_stat, source_order[,c("source", "rank")], by = "source")

# Order by group first, then by the rank within group
all_stat <- all_stat[order(all_stat$group, all_stat$rank, all_stat$type),]

all_stat$source <- factor(all_stat$source, levels = unique(all_stat$source))

all_stat$group_id <- as.numeric(all_stat$group)
all_stat$group_id <- all_stat$group_id %% 2  # Alternating 0 and 1 for different background colors

# Create a data frame for group labels
# First create the gapped plot without group labels
p_count_drugandsample_bg_no_labels <- ggplot(all_stat, aes(x = source, 
                                                           y = counts,
                                                           fill = type)) + 
  # Add alternating background
  geom_rect(data = subset(all_stat, group_id == 1),
            aes(xmin = as.numeric(source) - 0.5, 
                xmax = as.numeric(source) + 0.5,
                ymin = -Inf, ymax = Inf),
            fill = "gray90", alpha = 0.3,
            inherit.aes = FALSE) +
  # Add columns  
  geom_col(position = "dodge") + 
  # Add text labels
  geom_text(aes(label = counts), size = 4.5,
            position = position_dodge(0.9), 
            vjust = -0.5) +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    title = element_text(size = 17, color = "black"),
    axis.title.y = element_text(size = 17, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
    axis.title = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 17, colour = "black"),
    legend.text = element_text(size = 17),
    legend.position = "top"
  )

# Apply gg.gap
p_count_drugandsample_facet_with_gap <- gg.gap(
  plot = p_count_drugandsample_bg_no_labels,
  ylim = c(0, max(all_stat$counts) * 1.05),
  segments = list(c(2000, 54000)),
  tick_width = c(500, 5000),
  rel_heights = c(0.7, 0.05, 0.25)
)

# Plot all counts ----
datasource_stat <- data.frame(
  counts = c(length(unique(sample_anno$SampleID)),
             length(unique(drug_anno$DrugName))),
  type = c("Samples", "Drugs")
)

p_count_drugandsample_sum <- ggplot(datasource_stat, aes(x = type, 
                                                  y = counts,
                                                  fill = type)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = counts), size = 5,
            position = position_dodge(0.9), 
            vjust = -0.8) + 
  theme_bw() + 
  guides(fill = guide_legend(nrow = 1)) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    title = element_text(size = 17, color = "black"),
    axis.title.y = element_text(size = 17, color = "black"),
    axis.text = element_text(size = 17, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    # legend.position = "top"
    ) + 
  ylab("Count") + 
  ggtitle("") 

legend <- cowplot::get_legend(p_count_drugandsample_sum)

p_count_drugandsample_sum <- p_count_drugandsample_sum + theme(
  legend.position = "none"
)

p_count_drugandsample_sum_with_gap <- gg.gap(
  plot = p_count_drugandsample_sum + theme(legend.position = "none"),
  ylim = c(0, 70000),
  segments = list(c(4500, 55000)),
  tick_width = c(500, 5000),
  rel_heights = c(0.7, 0.05, 0.25)
)

p_count_drugandsample_sum_with_gap2 <- cowplot::plot_grid(
  legend,
  p_count_drugandsample_sum_with_gap,
  ncol = 1,
  rel_heights = c(0.2, 0.9)
)

# p_count_combined <- patchwork::wrap_plots(p_count_drugandsample_facet_with_gap, 
#                                           p_count_drugandsample_sum_with_gap + theme(axis.title.y = element_blank()), 
#                                           widths = c(4,10)) + 
#   plot_layout(guides = 'keep')

# overlap cell and drug ----
# Create lists of drugs and samples for all datasets
drug_list <- list(
  # PDC
  tavor = rownames(tavor_drug),
  PDTXBreast = rownames(PDTXBreast_drug),
  # CellLine
  ccle = rownames(ccle_drug),
  ctrp1 = rownames(ctrp1_drug),
  ctrp2 = rownames(ctrp2_drug),
  gdsc1 = rownames(gdsc1_drug),
  gdsc2 = rownames(gdsc2_drug),
  gCSI = rownames(gCSI_drug),
  prism = rownames(prism_drug),
  FIMM = rownames(FIMM_drug),
  UHNBreast = rownames(UHNBreast_drug),
  GRAY = rownames(GRAY_drug),
  NCI60 = rownames(NCI60_drug),
  # PDO
  UMPDO1 = rownames(UMPDO1_drug),
  UMPDO2 = rownames(UMPDO2_drug),
  UMPDO3 = rownames(UMPDO3_drug),
  # PDX
  Xeva = rownames(Xeva_drug)
)

sample_list <- list(
  # PDC
  tavor = colnames(tavor_drug),
  PDTXBreast = colnames(PDTXBreast_drug),
  # CellLine
  ccle = colnames(ccle_drug),
  ctrp1 = colnames(ctrp1_drug),
  ctrp2 = colnames(ctrp2_drug),
  gdsc1 = colnames(gdsc1_drug),
  gdsc2 = colnames(gdsc2_drug),
  gCSI = colnames(gCSI_drug),
  prism = colnames(prism_drug),
  FIMM = colnames(FIMM_drug),
  UHNBreast = colnames(UHNBreast_drug),
  GRAY = colnames(GRAY_drug),
  NCI60 = colnames(NCI60_drug),
  # PDO
  UMPDO1 = colnames(UMPDO1_drug),
  UMPDO2 = colnames(UMPDO2_drug),
  UMPDO3 = colnames(UMPDO3_drug),
  # PDX
  Xeva = colnames(Xeva_drug)
)
p_overlap_sample <- upset(fromList(sample_list), mainbar.y.label = "Sample Counts", text.scale = 2,
                        nsets = length(drug_list))
p_overlap_drug <- upset(fromList(drug_list), mainbar.y.label = "Drug counts", text.scale = 2,
                        nsets = length(drug_list))

# Plot molecular characteristics ----
# First, let's create a data frame from your list
tmp$fea_list2 <- fea_list
tmp$fea_list2$drug <- NULL
tmp$fea_list2$drug_raw <- NULL
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
                                   levels = rev(c( "ccle", "gdsc","Xeva", "gCSI",
                                                   "tavor", "NCI60",
                                                  "UMPDO1", "UMPDO2", "UMPDO3")))
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
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
