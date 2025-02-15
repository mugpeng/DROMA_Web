# DROMA
 Drug Response Omics Association Map.



DROMA-DB: A database that includes the largest published studies(PDC, PDO, PDX) investigating the cancer to chemical compound treatment, and the association between drug sensitivity and multi-omics(mRNA, CNV, protein, mutation, etc.).

![image-20250121100802484](images/image-20250121100802484.png)



Please cite when you used in your study:

Li, S., Peng, Y., Chen, M. et al. Facilitating integrative and personalized oncology omics analysis with UCSCXenaShiny. Commun Biol 7, 1200 (2024). https://doi.org/10.1038/s42003-024-06891-2



# GOALS

A database (DROMA-DB)
1) includes high-throughput cancer type PDO, PDC, PDX data associated between drug sensitivity and multi-omics, also Denglab in house data; 
2) Other useful gene-drug correlation information from FDA and papers.
3) Provide powerful tools for finding potential drug targets, repurposing drugs for target genes. 
4) Friendly methods to predict the patients treatment by providing omics data.
5) Supply both code package and web interface for people to use.



An AI (DROMA-AI)

1) Talk to her all your interests about genes or drugs like an omniscient friend/tool based on DROMA-DB and other training data.
2) Ask her to use any functions from DROMA-DB.



# TODO





# Tutorial

![](http://cos01.mugpeng.top/img/20250121100754.png)

The website consisted of three main sections (pages): 1) Drugs-omics Pairs Analysis; 2) Batch Features Associations Analysis; and 3) Statistics Information, which is the final page. The tour will start from this last page.



Following the principle of least change, AAC values from gCSI are rescaled such that a lower metric indicates higher sensitivity across all datasets. (AAC2 = max(AAC) - AAC)

![image-20240116150920858](http://cos01.mugpeng.top/img/image-20240116150920858.png)

Lower scores on metric scales mean stronger drug sensitivity.



## Statistics Information

The omics data we collected including mRNA expression, protein data, copy number variant(CNV), gene fusion, methylation, mutation(gene mutation, and a specific amino acid change).



In vitro cell line pharmacogenetics studies can be categorized into CCLE, GDSC, and other projects (produced by different insititutions). CCLE, GDSC projects produce the genomics data, then other experiments are conducted to generate drug data used the same cell lines one of the projects.

The cell lines from the Cancer Cell Line Encyclopedia (CCLE) were used to generate several drug sensitivity projects, including CTRP1, CTRP2, and PRISM. Meanwhile, the Genomics of Drug Sensitivity in Cancer (GDSC) project generated GDSC1 and GDSC2. And Genentech Cell Screening Initiative (gCSI) project has its own omics and drug response data.



The PRISM project tests the highest number of drugs, while GDSC1 focuses on testing the largest number of cells.

![image-20240116134858293](http://cos01.mugpeng.top/img/image-20240116134858293.png)

Projects have tested same drugs and same cells. When multiple projects test the same drugs and cells, there is typically a higher degree of overlap between projects that use the same cell lines.

![image-20240116140858989](http://cos01.mugpeng.top/img/image-20240116140858989.png)

![image-20240116140840856](http://cos01.mugpeng.top/img/image-20240116140840856.png)

The cell lines are mainly from lung cancer, colorectal cancer, and ovarian cancer:

![image-20240116141100546](http://cos01.mugpeng.top/img/image-20240116141100546.png)

Also there is annotation for cell and drug:

![image-20240116143433413](http://cos01.mugpeng.top/img/image-20240116143433413.png)

![image-20240116143440528](http://cos01.mugpeng.top/img/image-20240116143440528.png)



You can search drugs of interests.

![image-20240116143707243](http://cos01.mugpeng.top/img/image-20240116143707243.png)



## Main Function

### 1) Drugs-omics pairs analysis

This feature allowed user to explore the association between a selected drug resistance event and a certain omic. For continuous omics data like mRNA, methylation, copy number variant, protein, spearman correlation was calculated. While for discrete omics data such as mutation genes, mutation gene points or gene fusions, wilcoxon test is chosen for testing Signification.

![](http://cos01.mugpeng.top/img/20250121101140.png)



The title of each plot indicates the source of the omics and drug response data. For example, a plot titled `gdsc_ctrp1` would mean the omics data is from the GDSC project, while the drug sensitivity data is from the CTRP1 project, as mentioned in the initial *Statistics Information* section. Personally, I think comparing cells data from different organizations (e.g. GDSC vs CCLE) is reasonable for analyzing correlations, as we are primarily interested in examining the relationship between omic features and drug responses, regardless of the original data source. Combining data from multiple sources can provide a more comprehensive view of these relationships.

The upper figure shows that the gene expression of ABCC3 is significantly positively correlated with sensitivity to the drug YM-155 across all ten dataset combinations. This correlation could potentially be explained by the known function of the ABCC3 gene. Specifically, ABCC3 encodes an ATP-binding cassette transporter protein that is involved in exporting various molecules, including drugs, out of cells via active transport across cell membranes. Given its role in drug efflux, higher ABCC3 expression may correlate with increased efflux and reduced intracellular concentration of YM-155, resulting in greater resistance to the drug's effects. This might explain the observed positive correlation between ABCC3 expression levels and higher values for YM-155.



Another example is gene mutation of TP53 and drug AMG-232, it is obvious that the wild type TP53 has significant higher sensitivity:

![](http://cos01.mugpeng.top/img/20250121102231.png)

AMG-232 is an inhibitor of the p53-MDM2 interaction. Mutated tp53 may deactivate the suppressor program induced by AMG-232 through disruption of this interaction: [p53-family proteins and their regulators: hubs and spokes in tumor suppression | Cell Death & Differentiation (nature.com)](https://www.nature.com/articles/cdd201035)



This function displays a meta-analysis forest plot and a download option. The plot summarizes the results of multiple studies, with each row representing an individual study's effect size and confidence interval. The diamond at the bottom represents the overall pooled effect size (0.27, 95% CI: 0.23-0.30) across all studies. The highly significant p-value (p = 6.81e-53) confirms a strong overall effect. Because the confidence interval does not include zero, the effect is statistically significant. The plot also outputs an overall effect size (similar to an R coefficient) and a p-value. 

The direction of the overall effect is consistent with the individual study effect sizes (e.g., R coefficients in the mRNA example). A positive effect size indicates that higher gene expression is associated with drug resistance across studies, suggesting a potential drug resistance gene. Conversely, a negative effect size would indicate that lower gene expression (or wild type) is associated with reduced drug sensitivity, suggesting that the presence of the gene may enhance drug efficacy:

![](http://cos01.mugpeng.top/img/20250121101738.png)

You can save both the each test plot and summarized forest plot, or save the data used to plot to draw a new one by yourself.



### 2) Batch Features Associations Analysis

This analysis module helped people to conduct a significant test between a targeted feature(a drug or an omic) and all the features in a particular feature dataset grouped by their collected databases in a large scale.

The statistic test is calculated depend on the data types:

- For continuous features compared to continuous datasets (e.g. drug A levels versus all CNV features), the Pearson correlation test is used.
- For discrete features compared to discrete databases (e.g. TP53 mutation events versus all collected gene fusions), will use Wilcoxon test.
- For discrete features compared to discrete databases, Chi-squared test is the choice.

All the calculated results will be summarize by meta function with `metagen` in R package "meta". And produce a volcano plot:

![](http://cos01.mugpeng.top/img/20250121105624.png)

[Proteasome inhibitors block Ikaros degradation by lenalidomide in multiple myeloma - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5004433/)

![](http://cos01.mugpeng.top/img/20250121113016.png)



A feature-database pair will be considered statistically significant if both of the following criteria are met:

1. The absolute value of the effect size is greater than 0.2. 
2. The p-value is less than 0.001.
3. Top 5 in up and down will be marked.







We will find the potential related mRNA with Lapatinib as an example:

![image-20240116172531752](http://cos01.mugpeng.top/img/image-20240116172531752.png)



Frequency table has two columns, frequency col counts the number of pairs labeled as significant in all databases containing this pair. Proportion col is the fraction of significant pair in all pairs. You can choose the topmost to further examination with result table.

![image-20240117161051412](http://cos01.mugpeng.top/img/image-20240117161051412.png)

For example, we are interested about CDH1 gene, which is a classical cadherin of the cadherin superfamily, Mutations in this gene are correlated with gastric, breast, colorectal, thyroid and ovarian cancer. Loss of function of this gene is thought to contribute to cancer progression by increasing proliferation, invasion, and/or metastasis, described from genecode: [CDH1 Gene - GeneCards | CADH1 Protein | CADH1 Antibody](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CDH1&keywords=CDH1)

We can search this gene by the search box at the top right edge. The results indicate that higher expression of CDH1 is correlated with increased sensitivity to LAPATINIB. As drug resistance metrics are negatively correlated with drug sensitivity, a higher CDH1 expression level tends to predict greater LAPATINIB sensitivity.

![image-20240117161119979](http://cos01.mugpeng.top/img/image-20240117161119979.png)

A downloadable button is also provided, allowing users to access a CSV file containing the data. This CSV file can then be used for additional analyses or data processing as needed.

![image-20240117161344003](http://cos01.mugpeng.top/img/image-20240117161344003.png)



A simple online search reveals that CDH1 is related to ERBB2, and ERBB2 is a validated target of LAPATINIB. This suggests that CDH1 expression levels may help determine a tumor's response to LAPATINIB treatment, potentially through its relationship to the drug's primary target, ERBB2. 

![image-20240117160239539](http://cos01.mugpeng.top/img/image-20240117160239539.png)



By the way, we can also double check it through *Drugs-omics pairs analysis* module:

![image-20240117160637855](http://cos01.mugpeng.top/img/image-20240117160637855.png)



Besides, I also design some useful widgets.

A bar to indicate the progress:

![](http://cos01.mugpeng.top/img/20250121104517.png)

A jump box tell the analysis has finished:

![](http://cos01.mugpeng.top/img/20250121104503.png)



# Contact with me

Feel free to talk with me if you find any bugs or have any suggestions. :)

Email: mugpeng@foxmail.com, mugpeng@outlook.com, yc47680@um.edu.mo
