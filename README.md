# CAUTION
A version with a better design, cleaner structure, and AI enable will be release soon.
This repository will **NOT** be updated and maintained anymore. See you in the new version :). Bye. Peng, 20251117

[![Website](https://img.shields.io/website?url=https%3A//droma01.github.io/DROMA/)](https://droma01.github.io/DROMA/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-1.0+-green.svg)](https://shiny.rstudio.com/)
[![Version](https://img.shields.io/badge/version-0.3.0-blue.svg)](https://github.com/mugpeng/DROMA_Web)

# DROMA Web

A **Shiny web application** for **DROMA** (Drug Response Omics association MAp) - providing an interactive interface for cancer drug sensitivity and multi-omics association analysis.

[![image](https://github.com/user-attachments/assets/89b385e0-f5e4-4d8e-a33b-5f30bed039b2)](https://droma01.github.io/DROMA/)

It is a part of [DROMA project](https://github.com/mugpeng/DROMA). Visit the [official DROMA website](https://droma01.github.io/DROMA/) for comprehensive documentation and interactive examples.

## 🚀 Overview

DROMA Web provides a user-friendly web interface for analyzing drug-omics associations through integration with [DROMA.R](https://github.com/mugpeng/DROMA_R) and [DROMA.Set](https://github.com/mugpeng/DROMA_Set) packages.

### Key Features

- **🔗 Interactive Analysis**: Analyze drug-omics associations through an intuitive web interface
- **📊 Real-time Visualization**: Generate correlation plots, volcano plots, and forest plots on demand
- **📈 Drug-Omic Pair Analysis**: Explore relationships between drugs and molecular features (mRNA, CNV, methylation, mutations)
- **🗂️ Batch Feature Analysis**: Perform high-throughput association studies across multiple features
- **💾 Data Export**: Download results in PDF, CSV, and R formats
- **⚡ Parallel Processing**: Efficient computation with multi-core support
- **🔄 Z-score Normalization**: Real-time data normalization toggle
- **🎯 Multiple Data Types**: Support various data types (CellLine, PDC, PDO, PDX)

## 📦 Installation

### Prerequisites

- R 4.0+
- DROMA.Set and DROMA.R packages
- DROMA SQLite database (`droma.sqlite`)

### Required R Packages

```r
# Core packages
install.packages(c("shiny", "shinyWidgets", "shinyjs", "waiter", "DT"))

# Data manipulation
install.packages(c("dplyr", "data.table"))

# Statistical analysis
install.packages(c("meta", "metafor", "effsize"))

# Visualization
install.packages(c("ggplot2", "UpSetR", "ggpubr", "ggrepel",
                   "treemapify", "gridExtra", "patchwork", "cowplot"))

# Parallel processing
install.packages(c("future", "furrr", "snowfall", "parallel"))
```

### Install DROMA Packages

```r
# from GitHub
devtools::install_github("mugpeng/DROMA_Set")
devtools::install_github("mugpeng/DROMA_R")

# Install from local source
install.packages("../DROMA_Set", repos = NULL, type = "source")
install.packages("../DROMA_R", repos = NULL, type = "source")
```

## 🚀 Quick Start

### 1. Setup Database

Place your `droma.sqlite` database file in the `data/` directory:

```bash
DROMA_Web/
├── App.R
├── Modules/
├── data/
│   └── droma.sqlite  # <-- Place your database here
└── ...
```

### 2. Run the Application

```r
# Method 1: Using RStudio
# - Open the project in RStudio
# - Click "Run App" button

# Method 2: Using R console
setwd("path/to/DROMA_Web")
shiny::runApp()

# Method 3: Specify host and port
shiny::runApp(host = "0.0.0.0", port = 8080)
```

### 3. Access the Web Interface

Open your web browser and navigate to:
- Local: `http://localhost:3838` (default)
- Custom: `http://your-server:port` (if specified)

## 🖥️ Application Features

### Main Analysis Modules

1. **Drugs-omics Pairs Analysis**
   - Select drug and molecular feature pairs
   - Generate correlation and association plots
   - Perform meta-analysis across multiple studies

2. **Batch Features Associations Analysis**
   - High-throughput association screening
   - Volcano plot visualization
   - Export significant associations

3. **Drug Feature Analysis**
   - Explore drug sensitivity distributions
   - Compare across tumor types and sample attributes
   - Download annotated drug data

4. **Statistics and Annotations**
   - Dataset overview and statistics
   - Sample and drug annotation tables
   - Download annotation data

### Key Controls

- **Global Settings** (floating button on left):
  - Toggle Z-score normalization on/off
  - Affects all analyses in real-time

- **Data Filters**:
  - Filter by data type (CellLine, PDC, PDO, PDX)
  - Filter by tumor type

## 📊 Data Requirements

Your `droma.sqlite` database should contain:

- **Molecular Profiles**: mRNA, CNV, methylation, protein expression
- **Drug Response**: IC50, AUC, or other sensitivity metrics
- **Sample Annotations**: Tumor type, tissue origin, etc.
- **Drug Annotations**: Mechanism of action, target information

## 🔧 Configuration

### Environment Setup

The application uses a fixed database path: `data/droma.sqlite`

For custom database locations, modify `Modules/DataAdapter.R`:

```r
# In initializeDROMAData() function
db_path <- "your/custom/path/droma.sqlite"
```

### Performance Tuning

Adjust parallel processing in batch analysis:

```r
# In Modules/BatchFeature.R
used_core <- parallel::detectCores() - 1  # Use all but one core
```

## 🛠️ Troubleshooting

### Common Issues

**Database Not Found**:
```
Error: DROMA database not found at: data/droma.sqlite
```
- Solution: Ensure `droma.sqlite` is in the `data/` directory

**Package Loading Errors**:
```
Error: there is no package called 'DROMA.Set'
```
- Solution: Install DROMA packages before running

**Memory Issues**:
- Solution: Close unused datasets in browser tabs
- Consider using a subset of data for testing

**Port Already in Use**:
```r
shiny::runApp(port = 8080)  # Use different port
```

## 📝 Example Workflow

1. **Start the application**: `shiny::runApp()`

2. **Drug-Omic Analysis**:
   - Navigate to "Drugs-omics Pairs Analysis"
   - Select a drug (e.g., "Paclitaxel")
   - Select omic type (e.g., "mRNA Expression")
   - Select specific gene (e.g., "ABCB1")
   - Click "Analyze" to generate plots

3. **Batch Analysis**:
   - Navigate to "Batch Features Associations"
   - Select feature types (e.g., "Drug Sensitivity" vs "mRNA")
   - Select specific features
   - Adjust p-value threshold
   - Run analysis and view volcano plot

4. **Export Results**:
   - Use download buttons to save plots (PDF)
   - Export data tables (CSV)
   - Save R objects for further analysis

## 🤝 Contributing

Contributions are welcome! Please see our [Contributing Guidelines](../CONTRIBUTING.md) for details.

## 📄 License

This project is licensed under the MPL-2 License - see the [LICENSE](../LICENSE) file for details.

## 🔗 Related Projects

- [DROMA](https://github.com/mugpeng/DROMA) - Main DROMA project
- [DROMA.Set](https://github.com/mugpeng/DROMA_Set) - R package for data management
- [DROMA.R](https://github.com/mugpeng/DROMA_R) - R package for analysis functions
- [DROMA_MCP](https://github.com/mugpeng/DROMA_MCP) - MCP server for AI integration
- [Shiny](https://shiny.rstudio.com/) - R web application framework

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/mugpeng/DROMA_Web/issues)
- **Discussions**: [GitHub Discussions](https://github.com/mugpeng/DROMA_Web/discussions)
- **Email**: [Contact DROMA Team](mailto:yc47680@um.edu.mo)

## Citation

If you use DROMA Web in your research, please cite:

```
Li, S., Peng, Y., Chen, M. et al. Facilitating integrative and personalized oncology omics analysis with UCSCXenaShiny. Commun Biol 7, 1200 (2024). https://doi.org/10.1038/s42003-024-06891-2
```

---

**DROMA Web** - Interactive Cancer Pharmacogenomics Analysis 🧬💊🌐