# EVOmicsDB

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**EVOmicsDB: a clinically annotated computational framework for pan-cancer extracellular vesicle multi-omics analysis**

**Versioned archive:** [10.5281/zenodo.19534750](https://doi.org/10.5281/zenodo.19534750)

## What This Repository Contains

This repository contains the core analytical workflows, reproducibility scripts, and dataset metadata manifest associated with [EVOmicsDB](https://EVOmicsDB.com). It does not include the production web application source code. It is intended to support reproducibility and allow users to run the core analyses described in the accompanying manuscript locally.

**Included:**
- R scripts for differential analysis, enrichment, biomarker evaluation, visualization, network construction, and cross-omics integration
- Dataset metadata manifest (`metadata/EVOmicsDB_manifest.json`) describing all 102 curated datasets
- Reproducibility scripts that map directly to manuscript figures
- Minimal example data for verifying script execution

**Not included:**
- The EVOmicsDB web application source code (the production web platform is deployed and maintained separately at [https://EVOmicsDB.com](https://EVOmicsDB.com))
- Full curated expression matrices (available for download from the web portal and/or original public repositories)

## Repository Structure

```text
EVOmicsDB/
├── R/
│   ├── analysis/                        # Core statistical analysis
│   │   ├── differential_analysis.R      # limma-based differential expression
│   │   ├── enrichment_analysis.R        # GO/KEGG over-representation analysis
│   │   ├── gsea.R                       # Gene Set Enrichment Analysis
│   │   └── pca.R                        # Principal Component Analysis
│   ├── biomarker/                       # Biomarker evaluation
│   │   ├── roc.R                        # Single-feature ROC curves
│   │   ├── roc_logistic.R              # Multi-marker logistic regression panel
│   │   ├── roc_rf.R                     # Random Forest feature selection
│   │   └── roc_svm.R                    # SVM feature selection
│   ├── visualization/                   # Plotting
│   │   ├── volcano_plot.R
│   │   ├── heatmap.R
│   │   ├── go_bar_plot.R / go_bubble.R
│   │   └── kegg_bar_plot.R / kegg_bubble_plot.R
│   ├── network/                         # Network construction
│   │   ├── ppi.R                        # Protein–protein interaction (STRINGdb)
│   │   ├── regulatory_network.R         # miRNA→mRNA/Protein regulatory network
│   │   ├── lncrna_target_network.R      # lncRNA target intersection network
│   │   ├── go_network.R                 # Integrated GO network (mRNA + Protein)
│   │   └── kegg_network.R              # KEGG joint pathway network
│   └── cross_omics/
│       └── cross_omics_trend_analysis.R # Cross-layer correlation analysis
├── metadata/
│   ├── EVOmicsDB_manifest.json          # Metadata for all 102 curated datasets
│   └── metadata_fields.md              # Field definitions
├── example_data/
│   ├── demo_expression_matrix.tsv       # Minimal example data (expression matrix)
│   ├── demo_sample_metadata.tsv         # Minimal example data (group labels)
│   └── README_example_data.md
├── reproducibility/
│   ├── reproduce_figure4_crc_case.R     # CRC EV proteomics case study
│   ├── reproduce_figure5_hsp90aa1.R     # Pan-cancer HSP90AA1 evaluation
│   ├── reproduce_figure6_crossomics.R   # Cross-omics integration
│   └── README_reproducibility.md
├── docs/
│   ├── USAGE.md                         # Detailed usage & input/output reference
│   └── availability_and_implementation.md
├── sessionInfo.txt                      # R session snapshot for reproducibility
├── CITATION.cff
├── LICENSE
├── .gitignore
└── README.md
```

## System Requirements

- **R** ≥ 4.2.0 (developed and tested under R 4.3.x)
- **OS:** Linux, macOS, or Windows

## Installation

### Required R packages

**CRAN:**

```r
install.packages(c(
  "tidyverse", "argparse", "ggplot2", "ggpubr", "pheatmap",
  "pROC", "caret", "e1071", "randomForest",
  "igraph", "ggraph", "visNetwork", "jsonlite",
  "FactoMineR", "factoextra", "scales", "stringr", "dplyr"
))
```

**Bioconductor:**

```r
install.packages("BiocManager")

BiocManager::install(c(
  "limma", "edgeR", "impute", "imputeLCMD",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db",
  "multiMiR", "STRINGdb", "KEGGREST"
))
```

A complete session snapshot is provided in [`sessionInfo.txt`](sessionInfo.txt).

## Quick Start with Example Data

```bash
# 0. Create the output directory (not tracked by Git)
mkdir -p results

# 1. Run differential analysis on the demo matrix
Rscript R/analysis/differential_analysis.R \
  --counts_file example_data/demo_expression_matrix.tsv \
  --output_rda results/demo_DEG.rda \
  --output_csv results/demo_filtered_DEG.csv

# 2. Generate a volcano plot
Rscript R/visualization/volcano_plot.R \
  --input_rda results/demo_DEG.rda \
  --output_png results/demo_volcano.png

# 3. Single-feature ROC analysis
Rscript R/biomarker/roc.R \
  --input_rda results/demo_DEG.rda \
  --output_png results/demo_roc.png \
  --output_csv results/demo_roc.csv \
  --top_n 3
```

All scripts accept `--help` for full parameter documentation.


## Reproducing Manuscript Figures

| Script | Manuscript Figure | Dataset | Description |
|--------|-------------------|---------|-------------|
| `reproduce_figure4_crc_case.R` | Figure 4 | PXD047810 | CRC EV proteomics: DE → enrichment → PPI → ROC |
| `reproduce_figure5_hsp90aa1.R` | Figure 5 | Multiple (16 datasets) | Pan-cancer HSP90AA1 ROC evaluation |
| `reproduce_figure6_crossomics.R` | Figure 6 | GSE255775, GSE220445, PXD047810 | Cross-omics trend, regulatory network, pathway network |

**Prerequisites:** Download the relevant datasets from [EVOmicsDB](https://EVOmicsDB.com) or their original repositories (accession numbers listed in [`metadata/EVOmicsDB_manifest.json`](metadata/EVOmicsDB_manifest.json)).

## Data and Code Scope

This repository contains analytical code and metadata only. It does not include the production web application source code. Full curated data matrices are available through the [EVOmicsDB](https://EVOmicsDB.com) web portal and/or the original public repositories, subject to the corresponding source data policies. See [`docs/availability_and_implementation.md`](docs/availability_and_implementation.md) for details.

## Citation

If you use EVOmicsDB or these analytical workflows, please cite:

> Kang B†, Peng P†, Yang Q, Geng X, Zhang L, Bai L, Zhu Y, Ju Y\*, Min L\*. EVOmicsDB: a clinically annotated computational framework for pan-cancer extracellular vesicle multi-omics analysis. Manuscript submitted to *Bioinformatics*.

A machine-readable citation file is provided in [`CITATION.cff`](CITATION.cff).

## License

Code in this repository is released under the [MIT License](LICENSE). Dataset metadata in `metadata/` is provided for reference and traceability; original data are subject to the policies of their respective source repositories.

## Contact

- **Li Min** — minli@ccmu.edu.cn
- **Yingjiao Ju** — juyingjiao@outlook.com

Beijing Friendship Hospital, Capital Medical University, Beijing, China
