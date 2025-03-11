# R Functions for Phage Display Analysis Using Immcantation

This repository provides a collection of **R functions** for analyzing **phage display** data using the **Immcantation framework**. These functions help process, analyze, and visualize antibody repertoire data from phage display libraries.

---

## Features

- **Data preprocessing**: Clean and filter phage display data for analysis.
- **Gene usage analysis**: Analyze V(D)J gene segment usage patterns.
- **Clonal diversity and abundance**: Quantify and visualize diversity using Immcantation tools.
- **Sequence similarity and clustering**: Perform t-SNE and network analysis to group similar antibody sequences.
- **Statistical analysis**: Generate summary statistics and tables for repertoire analysis.

---

## Requirements

To use these functions, ensure you have the following:

- **R** (version 4.0 or higher)
- Required R packages:

  ```r
  install.packages(c(
    "alakazam", 
    "shazam", 
    "dplyr", 
    "ggplot2", 
    "igraph", 
    "tidyr", 
    "Rtsne", 
    "stringdist", 
    "plotly", 
    "ggfortify", 
    "RColorBrewer", 
    "rgl", 
    "DT", 
    "factoextra", 
    "tidyverse", 
    "reshape2", 
    "htmlwidgets", 
    "svglite", 
    "plotrix", 
    "plot3D", 
    "MHTdiscrete"
  ))

## Functions Overview

### Data Preprocessing
- **`load_data`**: Load and preprocess phage display data from AIRR-compatible files.
- **`data_summary`**: Summarize sequence frequencies and calculate standard errors.

### Gene Usage Analysis
- **`generate_vdj_usage_plot`**: Visualize V(D)J gene usage patterns.
- **`generate_pca_plot`**: Perform PCA to analyze gene usage trends.

### Clonal Diversity and Abundance
- **`generate_abundance_plot`**: Generate abundance curves for antibody sequences.
- **`generate_diversity_plot`**: Calculate and visualize clonal diversity.

### Sequence Similarity and Clustering
- **`generate_tsne_plot`**: Perform t-SNE for clustering antibody sequences by similarity.
- **`generate_network_plot`**: Create network graphs to analyze relationships between antibody clones.

### Statistical Analysis
- **`generate_statistics_table`**: Generate summary statistics for gene usage.
- **`generate_combinatorial_statistics_table`**: Analyze combinatorial patterns in V(D)J gene usage.
