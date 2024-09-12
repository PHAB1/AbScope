# Shiny Web Application for NGS B Cell Repertoire Analysis

This repository contains a web application developed in R using the **Shiny** framework. The application is specifically designed for the analysis of **Next-Generation Sequencing (NGS) B cell repertoire** data, offering tools to explore, visualize, and analyze immune repertoire datasets, with a focus on **gene usage PCA visualizations**.

## Features

- Interactive analysis of NGS B cell repertoire data
- **Principal Component Analysis (PCA) visualizations** for gene usage, allowing users to explore V(D)J gene segment patterns
- Tools for dimensionality reduction (e.g., t-SNE) and clustering
- Network analysis of B cell clones
- 3D visualizations and advanced plotting
- Summary statistics and data tables for repertoire analysis

## Requirements

To run this Shiny app locally, ensure you have the following:

- R (version 4.0 or higher)
- The following R packages:

  ```r
  install.packages(c(
    "shiny", 
    "alakazam", 
    "shazam", 
    "dplyr", 
    "scales", 
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
