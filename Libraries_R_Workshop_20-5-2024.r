# Define the required packages
required_packages <- c("dplyr", "edgeR", "limma", "ggplot2", "ggrepel", "plotly", "pheatmap", "reshape2", "wesanderson", "openxlsx", "htmlwidgets", "jsonlite", "GSEABase", "fgsea", "clusterProfiler", "org.Hs.eg.db", "msigdbr", "fastmap", "tibble", "tidyr", "magrittr", "stringr", "purrr")

# Load the libraries in one go
lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})
