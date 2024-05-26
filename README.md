# R Workshop 20-05-2024 (AP Dr Wong Kah Keng) #
---
## Step 1: Set Directory, Install and Load the Necessary Libraries
setwd("C:/Users/Wong/Desktop/RStudio_Workshop_2024/Selected_dataset/Limma_v2") # Set your directory using forward slashes (not backslashes)

### Function to install and load packages. Installs the package from CRAN or Bioconductor based on the bioc argument.
install_and_load <- function(package, bioc = FALSE, update_policy = "some") {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(paste("Installing package:", package))
    if (bioc) {
      BiocManager::install(package, update = update_policy)
    } else {
      install.packages(package)
    }
  } else {
    message(paste("Package already installed:", package))
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

#### Ensure Bioconductor is installed
install_and_load("BiocManager")
#### Set the update policy to 'some'
update_policy <- "some"

### For data manipulation
install_and_load("dplyr") # To process datasets in R e.g., filtering, selecting, and transforming data.

### For differential expression analysis
install_and_load("edgeR", bioc = TRUE, update_policy = update_policy) # Specifically for RNA-seq data
install_and_load("limma", bioc = TRUE, update_policy = update_policy) # For gene expression data including RNA-seq

### For plotting graphs
install_and_load("ggplot2") # To produce customisable plots using a grammar of graphics
install_and_load("ggrepel") # Prevents text labels from overlapping
install_and_load("plotly") # Creates interactive web-based graphs
install_and_load("pheatmap") # Creates heatmaps
install_and_load("reshape2") # To flexibly reshape data between wide and long formats
install_and_load("wesanderson") # Colour palettes

### For working with files
install_and_load("openxlsx") # Reading, writing, and editing Excel files
install_and_load("htmlwidgets") # Creating interactive web visualisations
install_and_load("jsonlite") # Reading and writing JSON data (JavaScript Object Notation; to store and export data)

### For Gene Set Enrichment Analysis (GSEA)
install_and_load("GSEABase", bioc = TRUE, update_policy = update_policy) # GSEA
install_and_load("fgsea", bioc = TRUE, update_policy = update_policy) # For fast preranked GSEA
install_and_load("clusterProfiler", bioc = TRUE, update_policy = update_policy) # For analysis and visualisation of functional profiles (gene clusters).
install_and_load("org.Hs.eg.db", bioc = TRUE, update_policy = update_policy) # Provides annotations (e.g., Entrez Gene IDs, gene symbols, chr locations. org = organism; eg = Entrez Gene; db: Indicates that the package is a database)
install_and_load("msigdbr") # Molecular Signatures Database (MSigDB) i.e., a collection of annotated gene sets for use with GSEA

### Additional required dependencies
install_and_load("fastmap") # Provides fast and memory-efficient key-value stores.
install_and_load("tibble") # An modern update of data frames in R.
install_and_load("tidyr") # Helps to tidy data by converting it into a more readable and analysable format.
install_and_load("magrittr") # Provides the pipe operator (`%>%`) for improving the readability and usability of codes.
install_and_load("stringr") # Simplifies string operations and manipulations in R.
install_and_load("purrr") # Enhances R's functional programming toolkit, making it easier to work with lists and other vectorised operations.


## Section 2: Retrieve and Load the RNA-seq Dataset to R, and Save the Session
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229136/suppl/GSE229136_G816_RNA-seq_counts.csv.gz" # URL leading to the file (GSE229136_G816_RNA-seq_counts.csv.gz) hosted on the NCBI GEO FTP server.

### Destination file path
destfile <- "GSE229136_G816_RNA-seq_counts.csv.gz" # Defines a variable destfile where the downloaded file will be saved. 

### Download the dataset
download.file(url, destfile, mode = "wb") # Downloads the file from the specified url and saves it to destfile. "wb" is to write data to a file in binary format i.e., data is written to the file exactly as it is. 

### Load the dataset
data <- read.csv(gzfile(destfile), header = TRUE, row.names = 1) # gzfile() function is used to open a connection to a Gzip-compressed file. destfile is the path to the Gzip-compressed CSV file. `header = TRUE` indicates that the first row of the CSV file contains column names. `row.names = 1` specifies that the first column of the CSV file should be used as row names for the data frame.
#### For many R scripts in the form of `xyz <- abc(def, ghi)` 
    # abc: The function being called.
    # def: The dataset or primary input to the function.
    # ghi: Options, parameters, or additional arguments/dependencies passed to the function.
    # xyz: The output or result of the function call, assigned to a variable.

# Alternatively, directly load it from own directory
# First, download the file from:  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229136/suppl/GSE229136_G816_RNA-seq_counts.csv.gz
# Then the following steps:
file_path <- "C:/Users/Wong/Desktop/RStudio_Workshop_2024/Selected_dataset/Limma_v2/GSE229136_G816_RNA-seq_counts.csv.gz" 
# Load the dataset
data <- read.csv(gzfile(file_path), header = TRUE, row.names = 1)

