# R Workshop 20-05-2024 (AP Dr Wong Kah Keng) #
---
```r
### R Workshop 20-05-2024 (AP Dr Wong Kah Keng) ###
### Step 1: Set Directory, Install and Load the Necessary Libraries ###
setwd("C:/Users/Wong/Desktop/RStudio_Workshop_2024/Selected_dataset/Limma_v2") # Set your directory using forward slashes (not backslashes)

# Function to install and load packages. Installs the package from CRAN or Bioconductor based on the bioc argument.
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

# Ensure Bioconductor is installed
install_and_load("BiocManager")
# Set the update policy to 'some'
update_policy <- "some"

# For data manipulation
install_and_load("dplyr") # To process datasets in R e.g., filtering, selecting, and transforming data.


