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

# For differential expression analysis
install_and_load("edgeR", bioc = TRUE, update_policy = update_policy) # Specifically for RNA-seq data
install_and_load("limma", bioc = TRUE, update_policy = update_policy) # For gene expression data including RNA-seq

# For plotting graphs
install_and_load("ggplot2") # To produce customisable plots using a grammar of graphics
install_and_load("ggrepel") # Prevents text labels from overlapping
install_and_load("plotly") # Creates interactive web-based graphs
install_and_load("pheatmap") # Creates heatmaps
install_and_load("reshape2") # To flexibly reshape data between wide and long formats
install_and_load("wesanderson") # Colour palettes

# For working with files
install_and_load("openxlsx") # Reading, writing, and editing Excel files
install_and_load("htmlwidgets") # Creating interactive web visualisations
install_and_load("jsonlite") # Reading and writing JSON data (JavaScript Object Notation; to store and export data)

# For Gene Set Enrichment Analysis (GSEA)
install_and_load("GSEABase", bioc = TRUE, update_policy = update_policy) # GSEA
install_and_load("fgsea", bioc = TRUE, update_policy = update_policy) # For fast preranked GSEA
install_and_load("clusterProfiler", bioc = TRUE, update_policy = update_policy) # For analysis and visualisation of functional profiles (gene clusters).
install_and_load("org.Hs.eg.db", bioc = TRUE, update_policy = update_policy) # Provides annotations (e.g., Entrez Gene IDs, gene symbols, chr locations. org = organism; eg = Entrez Gene; db: Indicates that the package is a database)
install_and_load("msigdbr") # Molecular Signatures Database (MSigDB) i.e., a collection of annotated gene sets for use with GSEA

# Additional required dependencies
install_and_load("fastmap") # Provides fast and memory-efficient key-value stores.
install_and_load("tibble") # An modern update of data frames in R.
install_and_load("tidyr") # Helps to tidy data by converting it into a more readable and analysable format.
install_and_load("magrittr") # Provides the pipe operator (`%>%`) for improving the readability and usability of codes.
install_and_load("stringr") # Simplifies string operations and manipulations in R.
install_and_load("purrr") # Enhances R's functional programming toolkit, making it easier to work with lists and other vectorised operations.


### Section 2: Retrieve and Load the RNA-seq Dataset to R, and Save the Session ###
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229136/suppl/GSE229136_G816_RNA-seq_counts.csv.gz" # URL leading to the file (GSE229136_G816_RNA-seq_counts.csv.gz) hosted on the NCBI GEO FTP server.

# Destination file path
destfile <- "GSE229136_G816_RNA-seq_counts.csv.gz" # Defines a variable destfile where the downloaded file will be saved. 

# Download the dataset
download.file(url, destfile, mode = "wb") # Downloads the file from the specified url and saves it to destfile. "wb" is to write data to a file in binary format i.e., data is written to the file exactly as it is. 

# Load the dataset
data <- read.csv(gzfile(destfile), header = TRUE, row.names = 1) # gzfile() function is used to open a connection to a Gzip-compressed file. destfile is the path to the Gzip-compressed CSV file. `header = TRUE` indicates that the first row of the CSV file contains column names. `row.names = 1` specifies that the first column of the CSV file should be used as row names for the data frame.
# For many R scripts in the form of `xyz <- abc(def, ghi)` 
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

# Data view in RStudio

```
<img src="https://github.com/kahkengwong/R_Workshop_2024/blob/main/Images/Image_1_Data.jpg?raw=true" alt="View Data Summary Screenshot" width="600">
```r

# Check the first few rows of the dataset in R’s console
head(data)

# Basic data summary
summary(data)
```
<img src="https://github.com/kahkengwong/R_Workshop_2024/blob/main/Images/Image_2_Summary.jpg?raw=true" alt="View Data Summary Screenshot" width="400">

```r

# Capture the summary output and export as Excel file
summary_data <- as.data.frame(summary(data)) 
write.xlsx(summary_data, file = "Summary_output.xlsx")

# Check the number of genes before filtering
cat("Number of genes before filtering: ", nrow(data), "/n") # 23708 

# Save the R session
save.image(file = "R_Workshop_20-5-2024.RData")

# Optional: After closing R, reload the R Session and the required libraries using:
setwd("C:/…") # Set directory containing the saved R Session file (.RData)
load("R_Workshop_20-5-2024.RData") # Load the .RData file
source("Libraries_R_Workshop_20-5-2024.R") # Load all the required libraries from a manually-created .R file containing the instructions to launch and reload the libraries


### Section 3: RNA-seq Data Filtering, Normalization, and Voom Analysis ###
# Assign samples to groups
groups <- data.frame(
  samples = colnames(data),
  group = factor(c(rep("Unt", 3), rep("TT", 3), rep("TMZ", 3)), 
  levels = c("Unt", "TT", "TMZ"))
)

# Filter lowly expressed genes
keep.exprs <- filterByExpr(data, group = groups$group) 

# Subset the data
data.filt <- data[keep.exprs, ] 

# Normalise the data. Normalisation is required to adjust for differences in sequencing depth [total number of reads (sequenced fragments) obtained for a sample in RNA-seq] to make samples comparable and improve the accuracy of analyses.
dge <- DGEList(counts = data.filt) 
dge <- calcNormFactors(dge) 
logcpm <- cpm(dge, log = TRUE) 

# Create design matrix, specifying how each sample is assigned to the experimental groups
design <- model.matrix(~0 + groups$group) 
colnames(design) <- levels(groups$group) 

# Fit the linear model
fit <- lmFit(logcpm, design) 

# Define contrasts
contrasts <- makeContrasts(
  TT_vs_Unt = TT - Unt,
  TMZ_vs_Unt = TMZ - Unt,
  levels = design
) 

# Fit contrasts
fit <- contrasts.fit(fit, contrasts) 
fit <- eBayes(fit) 

# Voom transformation from the limma package
voom_data <- voom(dge, design) 


### Section 4: RNA-seq Quality Control (QC) Steps ###
# Section 4.1: QC Step with 2D PCA Plot (to identify and remove batch effects)
pca <- prcomp(t(voom_data$E)) 

pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = groups$group)

percentVar <- round(100 * summary(pca)$importance[2, 1:2])

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste("PC1: ", percentVar[1], "% variance")) +
    ylab(paste("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of RNA-seq data") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text_repel(aes(label = rownames(pca_data)))

```
![View Data Screenshot](https://github.com/kahkengwong/R_Workshop_2024/blob/main/Images/Image_3_2D-PCA.jpg?raw=true)
```r
 
# Section 4.2: QC Step with 3D PCA Plot (to identify and remove batch effects)
transposed_data <- t(voom_data$E)
pca <- prcomp(transposed_data, scale. = TRUE)
scores <- pca$x[, 1:3]
scores <- apply(scores, 2, function(x) x / max(abs(x)))

groups <- factor(c(rep("Unt", 3), rep("TT", 3), rep("TMZ", 3)))
colors <- c("blue", "green", "red")
group_colors <- colors[groups]

pca_data <- data.frame(PC1 = scores[, 1], PC2 = scores[, 2], PC3 = scores[, 3], 
                       Sample = colnames(voom_data$E), GroupColor = group_colors)

plot <- plot_ly(data = pca_data, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'markers',
                marker = list(size = 10, color = ~GroupColor), text = ~Sample) %>%
    layout(title = "3D PCA Plot",
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
plot

```
![View Data Screenshot](https://github.com/kahkengwong/R_Workshop_2024/blob/main/Images/Image_4_3D-PCA.jpg?raw=true)
```r

# Save the 3D plot in HTML format
saveWidget(plot, '3D_Plot_PCA.html', selfcontained = TRUE)


# Section 4.3: QC Step with Samples Correlation
sample_correlation <- cor(voom_data$E, method = "spearman")

pheatmap(sample_correlation, clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", 
         color = colorRampPalette(c("#3399FF", "white", "#FF3333"))(50), display_numbers = TRUE)

print(head(sample_correlation))

 
# Section 4.4: QC Step with Boxplot of Normalised Counts
# Section 4.4.1: With `wesanderson` Package
palette_wes <- wes_palette("FantasticFox1", n = ncol(voom_data$E), type = "continuous")

boxplot(voom_data$E, las = 2, col = palette_wes, main = "Boxplot of Normalised Counts (FantasticFox1)",
        ylab = "Log2 counts per million", xlab = "Samples",
        outpch = 19, outcol = "black", outcex = 0.5) 

 
# Section 4.4.2: With `ggplot2` Package (alternative style)
normalised_counts <- voom_data$E 
df <- as.data.frame(normalised_counts)
samples <- colnames(df)
df$Gene <- rownames(df)
df_long <- melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression") 

# Box plot of normalised counts
p <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "Log-normalised counts", title = "Expression Values by Sample")

print(p)

# Alternative box plot without outliers using 'outlier.shape = NA' argument
p_alt1 <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, fatten = 1.5, color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank()) +
    labs(x = "Sample", y = "Log-normalised counts", title = "Expression Values by Sample")

print(p_alt1)

 
# Section 4.5: QC Step with Mean-variance Plot
mean_counts <- rowMeans(data.filt)
var_counts <- apply(data.filt, 1, var)
df_before_voom <- data.frame(Mean = mean_counts, Variance = var_counts)
df_before_voom <- df_before_voom[mean_counts > 0 & var_counts > 0, ]

mean_voom <- rowMeans(normalised_counts)
var_voom <- apply(normalised_counts, 1, var)
df_after_voom <- data.frame(Mean = mean_voom, Variance = var_voom)
df_after_voom <- df_after_voom[mean_voom > 0 & var_voom > 0, ]

ggplot() +
    geom_point(data = df_before_voom, aes(x = Mean, y = Variance), color = "red", alpha = 0.4) +
    geom_point(data = df_after_voom, aes(x = Mean, y = Variance), color = "blue", alpha = 0.4) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Mean-Variance Relationship", x = "Mean Expression", y = "Variance") +
    theme_minimal() +
    ggtitle("Mean vs. Variance (Red: Before Voom, Blue: After Voom)")

 
# Section 4.6: Additional QC Steps (optional)
# Section 4.6.1: Library Size Visualisation
# The "library size" for each sample refers to the total number of reads (or fragments) that were sequenced for that sample. This is calculated by summing up all the raw read counts across all genes for each sample. This is to check for any significant differences in sequencing depth between samples.
library_sizes <- colSums(data.filt)

df_library_sizes <- data.frame(Sample = colnames(data.filt), LibrarySize = library_sizes)

ggplot(df_library_sizes, aes(x = Sample, y = LibrarySize)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Library Sizes Across Samples", x = "Sample", y = "Library Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Section 4.6.2: BH-adjusted p-values Visualisation
# Extract BH-adjusted p-values
fit2 <- eBayes(fit)

tt_vs_untreated_pvals <- topTable(fit2, coef = "TT_vs_Unt", number = Inf, sort.by = "P")

tmz_vs_untreated_pvals <- topTable(fit2, coef = "TMZ_vs_Unt", number = Inf, sort.by = "P")

# Plot BH-adjusted p-values for inspection (TT vs untreated)
hist(tt_vs_untreated_pvals$adj.P.Val, breaks = 50, main = "BH-Adjusted p-values (TT vs Untreated)", xlab = "BH-Adjusted p-value", ylab = "Frequency", col = "lightgreen")

# Plot BH-adjusted p-values for inspection (TMZ vs untreated)
hist(tmz_vs_untreated_pvals$adj.P.Val, breaks = 50, main = "BH-Adjusted p-values (TMZ vs Untreated)", xlab = "BH-Adjusted p-value", ylab = "Frequency", col = "lightcoral")

# Section 4.6.3: Log2 Fold Changes Visualisation
# Extract log2 fold changes
tt_vs_untreated_lfc <- tt_vs_untreated_pvals$logFC
tmz_vs_untreated_lfc <- tmz_vs_untreated_pvals$logFC

# Plot distribution of log2 fold changes (TMZ vs untreated)
hist(tmz_vs_untreated_lfc, breaks = 50, main = "Log2 Fold Changes (TMZ vs Untreated)", xlab = "Log2 Fold Change", ylab = "Frequency", col = "lightcoral")

# Plot distribution of log2 fold changes (TT vs untreated)
hist(tt_vs_untreated_lfc, breaks = 50, main = "Log2 Fold Changes (TT vs Untreated)", xlab = "Log2 Fold Change", ylab = "Frequency", col = "lightgreen")

# Section 4.6.4: Gene Variances Visualisation
gene_variances <- apply(voom_data$E, 1, var)

hist(gene_variances, breaks = 50, main = "Distribution of Gene Variances", xlab = "Variance", ylab = "Frequency", col = "lightblue")
abline(h = 0, col = "blue")

# Save the updated R session (after all the QC steps)
save.image(file = "R_Workshop_20-5-2024.RData")


### Section 5: Differentially Expressed Genes (DEGs) Analysis ###
# Obtain List of DEGs for TT vs Untreated
tt_vs_unt <- topTable(fit, coef = "TT_vs_Unt", number = Inf, adjust.method = "BH")
tt_vs_unt$Gene <- rownames(tt_vs_unt)

# Obtain List of DEGs for TMZ vs Untreated
tmz_vs_unt <- topTable(fit, coef = "TMZ_vs_Unt", number = Inf, adjust.method = "BH")
tmz_vs_unt$Gene <- rownames(tmz_vs_unt)

# Reorder columns to put Gene as the first column
tt_vs_unt <- tt_vs_unt[, c("Gene", setdiff(colnames(tt_vs_unt), "Gene"))]
tmz_vs_unt <- tmz_vs_unt[, c("Gene", setdiff(colnames(tmz_vs_unt), "Gene"))]

# Save all genes to Excel
write.xlsx(list(TT_vs_Untreated = tt_vs_unt, TMZ_vs_Untreated = tmz_vs_unt), file = "All_Genes_TT_TMZ_vs_Untreated_v2.xlsx")
# In the Excel
1)	AveExpr (Average Expression)
2)	t (t-statistic): The difference in expression between the specified groups (e.g., TT vs. Untreated or TMZ vs. Untreated) relative to the variability in expression within groups
3)	B (log-odds): The log-odds that a gene is differentially expressed. A higher B value = greater confidence that the gene is differentially expressed. 


# Section 6: Volcano Plots of the Hallmark Gene Sets 
# Section 6.1: HALLMARK_TGF_BETA_SIGNALING 
# URL for the HALLMARK_TGF_BETA_SIGNALING gene set JSON file
json_url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_TGF_BETA_SIGNALING&fileType=json"

# Download the JSON file and read it directly into R
json_data <- fromJSON(json_url)

# Extract the list of gene symbols
hallmark_genes <- json_data$HALLMARK_TGF_BETA_SIGNALING$geneSymbols

# Print the first few genes to verify
print(head(hallmark_genes))

# Check the number of genes
length(hallmark_genes)

# Filter the limma results to include only genes in the hallmark_genes list
hallmark_results_TT_vs_Untreated <- tt_vs_unt[tt_vs_unt$Gene %in% hallmark_genes, ]
hallmark_results_TMZ_vs_Untreated <- tmz_vs_unt[tmz_vs_unt$Gene %in% hallmark_genes, ]

# Define colors
colors <- c("Not significant" = "grey", "Significant" = "navyblue")


### Minimal theme volcano plots
# Create volcano plots for TT vs Untreated
volcano_TT <- ggplot(hallmark_results_TT_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val), 
                                                           color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant"))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TT vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = c("Not significant" = "grey", "Significant" = "navyblue")) +
    theme(legend.position = "right")

# Display the plot
print(volcano_TT)


# Create volcano plots for TMZ vs Untreated
volcano_TMZ <- ggplot(hallmark_results_TMZ_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val), 
                                                             color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant"))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TMZ vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = c("Not significant" = "grey", "Significant" = "navyblue")) +
    theme(legend.position = "right")

# Display the plot
print(volcano_TMZ)


### Minimal theme but with gene labels
# Create volcano plots for TT vs Untreated
volcano_TT <- ggplot(hallmark_results_TT_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val), 
                                                           color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant"))) +
    geom_point(alpha = 0.5) +
    geom_text_repel(data = subset(hallmark_results_TT_vs_Untreated, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TT vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = c("Not significant" = "grey", "Significant" = "navyblue")) +
    theme(legend.position = "right")

# Print the plot
print(volcano_TT)

# Create volcano plots for TMZ vs Untreated
volcano_TMZ <- ggplot(hallmark_results_TMZ_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val), 
                                                             color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant"))) +
    geom_point(alpha = 0.5) +
    geom_text_repel(data = subset(hallmark_results_TMZ_vs_Untreated, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TMZ vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = c("Not significant" = "grey", "Significant" = "navyblue")) +
    theme(legend.position = "right")

# Print the plot
print(volcano_TMZ)


### With gene labels, border and fixed legends
# Create volcano plots for TT vs Untreated
volcano_TT <- ggplot(hallmark_results_TT_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")), alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_text_repel(data = subset(hallmark_results_TT_vs_Untreated, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TT vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = colors) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
# The aes(color = ...) inside the main ggplot call affected the overall aesthetics, leading to the legend issue. Fixed by moving the aes(color = ...) mapping to within the geom_point function.

# Display the plot
print(volcano_TT)
 
# Create volcano plots for TMZ vs Untreated
volcano_TMZ <- ggplot(hallmark_results_TMZ_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")), alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_text_repel(data = subset(hallmark_results_TMZ_vs_Untreated, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TGF_BETA_SIGNALING: TMZ vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = colors) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
# Display the plot
print(volcano_TMZ)


# Section 6.2: HALLMARK_TNFA_SIGNALING_VIA_NFKB
# URL for the HALLMARK_TNFA_SIGNALING_VIA_NFKB gene set JSON file
json_url_TNFa <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_TNFA_SIGNALING_VIA_NFKB&fileType=json"

# Download the JSON file and read it directly into R
json_data_TNFa <- fromJSON(json_url_TNFa)

# Extract the list of gene symbols
hallmark_genes_TNFa <- json_data_TNFa$HALLMARK_TNFA_SIGNALING_VIA_NFKB$geneSymbols

# Print the first few genes to verify
print(head(hallmark_genes_TNFa))

# Check the number of genes
length(hallmark_genes_TNFa)

# Filter the limma results to include only genes in the hallmark_genes_TNFa list
hallmark_results_TT_vs_Untreated_TNFa <- tt_vs_unt[tt_vs_unt$Gene %in% hallmark_genes_TNFa, ]
hallmark_results_TMZ_vs_Untreated_TNFa <- tmz_vs_unt[tmz_vs_unt$Gene %in% hallmark_genes_TNFa, ]

# Define colors
colors <- c("Not significant" = "grey", "Significant" = "navyblue")

# Create volcano plots for TT vs Untreated
volcano_TT_TNFa <- ggplot(hallmark_results_TT_vs_Untreated_TNFa, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")), alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_text_repel(data = subset(hallmark_results_TT_vs_Untreated_TNFa, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB: TT vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = colors) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# Display the plot
print(volcano_TT_TNFa)

# Create volcano plots for TT vs Untreated - Remove short lines that connect a gene name to its dot. The ggrepel package, specifically the geom_text_repel() function in ggplot2, includes those short lines by default to improve the clarity of the plot. 
volcano_TT_TNFa <- ggplot(hallmark_results_TT_vs_Untreated_TNFa, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")), alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_text_repel(data = subset(hallmark_results_TT_vs_Untreated_TNFa, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic", 
                    segment.color = NA) +
    theme_minimal() +
    labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB: TT vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = colors) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# Display the plot
print(volcano_TT_TNFa)


# Create volcano plots for TMZ vs Untreated
volcano_TMZ_TNFa <- ggplot(hallmark_results_TMZ_vs_Untreated_TNFa, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")), alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_text_repel(data = subset(hallmark_results_TMZ_vs_Untreated_TNFa, abs(logFC) > 1 & adj.P.Val < 0.05),
                    aes(label = Gene), size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, fontface = "italic") +
    theme_minimal() +
    labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB: TMZ vs Untreated",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value (BH)") +
    scale_color_manual(name = "Significance", values = colors) +
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# Display the plot
print(volcano_TMZ_TNFa)


### Section 7: Gene Set Enrichment Analysis (GSEA) ###
# Section 7.1: HALLMARK_TGF_BETA_SIGNALING
# Load the hallmark gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_tgf_beta <- hallmark_gene_sets[hallmark_gene_sets$gs_name == "HALLMARK_TGF_BETA_SIGNALING", ]$gene_symbol

# Create rank list for TT vs Untreated
tt_vs_unt_rank <- tt_vs_unt$logFC
names(tt_vs_unt_rank) <- tt_vs_unt$Gene
tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)

# Create rank list for TMZ vs Untreated
tmz_vs_unt_rank <- tmz_vs_unt$logFC
names(tmz_vs_unt_rank) <- tmz_vs_unt$Gene
tmz_vs_unt_rank <- sort(tmz_vs_unt_rank, decreasing = TRUE)

# Set seed and ensure proper ranking
set.seed(1) 
tt_vs_unt_rank <- tt_vs_unt$logFC + rnorm(length(tt_vs_unt$logFC), sd = 1e-5)
names(tt_vs_unt_rank) <- tt_vs_unt$Gene
tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)

tmz_vs_unt_rank <- tmz_vs_unt$logFC + rnorm(length(tmz_vs_unt$logFC), sd = 1e-5)
names(tmz_vs_unt_rank) <- tmz_vs_unt$Gene
tmz_vs_unt_rank <- sort(tmz_vs_unt_rank, decreasing = TRUE)

# Run GSEA for TT vs Untreated
fgsea_res_tt <- fgseaMultilevel(pathways = list(HALLMARK_TGF_BETA_SIGNALING = hallmark_tgf_beta), 
                                stats = tt_vs_unt_rank)

# Run GSEA for TMZ vs Untreated
fgsea_res_tmz <- fgseaMultilevel(pathways = list(HALLMARK_TGF_BETA_SIGNALING = hallmark_tgf_beta), 
                                 stats = tmz_vs_unt_rank)

# Function to plot GSEA results
plotGseaTable <- function(fgseaRes, stats, title) {
    # Ensure fgseaRes is not empty
    if (nrow(fgseaRes) > 0) {fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
        pathway <- hallmark_tgf_beta
        enrichmentScore <- fgseaRes$ES[1]
        nes <- fgseaRes$NES[1]
        pval <- fgseaRes$padj[1]
        
        gseaPlot <- plotEnrichment(pathway = pathway,stats = stats) + ggtitle(title) + theme_minimal() + labs(y = "Running Enrichment Score", x = "Rank in Ordered Gene List") + annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "/nFDR =", round(pval, 4)), hjust = 1.1, vjust = 1.1, size = 5, color = "blue", fontface = "bold")
        
        print(gseaPlot)
    } else {message("No significant enrichment found for ", title)}
}

# Plot GSEA results for TT vs Untreated
plotGseaTable(fgsea_res_tt, tt_vs_unt_rank, "TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable(fgsea_res_tmz, tmz_vs_unt_rank, "TMZ vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

## Use darker green line and FDR value in scientific notation using formatC function ##
# Function to plot GSEA results 
plotGseaTable <- function(fgseaRes, stats, title) {
    # Ensure fgseaRes is not empty
    if (nrow(fgseaRes) > 0) {
        fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
        pathway <- hallmark_tgf_beta  # Use the full gene set for the pathway
        enrichmentScore <- fgseaRes$ES[1]
        nes <- fgseaRes$NES[1]
        pval <- fgseaRes$padj[1]
        
        # Format the FDR value in scientific notation
        formatted_pval <- formatC(pval, format = "e", digits = 2) # The digits = 2 parameter ensures that the FDR value is displayed with two significant digits, but can adjust this as needed.
        
        # Extract plot object and modify the running score line color
        gseaPlot <- plotEnrichment(pathway = pathway, stats = stats) + 
                    ggtitle(title) + 
                    theme_minimal() + 
                    labs(y = "Running Enrichment Score", x = "Rank in Ordered Gene List") + 
                    annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "/nFDR =", formatted_pval), 
                             hjust = 1.1, vjust = 1.1, size = 5, color = "blue", fontface = "bold")
        
        # Modify the line color
        gseaPlot <- gseaPlot +
                    theme(legend.position = "none") +
                    geom_line(data = ggplot_build(gseaPlot)$data[[1]], aes(x = x, y = y), color = "darkgreen", size = 1)
        
        print(gseaPlot)
    } else {
        message("No significant enrichment found for ", title)
    }
}

# Plot GSEA results for TT vs Untreated
plotGseaTable(fgsea_res_tt, tt_vs_unt_rank, "TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable(fgsea_res_tmz, tmz_vs_unt_rank, "TMZ vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

## Add thin black box line surrounding the GSEA plot and use navyblue for the fonts ##
#Adjust the annotate function to modify the font size and color. Add a theme element to draw the box around the plot.
# Function to plot GSEA results
plotGseaTable <- function(fgseaRes, stats, title) {
    # Ensure fgseaRes is not empty
    if (nrow(fgseaRes) > 0) {
        fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
        pathway <- hallmark_tgf_beta  # Use the full gene set for the pathway
        enrichmentScore <- fgseaRes$ES[1]
        nes <- fgseaRes$NES[1]
        pval <- fgseaRes$padj[1]
        
        # Format the FDR value in scientific notation
        formatted_pval <- formatC(pval, format = "e", digits = 2)
        
        # Extract plot object and modify the running score line color
        gseaPlot <- plotEnrichment(pathway = pathway, stats = stats) + 
                    ggtitle(title) + 
                    theme_minimal() + 
                    labs(y = "Running Enrichment Score", x = "Rank in Ordered Gene List") + 
                    annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "\nFDR =", formatted_pval), 
                             hjust = 1.1, vjust = 1.1, size = 4, color = "navy", fontface = "bold") +
                    theme(
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                        axis.line = element_line(color = "black", linewidth = 0.5)
                    )
        
        # Modify the line color
        gseaPlot <- gseaPlot +
                    theme(legend.position = "none") +
                    geom_line(data = ggplot_build(gseaPlot)$data[[1]], aes(x = x, y = y), color = "darkgreen", linewidth = 1)
        
        print(gseaPlot)
    } else {
        message("No significant enrichment found for ", title)
    }
}

# Plot GSEA results for TT vs Untreated
plotGseaTable(fgsea_res_tt, tt_vs_unt_rank, "TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable(fgsea_res_tmz, tmz_vs_unt_rank, "TMZ vs Untreated: HALLMARK_TGF_BETA_SIGNALING")


# Section 7.2: HALLMARK_TNFA_SIGNALING_VIA_NFKB
# Load the hallmark gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_tnf_alpha <- hallmark_gene_sets[hallmark_gene_sets$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB", ]$gene_symbol

# Create rank list for TT vs Untreated
tt_vs_unt_rank_TNFa <- tt_vs_unt$logFC
names(tt_vs_unt_rank_TNFa) <- tt_vs_unt$Gene
tt_vs_unt_rank_TNFa <- sort(tt_vs_unt_rank_TNFa, decreasing = TRUE)

# Create rank list for TMZ vs Untreated
tmz_vs_unt_rank_TNFa <- tmz_vs_unt$logFC
names(tmz_vs_unt_rank_TNFa) <- tmz_vs_unt$Gene
tmz_vs_unt_rank_TNFa <- sort(tmz_vs_unt_rank_TNFa, decreasing = TRUE)

# Set seed and ensure proper ranking
set.seed(1)
tt_vs_unt_rank_TNFa <- tt_vs_unt$logFC + rnorm(length(tt_vs_unt$logFC), sd = 1e-5)
names(tt_vs_unt_rank_TNFa) <- tt_vs_unt$Gene
tt_vs_unt_rank_TNFa <- sort(tt_vs_unt_rank_TNFa, decreasing = TRUE)

tmz_vs_unt_rank_TNFa <- tmz_vs_unt$logFC + rnorm(length(tmz_vs_unt$logFC), sd = 1e-5)
names(tmz_vs_unt_rank_TNFa) <- tmz_vs_unt$Gene
tmz_vs_unt_rank_TNFa <- sort(tmz_vs_unt_rank_TNFa, decreasing = TRUE)

# Run GSEA for TT vs Untreated
fgsea_res_tt_TNFa <- fgseaMultilevel(pathways = list(HALLMARK_TNFA_SIGNALING_VIA_NFKB = hallmark_tnf_alpha), 
                                     stats = tt_vs_unt_rank_TNFa)

# Run GSEA for TMZ vs Untreated
fgsea_res_tmz_TNFa <- fgseaMultilevel(pathways = list(HALLMARK_TNFA_SIGNALING_VIA_NFKB = hallmark_tnf_alpha), 
                                      stats = tmz_vs_unt_rank_TNFa)

# Function to plot GSEA results
plotGseaTable_TNFa <- function(fgseaRes, stats, title) {
    # Ensure fgseaRes is not empty
    if (nrow(fgseaRes) > 0) {
        fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
        pathway <- hallmark_tnf_alpha  # Use the full gene set for the pathway
        enrichmentScore <- fgseaRes$ES[1]
        nes <- fgseaRes$NES[1]
        pval <- fgseaRes$padj[1]
        
        # Format the FDR value in scientific notation
        formatted_pval <- formatC(pval, format = "e", digits = 2)
        
        # Extract plot object and modify the running score line color
        gseaPlot <- plotEnrichment(pathway = pathway, stats = stats) + 
                    ggtitle(title) + 
                    theme_minimal() + 
                    labs(y = "Running Enrichment Score", x = "Rank in Ordered Gene List") + 
                    annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "/nFDR =", formatted_pval), 
                             hjust = 1.1, vjust = 1.1, size = 4, color = "navy", fontface = "bold") +
                    theme(
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                        axis.line = element_line(color = "black", linewidth = 0.5)
                    )
        
        # Modify the line color
        gseaPlot <- gseaPlot +
                    theme(legend.position = "none") +
                    geom_line(data = ggplot_build(gseaPlot)$data[[1]], aes(x = x, y = y), color = "darkgreen", linewidth = 1)
        
        print(gseaPlot)
    } else {
        message("No significant enrichment found for ", title)
    }
}

# Plot GSEA results for TT vs Untreated
plotGseaTable_TNFa(fgsea_res_tt_TNFa, tt_vs_unt_rank_TNFa, "TT vs Untreated: HALLMARK_TNFA_SIGNALING_VIA_NFKB")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable_TNFa(fgsea_res_tmz_TNFa, tmz_vs_unt_rank_TNFa, "TMZ vs Untreated: HALLMARK_TNFA_SIGNALING_VIA_NFKB")

# Save the R Session 
save.image(file = "R_Workshop_20-5-2024.RData")

# Reload the R Session and the Required Libraries
setwd("C:/Users/Wong/Desktop/RStudio_Workshop_2024/Selected_dataset/Limma_v2")
load("R_Workshop_20-5-2024.RData")
source("Libraries_R_Workshop_20-5-2024.R") 
