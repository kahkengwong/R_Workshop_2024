setwd("C:/Users/…") # Set your directory using forward slashes (not backslashes)
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
# `install_and_load`: Defines a function to install and load R packages.
# `package`: The name of the package to install and load.
# `bioc = FALSE`: Indicates whether the package should be installed from Bioconductor (default is FALSE).
# `update_policy = "some"`: Specifies the update policy for Bioconductor packages (default is "some").
    # `if (!requireNamespace(package, quietly = TRUE))`: Checks if the package is not already installed, quietly suppressing messages.
    # `message(paste("Installing package:", package))`: Prints a message indicating that the package is being installed.
      # `BiocManager::install(package, update = update_policy)`: Installs the package from Bioconductor with the specified update policy.
      # `install.packages(package)`: Installs the package from CRAN.
    # `message(paste("Package already installed:", package))`: Prints a message indicating that the package is already installed.
  # `suppressPackageStartupMessages(library(package, character.only = TRUE))`: Loads the package, suppressing startup messages.

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
install_and_load("ggplot2") # To produce customisable plots using a grammar of graphics (i.e., a layered approach, enabling us to precisely describe each component of a graphic)
install_and_load("ggrepel") # Prevents text labels from overlapping
install_and_load("plotly") # Creates interactive web-based graphs
install_and_load("pheatmap") # Creates heatmaps
install_and_load("reshape2") # To flexibly reshape data between wide and long formats
install_and_load("wesanderson") # Colour palettes

# For working with files
install_and_load("openxlsx") # Reading, writing, and editing Excel files
install_and_load("htmlwidgets") # Creating interactive web visualisations
install_and_load("jsonlite") # Reading and writing JSON data (JavaScript Object Notation; to store and export data)

# For gene set enrichment analysis
install_and_load("GSEABase", bioc = TRUE, update_policy = update_policy) # Gene Set Enrichment Analysis (GSEA)
install_and_load("fgsea", bioc = TRUE, update_policy = update_policy) # For fast preranked GSEA
install_and_load("clusterProfiler", bioc = TRUE, update_policy = update_policy) # For analysis and visualisation of functional profiles (gene clusters).
install_and_load("org.Hs.eg.db", bioc = TRUE, update_policy = update_policy) # Provides genome wide annotation using Entrez Gene identifiers
install_and_load("msigdbr") # Molecular Signatures Database (MSigDB) i.e., a collection of annotated gene sets for use with GSEA

# Additional required dependencies
install_and_load("fastmap") # Provides fast and memory-efficient key-value stores.
install_and_load("tibble") # An modern update of data frames in R.
install_and_load("tidyr") # Helps to tidy data by converting it into a more readable and analysable format.
install_and_load("magrittr") # Provides the pipe operator (`%>%`) for improving the readability and usability of codes.
install_and_load("stringr") # Simplifies string operations and manipulations in R.
install_and_load("purrr") # Enhances R's functional programming toolkit, making it easier to work with lists and other vectorised operations.

# URL of the dataset (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229136)
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229136/suppl/GSE229136_G816_RNA-seq_counts.csv.gz" # URL leading to the file (GSE229136_G816_RNA-seq_counts.csv.gz) hosted on the NCBI GEO FTP server.

# Destination file path
destfile <- "GSE229136_G816_RNA-seq_counts.csv.gz" # Defines a variable destfile where the downloaded file will be saved. 

# Download the dataset
download.file(url, destfile, mode = "wb") # Downloads the file from the specified url and saves it to destfile. "wb" is to write data to a file in binary format i.e., data is written to the file exactly as it is. This means that no translations or conversions are applied to the data, ensuring that the file's content remains unchanged. Called "binary" because it treats the file as a sequence of bytes, without interpreting the content. It is suitable for any type of file, especially those that are not plain text (e.g., images, executables, compressed files).

# Load the dataset
data <- read.csv(gzfile(destfile), header = TRUE, row.names = 1) # gzfile() function is used to open a connection to a Gzip-compressed file. destfile is the path to the Gzip-compressed CSV file. `header = TRUE` indicates that the first row of the CSV file contains column names. `row.names = 1` specifies that the first column of the CSV file should be used as row names for the data frame.
# For many R scripts in the form of `xyz <- abc(def, ghi)` 
    # abc: The function being called.
    # def: The dataset or primary input to the function.
    # ghi: Options, parameters, or additional arguments/dependencies passed to the function.
    # xyz: The output or result of the function call, assigned to a variable.

# Check the first few rows of the dataset
head(data)
# Basic data summary
summary(data)

# Capture the summary output and export as Excel file
summary_data <- as.data.frame(summary(data)) # Captures the summary output of the data frame `data` and converts it into another data frame `summary_data`. The summary() function generates the summary statistics, and as.data.frame() converts this output into a data frame format, making it easier to manipulate and export. 

write.xlsx(summary_data, file = "Summary_Output.xlsx") # Writes the `summary_data` data frame to an Excel file

# Check the number of genes before filtering
cat("Number of genes before filtering: ", nrow(data), "\n") # 23708 

# Save and reload the R session
save.image(file = "R_Workshop_20-5-2024.RData")
setwd("C:/Users/Wong/Desktop/RStudio_Workshop_2024/Selected_dataset/Limma_v2")
load("R_Workshop_20-5-2024.RData")
source("Libraries_R_Workshop_20-5-2024.R") # The lapply function uses this vector to load each package, but the vector itself remains in the environment as a result of the script execution.

# Assign samples to groups
groups <- data.frame(
  samples = colnames(data),
  group = factor(c(rep("Unt", 3), rep("TT", 3), rep("TMZ", 3)), 
  levels = c("Unt", "TT", "TMZ"))
)

# `groups <- data.frame(`: Initiates the creation of a data frame named groups.
# `samples = colnames(data)`: Extracts the column names of the data frame data. These column names are assigned to the `samples` column in the new groups data frame.
# `group = factor(c(rep("Unt", 3), rep("TT", 3), rep("TMZ", 3))`: Creates a vector of group labels for the samples. The rep function repeats each of the group labels three times. The vector c("Unt", "Unt", "Unt", "TT", "TT", "TT", "TMZ", "TMZ", "TMZ") is created, corresponding to the group labels for the samples.

# Filter lowly expressed genes
keep.exprs <- filterByExpr(data, group = groups$group) # `data`: The expression data matrix; `group = groups$group`: The grouping factor indicating the experimental groups for each sample; `filterByExpr`: From edgeR package, and to identify genes that are sufficiently expressed across the samples, filtering out lowly expressed genes. Returns a logical vector (keep.exprs) indicating which genes pass the expression filter.

# Subset the data
data.filt <- data[keep.exprs, ] # Subsets the original data matrix data to include only the genes that passed the expression filter i.e., retaining only the sufficiently expressed genes. In R, when subsetting a matrix or data frame, we can specify both rows and columns within the square brackets []. The `syntax data[rows, columns]` is used to select specific rows and columns. The space after the comma , means "select all columns." Summary: keep.exprs: Specifies which rows (genes) to keep based on the filter. Space After Comma: Indicates that all columns (samples) should be retained.

# Normalise the data. Normalisation is required to adjust for differences in sequencing depth [total number of reads (sequenced fragments) obtained for a sample in RNA-seq] to make samples comparable and improve the accuracy of analyses.
dge <- DGEList(counts = data.filt) # `DGEList`: Creates a DGEList object which is a container used in edgeR to store count data for differential expression analysis. `counts = data.filt` (argument): The filtered count data matrix. Hence, `dge` is now a `DGEList` object containing the count data.
dge <- calcNormFactors(dge) # `calcNormFactors`: Normalises data for differences in sequencing depth across samples. `dge`: The DGEList object. Result: `dge` is updated with normalisation factors to adjust the counts. 
logcpm <- cpm(dge, log = TRUE) # `cpm`: Computes counts per million (CPM). `log = TRUE`: Indicates that the CPM values should be log2-transformed. Result: logcpm is a matrix of log-transformed CPM values, which can be used for downstream analysis.

# Design matrix
design <- model.matrix(~0 + groups$group) # Creates a design matrix without an intercept, where each group gets its own column. This design matrix is crucial for setting up linear models for differential expression analysis, specifying how each sample is assigned to the experimental groups.
# `model.matrix`: Creates a design matrix for linear modelling, used in differential expression analysis to specify the experimental design. `~0`: Indicates no intercept in the model, meaning each level of the group factor will have its own column without including a column for the intercept. `groups$group`: The factor variable representing the group assignments of the samples (e.g., "Unt", "TT", "TMZ"). Result: `design` is a design matrix where each column corresponds to one of the groups, and each row corresponds to a sample.

colnames(design) <- levels(groups$group) # Sets the column names of the design matrix to the group names, ensuring clarity in the matrix representation.
# Assigns meaningful column names to the design matrix. `design` is the design matrix created in the previous step. `levels(groups$group)` is the factor levels of groups$group, which are the group names ("Unt", "TT", "TMZ"). Result: The columns of the design matrix are named according to the group levels, making it easier to interpret.

# Notes: Differential expression analysis often uses linear models to quantify the relationship between gene expression levels and experimental conditions or groups. Packages like limma and edgeR rely on design matrices to fit linear models and perform differential expression analysis.

# Fit the linear model
fit <- lmFit(logcpm, design) # `lmFit` (limma package): Fits a linear model to the expression data for each gene. `logcpm`: The log-transformed CPM matrix, which represents normalised expression values. `design`: The design matrix created earlier, specifying the experimental groups for each sample. Result: Produces a `fit` object with the results of the linear models for each gene.

# Define contrasts
contrasts <- makeContrasts(
  TT_vs_Unt = TT - Unt,
  TMZ_vs_Unt = TMZ - Unt,
  levels = design
) # `makeContrasts` (limma): Defines contrasts for comparing different experimental conditions in a linear model. `TT_vs_Unt`: Specifies a contrast comparing the "TT" group to the "Unt" (untreated). `levels = design`: Specifies the design matrix (design) that contains the levels of the experimental conditions. Result: `contrasts`: A matrix of contrast vectors that will be used for differential expression analysis.

# Fit contrasts
fit <- contrasts.fit(fit, contrasts) # `contrasts.fit` (limma): Applies the defined contrasts to the fitted linear models. `fit`: The fitted linear model object from `lmFit`. `contrasts`: The contrasts matrix created by `makeContrasts`. Result: Updates the fit object to incorporate the specified contrasts, allowing for comparison between conditions.
fit <- eBayes(fit) # `eBayes`: Applies empirical Bayes moderation. Empirical Bayes moderation shrinks variance estimates towards a global average, making them more stable. Empirical Bayes moderation is a standard and well-validated approach in differential expression analysis.

# Voom transformation
voom_data <- voom(dge, design) # `voom`: Applies a more sophisticated log2-CPM transformation with associated precision weights, preparing it for better linear modeling. `dge`: The DGEList object containing normalised count data. `design`: The design matrix specifying the experimental groups. Result: `voom_data` containing the normalised log2-CPM values after voom transformation (specifically, stored in voom_data$E).

# QC Step 1.1: 2D PCA Plot
pca <- prcomp(t(voom_data$E)) # `prcomp`: Performs PCA on the transposed (t) voom-transformed expression data (voom_data$E), reducing dimensionality.
# Purpose of PCA:
    #1) Identify Patterns: Detects major sources of variation in the data.
    #2) Check for Batch Effects: Reveals if samples cluster by biological condition or unwanted technical variation.
    #3) Visualise Data Structure: Helps visualise high-dimensional data in 2D or 3D plots.
# Reduces data dimensions by projecting it onto principal components, each capturing decreasing amounts of variance (PC1 > PC2 > PC3, etc.).

pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = groups$group)
# This extracts PCA results for plotting (by creating a data frame for plotting the first two principal components). 
# `pca` is the result of the prcomp function, which performs PCA. 
# `pca$x` contains the principal component scores for each sample. PC1 and PC2 are the first and second, respectively, principal component scores for all samples. 
# The space in `pca$x[,1]` indicates all rows and the first column.
# `condition` is the experimental group labels for each sample from groups$group.
# PCA reduces data dimensions by projecting it onto principal components, each capturing decreasing amounts of variance (PC1 > PC2 > PC3, etc.)

percentVar <- round(100 * summary(pca)$importance[2, 1:2])
# `summary(pca)$importance`: Extracts the proportion of variance explained by the principal components. It multiplies by 100 to convert to percentage and rounds to the nearest whole number for the first two components. Result: `percentVar` - A vector containing the percentage of variance explained by PC1 and PC2.

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste("PC1: ", percentVar[1], "% variance")) +
    ylab(paste("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of RNA-seq data") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text_repel(aes(label = rownames(pca_data)))

# `ggplot`: Initialises a ggplot object using pca_data. Aesthetics: Maps PC1 to x-axis, PC2 to y-axis, and colours points by condition.
# `geom_point(size = 3)`: Adds points of size 3 to the plot. 
# `xlab` and `ylab`: Labels x-axis and y-axis with the percentage of variance explained by PC1 and PC2.
# `ggtitle`: Title of the  graph
# `theme_minimal`: Applies a minimal theme for a clean plot appearance.
# `legend.position`: Right
# `geom_text_repel(aes(label = rownames(pca_data)))`: Adds labels to the points, avoiding overlap, using row names of pca_data.

# QC Step 1.2: 3D PCA Plot (Optional)
transposed_data <- t(voom_data$E)
pca <- prcomp(transposed_data, scale. = TRUE)
scores <- pca$x[, 1:3]
scores <- apply(scores, 2, function(x) x / max(abs(x)))

# `transposed_data`: Transposes the voom-transformed expression data matrix (Rows become columns and columns become rows, to prepare for 3D plot). 
# `prcomp`: Performs PCA on the transposed data with scaling. Result: `pca` object containing principal components.
# `scores <- pca$x[, 1:3]`: Extracts the scores (principal component values) for the first three principal components (PC1, PC2, PC3).
# `apply`: Normalises each column (PC1, PC2, PC3) by dividing by the maximum absolute value in that column. Result: `scores`: Scaled scores where each principal component's values are normalised.


groups <- factor(c(rep("Unt", 3), rep("TT", 3), rep("TMZ", 3)))
colors <- c("blue", "green", "red")
group_colors <- colors[groups]

# `groups`: Creates a factor vector with group labels "Unt", "TT", and "TMZ", repeated 3 times each.
# `colors`: Defines a vector of colors corresponding to the groups.
# `group_colors`: Assigns colors to each sample based on their group, resulting in a vector of colors matching the `groups` factor.


pca_data <- data.frame(PC1 = scores[, 1], PC2 = scores[, 2], PC3 = scores[, 3], 
                       Sample = colnames(voom_data$E), GroupColor = group_colors)
# `pca_data`: Creates a data frame containing the first three principal component scores (PC1, PC2, PC3), sample names, and group colors.

plot <- plot_ly(data = pca_data, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'markers',
                marker = list(size = 10, color = ~GroupColor), text = ~Sample) %>%
    layout(title = "3D PCA Plot",
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
# `plot_ly`: Creates an interactive 3D scatter plot using the `plotly` package.
# `x = ~PC1, y = ~PC2, z = ~PC3`: Maps PC1, PC2, and PC3 to the x, y, and z axes respectively. The `~` symbol (tilde) indicates formula syntax in R, mapping the column names to the plot's axes.
# `marker = list(size = 10, color = ~GroupColor)`: Sets the size and color of the markers. Colors are assigned based on the `GroupColor` column.
# `text = ~Sample`: Adds sample names as text labels to the markers.
# `%>%`: The pipe operator from the `magrittr` package, used to pass the plotly object to the next function.
# `layout`: Sets the layout of the plot, including the title and axis labels.

plot

saveWidget(plot, '3D_Plot_PCA.html', selfcontained = TRUE)
# `saveWidget`: Saves the interactive plot as an HTML file using the `htmlwidgets` package.
# `plot`: The plotly object created earlier.
# `'3D_Plot_PCA.html'`: The filename for the saved HTML file.
# `selfcontained = TRUE`: Indicates that the HTML file should be self-contained, embedding all necessary resources within the file.


# QC Step 2: Sample Correlation
sample_correlation <- cor(voom_data$E, method = "spearman")
pheatmap(sample_correlation, clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", 
         color = colorRampPalette(c("#3399FF", "white", "#FF3333"))(50), display_numbers = TRUE)
# `sample_correlation`: Computes the Spearman correlation matrix for the voom-transformed expression data.
# `cor`: Function to calculate the correlation matrix.
# `voom_data$E`: The log2-CPM expression data.
# `method = "spearman"`: Specifies the Spearman correlation method, which measures rank correlations (without assumption of normal distribution)
# `pheatmap`: Function from the `pheatmap` package to create a heatmap of the correlation matrix.
# `sample_correlation`: The correlation matrix computed earlier.
# `clustering_distance_rows = "euclidean"`: Uses Euclidean distance for hierarchical clustering of rows (commonly used in gene expression clustering)
# `clustering_distance_cols = "euclidean"`: Uses Euclidean distance for hierarchical clustering of columns.
# `clustering_method = "complete"`: Uses complete linkage for hierarchical clustering.
# `color = colorRampPalette(c("#3399FF", "white", "#FF3333"))(50)`: Defines a color gradient from blue to white to red for the heatmap.
# `display_numbers = TRUE`: Displays the correlation values in the heatmap cells.

print(head(sample_correlation))

# QC Step 3: Normalised counts for boxplot
# wesanderson FantasticFox1
library(wesanderson)
palette_wes <- wes_palette("FantasticFox1", n = ncol(voom_data$E), type = "continuous")
# `palette_wes`: Generates a colour palette from the "FantasticFox1" palette in the `wesanderson` package.
# `wes_palette("FantasticFox1", n = ncol(voom_data$E), type = "continuous")`: Creates a continuous colour palette with a number of colours equal to the number of columns (samples) in the `voom_data$E` matrix.

boxplot(voom_data$E, las = 2, col = palette_wes, main = "Boxplot of Normalised Counts (FantasticFox1)",
        ylab = "Log2 counts per million", xlab = "Samples",
        outpch = 19, outcol = "black", outcex = 0.5) # Show other palettes
# `boxplot`: Creates a boxplot of the voom-transformed expression data.
# `voom_data$E`: The log2-CPM expression data.
# `las = 2`: Rotates the x-axis labels to be perpendicular to the axis.
# `col = palette_wes`: Colours the boxplots using the `FantasticFox1` palette.
# `main`: Sets the main title of the plot to "Boxplot of Normalied Counts (FantasticFox1)".
# `ylab = "Log2 counts per million"`: Labels the y-axis as "Log2 counts per million".
# `xlab = "Samples"`: Labels the x-axis as "Samples".
# `outpch = 19`: Sets the plotting character for outliers to solid circles.
# `outcol = "black"`: Colours the outliers in black.
# `outcex = 0.5`: Sets the size of the outlier points.

# Alternative style with ggplot2 (Optional)
normalised_counts <- voom_data$E # Stores the voom-transformed expression data matrix.
df <- as.data.frame(normalised_counts) # Converts the normalised_counts matrix to a data frame.
samples <- colnames(df) # Retrieves the column names (sample names) from the data frame.
df$Gene <- rownames(df) # Adds a new column to the data frame with the row names (gene names) from the original matrix.
df_long <- melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression") # Converts the data frame from wide format to long format using the `melt` function from the `reshape2` package.
# `id.vars = "Gene"`: Specifies that the "Gene" column should remain as identifier variables.
# `variable.name = "Sample"`: Names the new column that will hold the former column names (samples).
# `value.name = "Expression"`: Names the new column that will hold the expression values.

# Box plot of normalised counts
p <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "Log-normalised counts", title = "Expression Values by Sample")

# `ggplot`: Initialises a ggplot object using `df_long` as the data source.
# `aes(x = Sample, y = Expression, fill = Sample)`: Maps the x-axis to the `Sample` column, the y-axis to the `Expression` column, and uses `Sample` to fill the boxplots with colours.
# `geom_boxplot()`: Adds a boxplot layer to the plot, displaying the distribution of expression values for each sample.
# `theme(axis.text.x = element_text(angle = 45, hjust = 1))`: Customizes the plot's theme by rotating the x-axis text labels 45 degrees and adjusting their horizontal justification for better readability.
# `labs(...)`: Adds axis labels and a title to the plot.

print(p)

# Alternative box plot without outliers using 'outlier.shape = NA' argument
p_alt1 <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, fatten = 1.5, color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank()) +
    labs(x = "Sample", y = "Log-normalised counts", title = "Expression Values by Sample")

print(p_alt1)

# QC Step 4: Mean-Variance plot
mean_counts <- rowMeans(data.filt)
# `mean_counts`: Computes the mean expression level for each gene across all samples.
# `rowMeans(data.filt)`: Calculates the mean of each row (gene) in the filtered expression data matrix `data.filt`.

var_counts <- apply(data.filt, 1, var)
# `var_counts`: Computes the variance of expression levels for each gene across all samples.
# `apply(data.filt, 1, var)`: Applies the `var` function to each row (gene) in the filtered expression data matrix `data.filt`.

df_before_voom <- data.frame(Mean = mean_counts, Variance = var_counts)
# `df_before_voom`: Creates a data frame with the mean and variance of expression levels for each gene.
# `data.frame(Mean = mean_counts, Variance = var_counts)`: Constructs a data frame with two columns: `Mean` and `Variance`.

df_before_voom <- df_before_voom[mean_counts > 0 & var_counts > 0, ]
# `df_before_voom`: Filters the data frame to include only genes with positive mean and variance.
# `[mean_counts > 0 & var_counts > 0, ]`: Subsets the data frame, keeping only rows where both mean and variance are greater than zero.

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


### Additional checks ###
# Library size visualisation. The "library size" for each sample refers to the total number of reads (or fragments) that were sequenced for that sample. This is calculated by summing up all the raw read counts across all genes for each sample. This is an important quality control step to check for any significant differences in sequencing depth between samples.
library_sizes <- colSums(data.filt)
# `library_sizes`: Computes the total count of reads for each sample.
# `colSums(data.filt)`: Calculates the sum of each column (sample) in the filtered expression data matrix `data.filt`.

df_library_sizes <- data.frame(Sample = colnames(data.filt), LibrarySize = library_sizes)
# `df_library_sizes`: Creates a data frame with sample names and their corresponding library sizes.
# `data.frame(Sample = colnames(data.filt), LibrarySize = library_sizes)`: Constructs a data frame with two columns: `Sample` and `LibrarySize`.

ggplot(df_library_sizes, aes(x = Sample, y = LibrarySize)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Library Sizes Across Samples", x = "Sample", y = "Library Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# `ggplot`: Initialises a ggplot object using `df_library_sizes` as the data source.
# `aes(x = Sample, y = LibrarySize)`: Maps the x-axis to the `Sample` column and the y-axis to the `LibrarySize` column.
# `geom_bar(stat = "identity")`: Adds a bar plot layer to the plot, using the raw values for bar heights.
# `theme_minimal()`: Applies a minimal theme for a clean plot.
# `labs(...)`: Adds axis labels and a title to the plot.
# `theme(axis.text.x = element_text(angle = 45, hjust = 1))`: Customises the plot's theme by rotating the x-axis text labels 45 degrees and adjusting their horizontal justification for better readability.

    
# Extract BH-adjusted p-values
fit2 <- eBayes(fit)
# `fit2`: Applies empirical Bayes moderation to the linear model fits.
# `eBayes(fit)`: Uses the `eBayes` function from the `limma` package to compute moderated statistics for the linear model fits stored in `fit`.

tt_vs_untreated_pvals <- topTable(fit2, coef = "TT_vs_Unt", number = Inf, sort.by = "P")
# `tt_vs_untreated_pvals`: Extracts the top table of differentially expressed genes for the contrast "TT vs. Untreated".
# `topTable(fit2, coef = "TT_vs_Unt", number = Inf, sort.by = "P")`: Uses the `topTable` function from the `limma` package to get the results for the specified contrast.
# `fit2`: The object containing the empirical Bayes moderated statistics.
# `coef = "TT_vs_Unt"`: Specifies the contrast to extract the results.
# `number = Inf`: Requests all genes (no limit on the number).
# `sort.by = "P"`: Sorts the results by p-value.

tmz_vs_untreated_pvals <- topTable(fit2, coef = "TMZ_vs_Unt", number = Inf, sort.by = "P")
# `tmz_vs_untreated_pvals`: Extracts the top table of differentially expressed genes for the contrast "TMZ vs. Untreated".
# `topTable(fit2, coef = "TMZ_vs_Unt", number = Inf, sort.by = "P")`: Uses the `topTable` function from the `limma` package to get the results for the specified contrast.
# `fit2`: The object containing the empirical Bayes moderated statistics.
# `coef = "TMZ_vs_Unt"`: Specifies the contrast to extract the results for.
# `number = Inf`: Requests all genes (no limit on the number).
# `sort.by = "P"`: Sorts the results by p-value.

# Plot raw BH-adjusted p-values for inspection (TT vs untreated)
hist(tt_vs_untreated_pvals$adj.P.Val, breaks = 50, main = "BH-Adjusted p-values (TT vs Untreated)", xlab = "BH-Adjusted p-value", ylab = "Frequency", col = "lightgreen")
# `hist`: Creates a histogram to visualise the distribution of BH-adjusted p-values.
# `tt_vs_untreated_pvals$adj.P.Val`: Specifies the BH-adjusted p-values for the "TT vs Untreated" contrast.
# `breaks = 50`: Sets the number of bins in the histogram to 50.

# Plot raw BH-adjusted p-values for inspection (TMZ vs untreated)
hist(tmz_vs_untreated_pvals$adj.P.Val, breaks = 50, main = "BH-Adjusted p-values (TMZ vs Untreated)", xlab = "BH-Adjusted p-value", ylab = "Frequency", col = "lightcoral")

# Extract log2 fold changes
tt_vs_untreated_lfc <- tt_vs_untreated_pvals$logFC
# `tt_vs_untreated_lfc`: Extracts log2 fold changes for the "TT vs Untreated" contrast.
# `tt_vs_untreated_pvals$logFC`: Retrieves the log2 fold change values from the differential expression results for the "TT vs Untreated" contrast.
tmz_vs_untreated_lfc <- tmz_vs_untreated_pvals$logFC

# Plot distribution of log2 fold changes (TMZ vs untreated)
hist(tmz_vs_untreated_lfc, breaks = 50, main = "Log2 Fold Changes (TMZ vs Untreated)", xlab = "Log2 Fold Change", ylab = "Frequency", col = "lightcoral")

# Plot distribution of log2 fold changes (TT vs untreated)
hist(tt_vs_untreated_lfc, breaks = 50, main = "Log2 Fold Changes (TT vs Untreated)", xlab = "Log2 Fold Change", ylab = "Frequency", col = "lightgreen")

# Plot gene variances (not needed to show)
gene_variances <- apply(voom_data$E, 1, var)
# `gene_variances`: Computes the variance of expression levels for each gene.
# `apply(voom_data$E, 1, var)`: Applies the `var` function to each row (gene) in the voom-transformed expression data matrix `voom_data$E`.
# `1`: Indicates that the function should be applied to rows (genes).
# `var`: The function used to calculate the variance.
hist(gene_variances, breaks = 50, main = "Distribution of Gene Variances", xlab = "Variance", ylab = "Frequency", col = "lightblue")
abline(h = 0, col = "blue")



### DEGs List ###
# Obtain results for TT vs Untreated
tt_vs_unt <- topTable(fit, coef = "TT_vs_Unt", number = Inf, adjust.method = "BH")
tt_vs_unt$Gene <- rownames(tt_vs_unt)
# `tt_vs_unt`: Extracts the top table of differentially expressed genes for the "TT vs Untreated" contrast.
# `topTable(fit, coef = "TT_vs_Unt", number = Inf, adjust.method = "BH")`: Uses the `topTable` function from the `limma` package to get the results for the specified contrast.
# `fit`: The object containing the fitted linear models.
# `coef = "TT_vs_Unt"`: Specifies the contrast to extract the results for.
# `number = Inf`: Requests all genes (no limit on the number).
# `adjust.method = "BH"`: Adjusts p-values using the Benjamini-Hochberg method.

# Obtain results for TMZ vs Untreated
tmz_vs_unt <- topTable(fit, coef = "TMZ_vs_Unt", number = Inf, adjust.method = "BH")
tmz_vs_unt$Gene <- rownames(tmz_vs_unt)

# Reorder columns to put Gene as the first column
tt_vs_unt <- tt_vs_unt[, c("Gene", setdiff(colnames(tt_vs_unt), "Gene"))]
tmz_vs_unt <- tmz_vs_unt[, c("Gene", setdiff(colnames(tmz_vs_unt), "Gene"))]
# `tt_vs_unt`: Reorders the columns of the `tt_vs_unt` data frame to place the "Gene" column first.
# `tt_vs_unt[, c("Gene", setdiff(colnames(tt_vs_unt), "Gene"))]`: Uses column indexing to reorder columns.
# `c("Gene", setdiff(colnames(tt_vs_unt), "Gene"))`: Creates a vector with "Gene" as the first element followed by all other column names except "Gene".

# Save all genes to Excel
write.xlsx(list(TT_vs_Untreated = tt_vs_unt, TMZ_vs_Untreated = tmz_vs_unt), file = "All_Genes_TT_TMZ_vs_Untreated_v2.xlsx")


### Volcano Plot ###
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

# Alternative: Manually define the genes
hallmark_genes_TGFb <- c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP") 

# Filter the limma results to include only genes in the hallmark_genes list
hallmark_results_TT_vs_Untreated <- tt_vs_unt[tt_vs_unt$Gene %in% hallmark_genes, ]
hallmark_results_TMZ_vs_Untreated <- tmz_vs_unt[tmz_vs_unt$Gene %in% hallmark_genes, ]
# `hallmark_results_TT_vs_Untreated`: Filters the `tt_vs_unt` data frame to include only rows where the "Gene" column is in the `hallmark_genes` list.
# `tt_vs_unt$Gene %in% hallmark_genes`: Creates a logical vector indicating whether each gene in `tt_vs_unt$Gene` is present in `hallmark_genes`.
# `tt_vs_unt[tt_vs_unt$Gene %in% hallmark_genes, ]`: Subsets `tt_vs_unt` to include only the rows where the "Gene" is in `hallmark_genes`.
# `%in%`: An operator that checks if elements of the left operand are present in the right operand, returning a logical vector. In this context, it checks if each gene in the "Gene" column is in the hallmark_genes list.

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

# `volcano_TT`: Initialises a ggplot object to create a volcano plot for the "TT vs Untreated" results.
# `ggplot(hallmark_results_TT_vs_Untreated, aes(x = logFC, y = -log10(adj.P.Val), color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")))`: Sets up the plot with the data and aesthetics.
# `hallmark_results_TT_vs_Untreated`: Uses the filtered results for the "TT vs Untreated" contrast that include only hallmark genes.
# `aes(x = logFC, y = -log10(adj.P.Val))`: Maps the log fold change (`logFC`) to the x-axis and the negative log10 adjusted p-value (`-log10(adj.P.Val)`) to the y-axis.
# `color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not significant")`: Colors the points based on their significance, where points with an absolute log fold change greater than 1 and an adjusted p-value less than 0.05 are labeled "Significant", and others are labeled "Not significant".

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

# With gene labels
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

# Check the number of genes
length(hallmark_genes)

# Filter the limma results to include only genes in the hallmark_genes list
hallmark_results_TT_vs_Untreated <- tt_vs_unt[tt_vs_unt$Gene %in% hallmark_genes, ]
hallmark_results_TMZ_vs_Untreated <- tmz_vs_unt[tmz_vs_unt$Gene %in% hallmark_genes, ]

# Define colors
colors <- c("Not significant" = "grey", "Significant" = "navyblue")

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


### Hallmark TNFa ###
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

# Alternative: Manually define the hallmark genes 
hallmark_genes_TNFa <- c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIGI","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")

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


### GSEA ###
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)

# Load the hallmark gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
# `hallmark_gene_sets`: Loads the hallmark gene sets for Hs from the MSigDB database.
# `msigdbr(species = "Homo sapiens", category = "H")`: Uses the `msigdbr` package to retrieve hallmark gene sets for the specified species and category.

hallmark_tgf_beta <- hallmark_gene_sets[hallmark_gene_sets$gs_name == "HALLMARK_TGF_BETA_SIGNALING", ]$gene_symbol
# `hallmark_tgf_beta`: Extracts the gene symbols for the "HALLMARK_TGF_BETA_SIGNALING" gene set.
# `hallmark_gene_sets[hallmark_gene_sets$gs_name == "HALLMARK_TGF_BETA_SIGNALING", ]`: Filters the `hallmark_gene_sets` data frame to include only rows where the `gs_name` column matches "HALLMARK_TGF_BETA_SIGNALING".
# `$gene_symbol`: Selects the `gene_symbol` column from the filtered data frame.

# Alternative: Manually define the hallmark genes 
hallmark_tgf_beta <- c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP") 

# Create rank list for TT vs Untreated
tt_vs_unt_rank <- tt_vs_unt$logFC
names(tt_vs_unt_rank) <- tt_vs_unt$Gene
tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)


# Create rank list for TMZ vs Untreated
tmz_vs_unt_rank <- tmz_vs_unt$logFC
# `tt_vs_unt_rank`: Extracts the log fold change values from the `tt_vs_unt` data frame.
# `tt_vs_unt$logFC`: Selects the log fold change column from the `tt_vs_unt` data frame.
names(tmz_vs_unt_rank) <- tmz_vs_unt$Gene
# `names(tt_vs_unt_rank)`: Assigns gene names to the log fold change values.
# `tt_vs_unt$Gene`: Uses the gene names from the `tt_vs_unt` data frame to name the elements of `tt_vs_unt_rank`.
tmz_vs_unt_rank <- sort(tmz_vs_unt_rank, decreasing = TRUE)
# `tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)`: Sorts the log fold change values in decreasing order.
# `sort(tt_vs_unt_rank, decreasing = TRUE)`: Arranges the log fold change values from highest to lowest.

# Set seed and ensure proper ranking
set.seed(1) # Sets the seed for random number generation to ensure reproducibility of results.
tt_vs_unt_rank <- tt_vs_unt$logFC + rnorm(length(tt_vs_unt$logFC), sd = 1e-5)
# `tt_vs_unt_rank`: Adds a small random noise to the log fold change values to avoid ties.
# `tt_vs_unt$logFC`: The original log fold change values.
# `rnorm(length(tt_vs_unt$logFC), sd = 1e-5)`: Generates random normal values with a standard deviation of 1e-5, the same length as the log fold change values.
names(tt_vs_unt_rank) <- tt_vs_unt$Gene
# `names(tt_vs_unt_rank)`: Assigns gene names to the log fold change values with added noise.
# `tt_vs_unt$Gene`: Uses the gene names from the `tt_vs_unt` data frame to name the elements of `tt_vs_unt_rank`.
tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)
# `tt_vs_unt_rank <- sort(tt_vs_unt_rank, decreasing = TRUE)`: Sorts the log fold change values with added noise in decreasing order.
# `sort(tt_vs_unt_rank, decreasing = TRUE)`: Arranges the log fold change values from highest to lowest.

tmz_vs_unt_rank <- tmz_vs_unt$logFC + rnorm(length(tmz_vs_unt$logFC), sd = 1e-5)
names(tmz_vs_unt_rank) <- tmz_vs_unt$Gene
tmz_vs_unt_rank <- sort(tmz_vs_unt_rank, decreasing = TRUE)

# Run GSEA for TT vs Untreated
fgsea_res_tt <- fgseaMultilevel(pathways = list(HALLMARK_TGF_BETA_SIGNALING = hallmark_tgf_beta), 
                                stats = tt_vs_unt_rank)
# `fgsea_res_tt`: Stores the results of the Gene Set Enrichment Analysis (GSEA) for the "TT vs Untreated" contrast.
# `fgseaMultilevel(...)`: Runs the multilevel GSEA analysis using the `fgsea` package.
# `pathways = list(...)`: Specifies the gene sets to test, in this case, the "HALLMARK_TGF_BETA_SIGNALING" gene set.
# `stats = tt_vs_unt_rank`: Provides the ranked list of log fold changes for the "TT vs Untreated" contrast.

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
        
        gseaPlot <- plotEnrichment(pathway = pathway,stats = stats) + ggtitle(title) + theme_minimal() + labs(y = "Running Enrichment Score", x = "Rank in Ordered Gene List") + annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "\nFDR =", round(pval, 4)), hjust = 1.1, vjust = 1.1, size = 5, color = "blue", fontface = "bold")
        
        print(gseaPlot)
    } else {message("No significant enrichment found for ", title)}
}

        # `fgseaRes <- fgseaRes[order(fgseaRes$padj), ]`: Sorts the GSEA results by adjusted p-value in ascending order.
        # `pathway`: Uses the full "HALLMARK_TGF_BETA_SIGNALING" gene set for plotting.
        # `enrichmentScore`: Extracts the enrichment score for the top-ranked pathway.
        # `nes`: Extracts the normalised enrichment score for the top-ranked pathway.
        # `pval`: Extracts the adjusted p-value for the top-ranked pathway.
        
        # `plotEnrichment(pathway = pathway, stats = stats)`: Plots the enrichment score for the specified pathway and ranked statistics.
        # `ggtitle(title)`: Adds a title to the plot.
        # `theme_minimal()`: Applies a minimal theme to the plot.
        # `labs(...)`: Labels the y-axis and x-axis.
        # `annotate(...)`: Adds text annotations to the plot with the NES and adjusted p-value.
        
        # `print(gseaPlot)`: Displays the plot.
        # `message("No significant enrichment found for ", title)`: Prints a message if no significant enrichment is found.

# Plot GSEA results for TT vs Untreated
plotGseaTable(fgsea_res_tt, tt_vs_unt_rank, "TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING")
# `plotGseaTable(...)`: Calls the `plotGseaTable` function to plot the GSEA results for the "TT vs Untreated" contrast.
# `fgsea_res_tt`: The GSEA results for the "TT vs Untreated" contrast.
# `tt_vs_unt_rank`: The ranked list of log fold changes for the "TT vs Untreated" contrast.
# `"TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING"`: The title for the plot.

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
                    annotate("text", x = Inf, y = Inf, label = paste("NES =", round(nes, 2), "\nFDR =", formatted_pval), 
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
                    geom_line(data = ggplot_build(gseaPlot)$data[[1]], aes(x = x, y = y), color = "#009999", linewidth = 1)
        
        print(gseaPlot)
    } else {
        message("No significant enrichment found for ", title)
    }
}

# Plot GSEA results for TT vs Untreated
plotGseaTable(fgsea_res_tt, tt_vs_unt_rank, "TT vs Untreated: HALLMARK_TGF_BETA_SIGNALING")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable(fgsea_res_tmz, tmz_vs_unt_rank, "TMZ vs Untreated: HALLMARK_TGF_BETA_SIGNALING")



### GSEA TNFa ###
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
plotGseaTable_TNFa(fgsea_res_tt_TNFa, tt_vs_unt_rank_TNFa, "TT vs Untreated: HALLMARK_TNFA_SIGNALING_VIA_NFKB")

# Plot GSEA results for TMZ vs Untreated
plotGseaTable_TNFa(fgsea_res_tmz_TNFa, tmz_vs_unt_rank_TNFa, "TMZ vs Untreated: HALLMARK_TNFA_SIGNALING_VIA_NFKB")