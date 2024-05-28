### Codes to generate template 3D plot ###
# Set your directory
setwd("C:/Users/...") # Set your directory using forward slashes (not backslashes)

# Load required libraries
library(readxl)
library(plotly)
library(htmlwidgets)

# Load the Excel file (change the values as desired to suit your data)
file_path <- "Template_Data_for_3D_Plot.xlsx"
data <- read_excel(file_path)

# Create a 3D scatter plot (modify the colors and size of the bubbles as desired)
plot <- plot_ly(data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~GroupColor, colors = c("blue", "darkgreen", "red"), text = ~Sample,
                marker = list(size = 10)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'x-axis'),
                        yaxis = list(title = 'y-axis'),
                        zaxis = list(title = 'z-axis')))

plot # This shows the plot in RStudio (in the Viewer tab)

# Save the plot as an HTML file
saveWidget(plot, '3D_Plot_Template.html', selfcontained = TRUE)


### Codes to generate template boxplot with wesanderson palette ###
# Set your directory
setwd("C:/Users/...") # Set your directory using forward slashes (not backslashes)

# Load required libraries
library(wesanderson)
library(openxlsx)

# Load the data (modify according to your data as desired)
boxplot_data <- read.xlsx("Template_Data_for_Boxplot.xlsx")

# Transform data into a matrix format for the boxplot
voom_matrix <- matrix(boxplot_data$LogCounts, ncol = 9, byrow = TRUE)
colnames(voom_matrix) <- c("Unt1", "Unt2", "Unt3", "TT1", "TT2", "TT3", "TMZ1", "TMZ2", "TMZ3")

# Boxplot using the wesanderson palette
palette_wes <- wes_palette("FantasticFox1", n = ncol(voom_matrix), type = "continuous")
boxplot(voom_matrix, las = 2, col = palette_wes, main = "Enter the title of your graph (optional)",
        ylab = "y-axis label", xlab = "x-axis label",
        outpch = 19, outcol = "black", outcex = 0.5)

# Save the plot in PDF, JPG, PNG or other format from RStudio console (the 'Export' button) 