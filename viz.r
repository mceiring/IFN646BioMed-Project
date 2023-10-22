library(ggplot2)
library(dplyr)
library(readr)

# List all the CSV files in the 'results' directory
file_paths <- list.files(path="results", pattern="*.csv", full.names=TRUE)

# Read all CSV files into a list of data frames
datasets <- lapply(file_paths, read_csv)

# Assuming you want similar visualizations for all datasets, we can loop through each dataset
for(i in seq_along(datasets)) {
  
  current_data <- datasets[[i]]
  
  # Load the CSV data
  df <- read_csv("data.csv")
  
  # Visualization 1: Bar Chart for Log2 Fold Change
  ggplot(df, aes(x=Identifier, y=log2FoldChange)) + 
    geom_bar(stat="identity", fill="skyblue", color="black") + 
    labs(title="Log2 Fold Change by Gene Identifier", y="Log2 Fold Change") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Visualization 2: Volcano Plot
  ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue), color=padj < 0.05)) + 
    geom_point() + 
    scale_color_manual(values=c("gray", "red")) +
    labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10(P-value)") + 
    theme_minimal()
  
  # Visualization 3: Error Bar Chart for Log2 Fold Change with Standard Error
  ggplot(df, aes(x=Identifier, y=log2FoldChange)) + 
    geom_bar(stat="identity", fill="skyblue", color="black") + 
    geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=0.2) +
    labs(title="Log2 Fold Change with Standard Error by Gene Identifier", y="Log2 Fold Change") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Visualization 4: Bar Chart for Base Mean Values
  ggplot(df, aes(x=Identifier, y=baseMean)) + 
    geom_bar(stat="identity", fill="skyblue", color="black") + 
    labs(title="Base Mean Value by Gene Identifier", y="Base Mean Value") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  }
