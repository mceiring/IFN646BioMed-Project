# Load required libraries
library(DESeq2)

# Import counts.txt into R
count_data <- read.table("counts.txt", header=TRUE, row.names=1)

## Remove columns that are not sample counts (like Chr, Start, End, Strand, Length)
count_data <- count_data[,6:ncol(count_data)]

# Create a sample information table (sampleTable). This is an example and should be customized based on your experiment design.
#sample_table <- data.frame(
#  sampleName = colnames(count_data),
#  fileName = colnames(count_data),
#  condition = c("control", "treatment", ...), # Replace with your actual conditions
#  row.names = colnames(count_data)
#)

## Create the DESeqDataSet
#dds <- DESeqDataSetFromMatrix(
#  countData = count_data,
#  colData = sample_table,
#  design = ~ condition # Replace with your actual design formula
#)