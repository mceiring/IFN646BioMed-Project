# ~~~ Setup ~~~ #
# Load required libraries
library(DESeq2)
library(ggplot2)
library(VennDiagram)

# ~~~ Data ~~~ #
# Category lists
conditions = c("control", "treatment")
error_levels = c("low", "medium", "high")
replicates = c("01", "02", "03")

# Inisialise vars
count_list_alignment_free = list()
count_list_alignment_based = list()
target_ids_alignment_free = NULL
target_ids_alignment_based = NULL

# Read in alignment-free data
for(error in error_levels){
  for(rep in replicates){
    for(cond in conditions){
      path = paste0("alignment_free/", cond, "_", error, "_", rep, "/abundance.tsv")
      data_alignment_free = read.table(path, header=TRUE, sep="\t")

      if (is.null(target_ids_alignment_free)){
        target_ids_alignment_free = data_alignment_free$"target_id"
      }

      count_list_alignment_free[[paste(cond, error, rep, sep="_")]] = data_alignment_free$est_counts
    }
  }
}

# Read in alignment-based data
for(error in error_levels){
  for(rep in replicates){
    for(cond in conditions){
      path = paste0("alignment_based/", cond, "_", error, "_", rep, "_counts.txt")
      data_alignment_based = read.table(path, header=TRUE, sep="\t", skip=1, check.names=FALSE)

      if (is.null(target_ids_alignment_based)){
        target_ids_alignment_based = data_alignment_based$"Chr"
      }
      
      # Dynamically get the last column for counts (file path but .sam)
      count_column = names(data_alignment_based)[ncol(data_alignment_based)]
      count_list_alignment_based[[paste(cond, error, rep, sep="_")]] = data_alignment_based[[count_column]]
    }
  }
}

# Convert the list to a matrix
count_matrix_alignment_free = do.call(cbind, count_list_alignment_free)
count_matrix_alignment_based = do.call(cbind, count_list_alignment_based)

# Add target_ids as rownames
rownames(count_matrix_alignment_free) = target_ids_alignment_free
rownames(count_matrix_alignment_based) = target_ids_alignment_based

# tximport not working, must round for alignment-free; doing also for alignment-based (not sure if float)
count_matrix_alignment_free = round(count_matrix_alignment_free)  
count_matrix_alignment_based = round(count_matrix_alignment_based)

# ~~~ Prepare for Differential Expression Analysis ~~~ #
condition_vector = rep(rep(conditions, each=length(replicates)), times=length(error_levels))
error_vector = rep(error_levels, each=length(replicates)*length(conditions))
replicate_vector = rep(replicates, times=length(conditions)*length(error_levels))

col_data = data.frame(condition=condition_vector, error=error_vector, replicate=replicate_vector)

# ~~~ Differential Expression Analysis ~~~ #
results_alignment_free = list()
results_alignment_based = list()

for(error in error_levels){
  # Alignment-Free Analysis
  sub_dds_alignment_free = count_matrix_alignment_free[,error_vector == error]
  sub_col_data_alignment_free = col_data[error_vector == error,]
  
  dds_alignment_free = DESeqDataSetFromMatrix(countData = sub_dds_alignment_free,
                               colData = sub_col_data_alignment_free,
                               design = ~ condition)
  
  dds_alignment_free = DESeq(dds_alignment_free)
  res_alignment_free = results(dds_alignment_free)
  results_alignment_free[[error]] = res_alignment_free
  
  # Alignment-Based Analysis
  sub_dds_alignment_based = count_matrix_alignment_based[,error_vector == error]
  sub_col_data_alignment_based = col_data[error_vector == error,]
  
  dds_alignment_based = DESeqDataSetFromMatrix(countData = sub_dds_alignment_based,
                                               colData = sub_col_data_alignment_based,
                                               design = ~ condition)
  
  dds_alignment_based = DESeq(dds_alignment_based)
  res_alignment_based = results(dds_alignment_based)
  results_alignment_based[[error]] = res_alignment_based
}

sig_genes_alignment_free = list()
sig_genes_alignment_based = list()

for(error in error_levels){
  # Alignment-Free Analysis
  res_alignment_free = results_alignment_free[[error]]
  sig_genes_alignment_free[[error]] = subset(res_alignment_free, padj < 0.05)
  
  # Alignment-Based Analysis
  res_alignment_based = results_alignment_based[[error]]
  sig_genes_alignment_based[[error]] = subset(res_alignment_based, padj < 0.05)
  
  # Writing results
  output_path_free = paste0("results/sig_genes_alignment_free_", error, ".csv")
  write.csv(sig_genes_alignment_free[[error]], file = output_path_free, row.names = TRUE)
  output_path_based = paste0("results/sig_genes_alignment_based_", error, ".csv")
  write.csv(sig_genes_alignment_based[[error]], file = output_path_based, row.names = TRUE)
}

# ~~~ Evaluation ~~~ #
# 1. Reproducibility (CoV)

# 2. Robustness to Noise (Spearman)

# 3. Alignment-Based vs Alignment-Free (Jaccard)
