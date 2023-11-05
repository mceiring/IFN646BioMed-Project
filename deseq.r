# ~~~ Setup ~~~ #
# Load required libraries
library(DESeq2)
library(ggplot2)

# ~~~ Data ~~~ #
# Category lists
conditions = c("control", "treatment")
error_levels = c("low", "medium", "high")
replicates = c("01", "02", "03")

# Inisialise vars
count_list_alignment_free = list()
count_list_alignment_based = list()
target_ids_alignment_free = c()
target_ids_alignment_based = c()

# Read in alignment-free data
for(error in error_levels){
  for(replicate in replicates){
    for(cond in conditions){
      # Read in data
      path = paste0("alignment_free/", cond, "_", error, "_", replicate, "/abundance.tsv")
      data_alignment_free = read.table(path, header=TRUE, sep="\t")

      # Build out count frame and unique IDs
      count_list_alignment_free[[paste(cond, error, replicate, sep="_")]] = data_alignment_free$est_counts
      target_ids_alignment_free = union(target_ids_alignment_free, data_alignment_free$"target_id")
    }
  }
}

# Read in alignment-based data
for(error in error_levels){
  for(replicate in replicates){
    for(cond in conditions){
      # Read in data
      path = paste0("alignment_based/", cond, "_", error, "_", replicate, "_counts.txt")
      data_alignment_based = read.table(path, header=TRUE, sep="\t", skip=1, check.names=FALSE)
      
      # Define the count_column variable here, assuming the last column contains the count data
      count_column = names(data_alignment_based)[ncol(data_alignment_based)]
      
      # Split the concatenated gene IDs and replicate the counts for each gene ID
      split_ids = strsplit(as.character(data_alignment_based$"Chr"), ";")
      split_counts = mapply(rep, data_alignment_based[[count_column]], sapply(split_ids, length))

      # Flatten the list to a vector
      flat_ids = unlist(split_ids)
      flat_counts = unlist(split_counts)

      # Create a new data frame with the separated gene IDs and their counts
      separated_data = data.frame(gene_id = flat_ids, count = flat_counts)

      # Aggregate counts by gene_id in case there are duplicates
      aggregated_data = aggregate(count ~ gene_id, data = separated_data, sum)

      # Build out count frame and unique IDs
      count_list_alignment_based[[paste(cond, error, replicate, sep="_")]] = aggregated_data$count
      target_ids_alignment_based = union(target_ids_alignment_based, aggregated_data$gene_id)
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

# # ~~~ Prepare for Differential Expression Analysis ~~~ #
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
  sig_genes_alignment_free[[error]] = subset(res_alignment_free, padj < 0.25)
  
  # Alignment-Based Analysis
  res_alignment_based = results_alignment_based[[error]]
  sig_genes_alignment_based[[error]] = subset(res_alignment_based, padj < 0.25)
  
  # Writing results
  output_path_free = paste0("results/sig_genes_alignment_free_", error, ".csv")
  write.csv(sig_genes_alignment_free[[error]], file = output_path_free, row.names = TRUE)
  output_path_based = paste0("results/sig_genes_alignment_based_", error, ".csv")
  write.csv(sig_genes_alignment_based[[error]], file = output_path_based, row.names = TRUE)
}

# ~~~ Evaluation ~~~ #
# 1. Reproducibility (CoV)
# Calculate CoV for alignment-free data
cov_alignment_free <- apply(count_matrix_alignment_free, 1, function(x) sd(x) / mean(x))
mean_cov_alignment_free <- mean(cov_alignment_free, na.rm = TRUE)

# Calculate CoV for alignment-based data
cov_alignment_based <- apply(count_matrix_alignment_based, 1, function(x) sd(x) / mean(x))
mean_cov_alignment_based <- mean(cov_alignment_based, na.rm = TRUE)

write.csv(cov_alignment_free, "results/cov_alignment_free.csv", row.names = TRUE)
write.csv(cov_alignment_based, "results/cov_alignment_based.csv", row.names = TRUE)

# 2. Robustness to Noise (Spearman)
# true_counts_matrix <- read.csv("polyester_reads/fold_changes.csv")
# true_counts_matrix <- true_counts_matrix[!duplicated(true_counts_matrix$gene_id), ]
# row.names(true_counts_matrix) <- true_counts_matrix$gene_id
# true_counts_matrix <- true_counts_matrix[,-1]


# spearman_correlation_free <- apply(count_matrix_alignment_free, 1, function(x) cor(x, true_counts_matrix, method = "spearman"))
# spearman_correlation_based <- apply(count_matrix_alignment_based, 1, function(x) cor(x, true_counts_matrix, method = "spearman"))

# write.csv(spearman_correlation_free, "results/spearman_correlation_free.csv", row.names = TRUE)
# write.csv(spearman_correlation_based, "results/spearman_correlation_based.csv", row.names = TRUE)

# 3. Alignment-Based vs Alignment-Free (Jaccard)
# Calculate Jaccard Index for each error level
jaccard_index <- sapply(error_levels, function(error) {
  sig_free <- names(sig_genes_alignment_free[[error]])
  sig_based <- names(sig_genes_alignment_based[[error]])
  length(intersect(sig_free, sig_based)) / length(union(sig_free, sig_based))
})

write.csv(data.frame(Jaccard_Index = jaccard_index), "results/jaccard_index.csv", row.names = FALSE)