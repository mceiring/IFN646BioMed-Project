# ~~~ Setup ~~~ #
# Load required libraries
library(polyester)
library(Biostrings)
library(rtracklayer)
library(ShortRead)

seed_number <- 15
set.seed(seed_number) # for reproducibility

# ~~~ Data ~~~ #
# Paths
path_sequences <- "data/sequences.fna"
path_annotations <- "data/annotation.gtf"

# Loading sequence data
genome_sequence <- readDNAStringSet(path_sequences)
names(genome_sequence) <- gsub(" .*", "", names(genome_sequence))  # need to split off id

# Loading annotation data
annotations <- read.table(path_annotations, comment.char="#", quote="",
                         col.names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
                         fill = TRUE, sep="\t")

exon_annotations <- annotations[annotations$feature == "exon", ] # 35384

# Sampling annotations (cumbersome to use all 35384)
sample_indices <- sample(nrow(exon_annotations), 5000)
sample_exons <- exon_annotations[sample_indices, ]

# ~~~ Read Generation ~~~ #
# Error rates for low, medium, and high scenarios
error_rates <- list(low=0.015, medium=0.02, high=0.025)

# Define output directory
output_dir <- "polyester_reads"

# Fold changes for genes
# use random fold changes between 0.5 (downregulation) and 2 (upregulation)
# for the treatment condition compared to control.
fold_changes <- runif(n = length(sample_exons$seqname), min = 0.5, max = 2)
fold_changes_df <- data.frame(gene_id = sample_exons$seqname, fold_change = fold_changes)
write.csv(fold_changes_df, file = paste0(output_dir, "/fold_changes.csv"), row.names = FALSE)   # for comparrison later

# Set fold changes to 1 for the control group to simulate no differential expression (NONE argument not working...)
fold_changes_control <- rep(1, length(sample_exons$seqname))

# Simulate reads for control and treatment under different error scenarios
for (error_level in names(error_rates)) {
  # Control group - no fold changes, equivalent to fold change of 1 for all genes
  simulate_experiment(
    seqpath = genome_sequence,
    gtf = sample_exons,
    fold_changes = fold_changes_control,
    num_reps = 3,
    reads_per_transcript = 200,
    seq_len = 100,
    error_rate = error_rates[[error_level]],
    outdir = sprintf("%s/control_%s", output_dir, error_level),
    seqformat = "fastq"
  )
  
  # Treatment group - apply fold changes
  simulate_experiment(
    seqpath = genome_sequence,
    gtf = sample_exons,
    fold_changes = fold_changes, # Apply the fold changes here
    num_reps = 3,
    reads_per_transcript = 200,
    seq_len = 100,
    error_rate = error_rates[[error_level]],
    outdir = sprintf("%s/treatment_%s", output_dir, error_level),
    seqformat = "fastq"
  )
}

# ~~~ Read Shuffle ~~~ #
shuffle_reads_fasta <- function(file_path_1, file_path_2, shuffle_percentage=1) {
  # Read the FASTA files
  reads_1 <- readDNAStringSet(file_path_1)
  reads_2 <- readDNAStringSet(file_path_2)
  
  # Ensure that both files have the same number of reads
  if (length(reads_1) != length(reads_2)) {
    stop("The number of reads in the paired-end files do not match.")
  }
  
  # Determine the number of reads to shuffle
  num_to_shuffle <- round(length(reads_1) * shuffle_percentage)
  
  # Randomly select indices to shuffle
  indices_to_shuffle <- sample(1:length(reads_1), num_to_shuffle)
  
  # Shuffle the selected indices
  shuffled_indices <- sample(indices_to_shuffle)
  
  # Apply the shuffled order to both ends
  reads_1[indices_to_shuffle] <- reads_1[shuffled_indices]
  reads_2[indices_to_shuffle] <- reads_2[shuffled_indices]
  
  # Return the shuffled reads
  return(list(reads_1, reads_2))
}

# Loop through conditions, error scenarios, and replicates
for (condition in c("control", "treatment")) {
  for (error_level in names(error_rates)) {
    for (rep in 1:3) {
      # Define input file paths
      input_file_1 <- sprintf("%s/%s_%s/sample_%02d_1.fasta", output_dir, condition, error_level, rep)
      input_file_2 <- sprintf("%s/%s_%s/sample_%02d_2.fasta", output_dir, condition, error_level, rep)
      
      # Shuffle reads and write back to FASTA
      set.seed(seed_number*rep*2)  # adjusting seed each loop so that replicates are distinct
      shuffled_reads <- shuffle_reads_fasta(input_file_1, input_file_2)  # assumes 100% shuffle rate
      output_file_1 <- sprintf("%s/%s_%s/shuffled_sample_%02d_1.fasta", output_dir, condition, error_level, rep)
      output_file_2 <- sprintf("%s/%s_%s/shuffled_sample_%02d_2.fasta", output_dir, condition, error_level, rep)
      writeXStringSet(shuffled_reads[[1]], output_file_1)
      writeXStringSet(shuffled_reads[[2]], output_file_2)
    }
  }
}