# ~~~ Setup ~~~ #
# Load required libraries
library(polyester)
library(Biostrings)
library(rtracklayer)
library(ShortRead)

set.seed(15) # for reproducibility

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
error_rates <- list(low=0.001, medium=0.005, high=0.01)

# Fold changes for genes
# use random fold changes between 0.5 (downregulation) and 2 (upregulation)
# for the treatment condition compared to control.
fold_changes <- runif(n = length(sample_exons$start), min = 0.5, max = 2)

# Define output directory
output_dir <- "polyester_reads"

# Simulate reads for control and treatment under different error scenarios
for (condition in c("control", "treatment")) {
  for (error_level in names(error_rates)) {
    simulate_experiment(
      seqpath = genome_sequence,
      gtf = sample_exons,
      fold_changes = fold_changes,
      num_reps = 3,
      reads_per_transcript = 50,
      seq_len = 100,
      error_rate = error_rates[[error_level]],
      outdir = sprintf("%s/%s_%s", output_dir, condition, error_level)
      )
  }
}

# ~~~ Read Shuffle ~~~ #
shuffle_reads_fasta <- function(file_path, shuffle_percentage) {
  # Read the FASTA file
  reads <- readDNAStringSet(file_path)
  
  # Determine the number of reads to shuffle
  num_to_shuffle <- round(length(reads) * shuffle_percentage)
  
  # Randomly select reads to shuffle
  indices_to_shuffle <- sample(1:length(reads), num_to_shuffle)
  
  # Shuffle the selected reads
  shuffled_reads <- sample(reads[indices_to_shuffle])
  
  # Replace the original reads with the shuffled ones
  reads[indices_to_shuffle] <- shuffled_reads
  
  # Return the shuffled reads
  return(reads)
}

# Define shuffling levels and their respective percentages
shuffle_levels <- list(low = 0.05, medium = 0.25, high = 0.5)

# Loop through conditions, error scenarios, and replicates
for (condition in c("control", "treatment")) {
  for (error_level in names(error_rates)) {
    for (rep in 1:3) {
      for (end_num in 1:2) { # Loop over the paired-end reads
        # Define input file path
        input_file <- sprintf("%s/%s_%s/sample_%02d_%d.fasta", output_dir, condition, error_level, rep, end_num)
        
        # Shuffle reads based on the shuffling level
        shuffled_percentage <- shuffle_levels[[error_level]]
        shuffled_reads <- shuffle_reads_fasta(input_file, shuffled_percentage)
        
        # Write shuffled reads back to FASTA
        output_file <- sprintf("%s/%s_%s/shuffled_sample_%02d_%d.fasta", output_dir, condition, error_level, rep, end_num)
        writeXStringSet(shuffled_reads, output_file)
      }
    }
  }
}