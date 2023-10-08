# Load required libraries
library(polyester)
library(Biostrings)
library(rtracklayer)

# File paths
path_sequences <- "data/nucleotide_sequences.fna"
path_annotations <- "data/annotated_genome.gtf"

# Loading data
genome_sequence <- readDNAStringSet(path_sequences)
names(genome_sequence) <- gsub(" .*", "", names(genome_sequence))  # need to split off id # nolint

annotations <- read.table(path_annotations, comment.char="#", quote="",  # nolint
                         col.names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), # nolint
                         fill = TRUE, sep="\t") # nolint

exon_annotations <- annotations[annotations$feature == "exon", ]
fold_changes_vector <- c(rep(1, 500), rep(2, 50))

# Synthetic read generation
simulate_experiment(seqpath = genome_sequence,
                    gtf = annotations,
                    num_reps = 5, # Number of replicates
                    readlen = 10, # Read length
                    numreads = 50, # Number of reads
                    outdir = "read_outputs/",
                    strand_specific = FALSE,
                    fold_changes = fold_changes_vector)