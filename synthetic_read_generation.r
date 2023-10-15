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

exon_annotations <- annotations[annotations$feature == "exon", ] # 35384
fold_changes_vector <- c(rep(1, 1), rep(1, 1))

# Synthetic read generation
simulate_experiment(seqpath = genome_sequence,
                    gtf = exon_annotations,
                    num_reps = 1, # sets of reads
                    reads_per_transcript = 10, # reads per set
                    readlen = 50, # chars per read
                    fraglen = 100, # fragment length mean
                    fragsd = 10, # fragment length st. dev.
                    error_rate = 0.0, # % incorrect
                    fold_changes = fold_changes_vector, # diff. exp. controller
                    paired = TRUE, # both sides of helix
                    size = NULL, # size of reads generated w/ -binomial dist.
                    write_info = FALSE, # write diff. exp. transcripts
                    transcriptid = NULL, # ids for transcripts
                    seed = NULL, # rng seed
                    outdir = "polyester_reads/")


# multiple experiment setup
# read_lengths <- c(50, 75, 100)
# error_rates <- c(0.0, 0.01, 0.02)
# for (readlen in read_lengths) {
#     for (error in error_rates) {
#         simulate_experiment(seqpath = genome_sequence,
#                             gtf = exon_annotations,
#                             num_reps = 1,
#                             reads_per_transcript = 10,
#                             readlen = readlen,
#                             fraglen = 100,
#                             fragsd = 10,
#                             error_rate = error,
#                             fold_changes = fold_changes_vector,
#                             paired = TRUE,
#                             outdir = paste0("polyester_reads/readlen_", readlen, "_error_", error, "/"))
#     }
# }


# seed arg helpful for consistency later
