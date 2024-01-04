# using id's of sampling files to make a fast file that corrosponds to the meta data 
library(Biostrings)
library(readr)
directory <- "/Users/siavashriazi/Desktop/SFU/Rt/Rt codes/Chris/data/siavash-actualScens/"
setwd(directory)

fasta_file_name <- "Scen02A_popsz10K_initSus15_wgs_full.fasta"
sample_file_name <- "scen02B_keep75tokeep25percentIDS.csv"

# Read the FASTA file
fasta_file <- readDNAStringSet(fasta_file_name)

# Read the CSV file with sample IDs
sample_ids_df <- read_csv(sample_file_name)

# Extract sample IDs and convert them to numeric
sample_ids <- as.numeric(sample_ids_df$x)

# Function to extract ID from FASTA header
extract_id <- function(header) {
  as.numeric(strsplit(substring(header, 4), "_")[[1]][1])
}

# Filter the sequences based on IDs
filtered_sequences <- fasta_file[sapply(names(fasta_file), extract_id) %in% sample_ids]

# Write the filtered sequences to a new FASTA file
writeXStringSet(filtered_sequences, "filtered_sequences.fasta")
