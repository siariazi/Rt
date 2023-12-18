# adding date to fasta file generated from OOPidemic package using the meta data csv file
# with this version in beast you need to choose "before the present" in dates
library(Biostrings)
directory <- "/Users/siavashriazi/Desktop/SFU/Rt/Rt codes/Chris/data"

# set working directory 
setwd(directory)

meta_file = "seir_v1_wgs_metadata"
fasta_file = "seir_v1_wgs"
output_fasta = paste(fasta_file,"2.fasta",sep = "")
meta = read.csv(paste(meta_file,".csv",sep = ""))
# note that this is working for natural numbers, if I have decimal number it should be changed
meta$time = max(meta$collection_time)-meta$collection_time 
name_to_time <- setNames(meta$time, meta$name)

# this is for age argument of LPhy, I don't know if it's going to be useful 
time_str = paste(rev(meta$time),collapse =",")
time_str2 = paste(rev(unique(meta$time)),collapse =",")

# Read the fasta file
fasta_seq <- readDNAStringSet(paste(fasta_file,".fasta",sep = ""))
# Modify the fasta headers
names(fasta_seq) <- sapply(names(fasta_seq), function(name) {
  paste(name, name_to_time[name], sep="|")
})

# Write the modified fasta to a new file
writeXStringSet(fasta_seq, output_fasta)
