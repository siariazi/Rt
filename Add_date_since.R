# adding date to fasta file generated from OOPidemic package using the meta data csv file
# this is good when in beast in set date, I choose "since some time in the past"
library(Biostrings)
directory <- "/Users/siavashriazi/Desktop/SFU/Rt/Rt codes/Chris/data/siavash-actualScens/"

# set working directory 
setwd(directory)

meta_file = "Scen02A_popsz10K_initSus15_wgs_metadata"
fasta_file = "Scen02A_popsz10K_initSus15_wgs"
output_fasta = paste(fasta_file,"3.fasta",sep = "")
meta = read.csv(paste(meta_file,".csv",sep = ""))
# note that this is working for natural numbers, if I have decimal number it should be changed
#meta$time = max(meta$collection_time)-meta$collection_time 
name_to_time <- setNames(meta$collection_time, meta$name)

# this is for age argument of LPhy, I don't know if it's going to be useful 
#time_str = paste(rev(meta$time),collapse =",")
#time_str2 = paste(rev(unique(meta$time)),collapse =",")

# Read the fasta file
fasta_seq <- readDNAStringSet(paste(fasta_file,".fasta",sep = ""))
# Modify the fasta headers
names(fasta_seq) <- sapply(names(fasta_seq), function(name) {
  paste(name, name_to_time[name], sep="|")
})

# Write the modified fasta to a new file
writeXStringSet(fasta_seq, output_fasta)
