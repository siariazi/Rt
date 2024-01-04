# sampling meta file with a defined limit and proportion to come up with a sampled meta and corrosponding fasta file
# the ouptput files should be later be oppened with Add_date_since.R
library(Biostrings)
library(readr)
library(dplyr)

directory <- "/Users/siavashriazi/SFU/Rt/Rt codes/Chris/data/siavash-actualScens/"
setwd(directory)

fasta_file_name <- "Scen02A_popsz10K_initSus15_wgs"
meta_data <- "Scen02A_popsz10K_initSus15_wgs_metadata"
# Read the FASTA file
fasta_file <- readDNAStringSet(paste(fasta_file_name,".fasta",sep = ""))

# Read the CSV file with sample IDs
meta_data_df <- read_csv(paste(meta_data,".csv",sep=""))

# Create a dataframe from meta data (all data) that shows the cumulative recurrence in each time 
meta_check <- meta_data_df %>%
  group_by(collection_time) %>%
  summarize(recurrence = n())

#write.csv(meta_check,paste(meta_data,'_all_recurrence.csv',sep = ""))

# proportion of sampling
prop = 0.2

# creating the meta data of sampled population 
sampled_meta <- meta_data_df %>% sample_frac(prop)

# writing the meta data of sampled recurrence in a file 
write.csv(sampled_meta,paste(meta_data,'_samp_prop_',prop,".csv",sep = ""))

# Create the new dataframe to calculate the summary of recurrence in sampled 
sampled_check <- sampled_meta %>%
  group_by(collection_time) %>%
  summarize(recurrence = n())

# writing the summary of sampled recurrence in a file 
write.csv(sampled_check,paste(meta_data,'_samp_recur_prop_',prop,".csv",sep = ""))

plot(recurrence~collection_time,meta_check,type='l',xlab="time")
points(recurrence~collection_time,sampled_check)
title(main = paste("sample prop:",prop))

# Extract the IDs from the 'name' column
ids_to_extract <- sampled_meta$name
sampled_fasta <- fasta_file[names(fasta_file) %in% ids_to_extract] 

# Write the filtered sequences to a new FASTA file
writeXStringSet(sampled_fasta, paste(fasta_file_name,"_prop_",prop,".fasta",sep = ""))

############### sampling a fix number
# Set the fixed number of samples for overrepresented times
lim <- 10 

# Identify time points with more than limit samples
overrepresented_times <- meta_data_df %>%
  group_by(collection_time) %>%
  summarize(count = n()) %>%
  filter(count > lim) %>%
  pull(collection_time)

# Sample a fixed number of rows for times that are overrepresented, and combine the results
sampled_meta2 <- meta_data_df %>%
  group_by(collection_time) %>%
  mutate(
    sampled = if_else(
      collection_time %in% overrepresented_times,
      row_number() %in% sample(row_number(), min(lim, n())),
      TRUE
    )
  ) %>%
  filter(sampled == TRUE) %>%
  select(-sampled)

# writing the meta data of sampled recurrence in a file 
write.csv(sampled_meta2,paste(meta_data,'_samp_limit',lim,".csv",sep = ""))

# Create the new dataframe to calculate the summary of recurrence in sampled 
sampled_check2 <- sampled_meta2 %>%
  group_by(collection_time) %>%
  summarize(recurrence = n())

# writing the summary of sampled recurrence in a file 
write.csv(sampled_check2,paste(meta_data,'_samp_recur_limit',lim,".csv",sep = ""))

plot(recurrence~collection_time,meta_check,type='l',xlab="time")
points(recurrence~collection_time,sampled_check2)
title(main = paste("limit:",lim))


# Extract the IDs from the 'name' column
ids_to_extract2 <- sampled_meta2$name
sampled_fasta2 <- fasta_file[names(fasta_file) %in% ids_to_extract2] 

# Write the filtered sequences to a new FASTA file
writeXStringSet(sampled_fasta2, paste(fasta_file_name,"_limit_",lim,".fasta",sep = ""))
####################
# Third scenario where there is two sampling proportions
samp_time <- 20 # this is where sampling changes from prop1 to prop2
prop1 = 0.2
prop2 = 0.05

# Create sampled_meta3 with two sampling proportions
sampled_meta3 <- meta_data_df %>%
  mutate(
    sampled = if_else(collection_time <= samp_time, row_number() %in% sample(row_number(), n() * prop1), row_number() %in% sample(row_number(), n() * prop2))
  ) %>%
  filter(sampled) %>%
  select(-sampled)

# Create the new dataframe to calculate the summary of recurrence in sampled 
sampled_check3 <- sampled_meta3 %>%
  group_by(collection_time) %>%
  summarize(recurrence = n())

plot(recurrence~collection_time,meta_check,type='l',xlab="time")
points(recurrence~collection_time,sampled_check3)
title(main = paste("cut-off time: ",samp_time,", ","rates: ",prop1,", ",prop2,sep = ""))

# Extract the IDs from the 'name' column
ids_to_extract3 <- sampled_meta3$name
sampled_fasta3 <- fasta_file[names(fasta_file) %in% ids_to_extract3] 
# writing the meta data of sampled recurrence in a file 
write.csv(sampled_meta3,paste(meta_data,'_time_',samp_time,"_rates_",prop1,"_",prop2,".csv",sep = ""))

# writing the summary of sampled recurrence in a file 
write.csv(sampled_check3,paste(meta_data,"_samp_recur_time_",samp_time,"_rates_",prop1,"_",prop2,".csv",sep = ""))

# Write the filtered sequences to a new FASTA file
writeXStringSet(sampled_fasta3, paste(fasta_file_name,'_time_',samp_time,"_rates_",prop1,"_",prop2,".fasta",sep = ""))
####################
# limit of counts, so for any number larger than that, I take a fraction of that time
lim = 5 # 10 is good with prop = 0.2

# Identify time points with more than limit samples
overrepresented_times <- meta_data_df %>%
  group_by(collection_time) %>%
  summarize(count = n()) %>%
  filter(count > lim) %>%
  pull(collection_time)

# Sample proportion (prop) for each overrepresented time point, and combine the results
sampled_meta <- meta_data_df %>%
  group_by(collection_time) %>%
  mutate(sampled = if_else(collection_time %in% overrepresented_times, row_number() %in% sample(row_number(), n() * prop), TRUE)) %>%
  filter(sampled == TRUE) %>%
  select(-sampled)

# writing the new sampled meta data
write.csv(sampled_meta,paste(meta_data,"_lim:",lim,"_prop:",prop,".csv",sep = ""))
