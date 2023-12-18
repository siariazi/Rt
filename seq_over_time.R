library(Biostrings)
library(ggplot2)
library(lubridate)

# Read the fasta file
fasta_seq <- readDNAStringSet("/Users/siavashriazi/Desktop/SFU/BEAST files/Rt Beast/caroline2.fasta")

# Extract the dates from the names of the fasta sequences
dates <- sapply(names(fasta_seq), function(x) {
  unlist(strsplit(x, "\\|"))[3]
})

# Count the number of sequences per date
date_counts <- table(dates)

# Convert to a dataframe
df <- data.frame(date = as.Date(names(date_counts)), count = as.numeric(date_counts))

# Plot
ggplot(df, aes(x = date, y = count)) +
  geom_point() +
  labs(title = "Number of sequences over time", x = "Date", y = "Number of sequences") + theme_bw()


################
# Bin by week
# Add a week starting date for binning
df$week_start <- floor_date(df$date, unit = "week")

# Count the number of sequences per week
date_counts <- table(df$week_start)

# Convert to a dataframe
df_weekly <- data.frame(date = as.Date(names(date_counts)), count = as.numeric(date_counts))

# Plot
ggplot(df_weekly, aes(x = date, y = count)) +
  geom_point() +
  geom_line(group = 1) + theme_bw() +
  labs(title = "Number of sequences over time (binned by week)", x = "Week starting", y = "Number of sequences")

