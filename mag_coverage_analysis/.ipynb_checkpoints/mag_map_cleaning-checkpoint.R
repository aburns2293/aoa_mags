
library(readr)
library(dplyr)
library(stringr)

# import aoa sequence ids

aoa_sequence_ids <- read_csv("aoa_sequence_ids.txt", 
                             col_names = FALSE)
aoa_sequence_ids$ID <- sapply(strsplit(aoa_sequence_ids$X1, " "), `[`, 1)
aoa_sequence_ids$Function <- str_replace(aoa_sequence_ids$X1, "\\bWP\\w*.1", '')
aoa_sequence_ids$Function <- str_replace(aoa_sequence_ids$Function, "MULTISPECIES: ", '')
aoa_sequence_ids <- rename(aoa_sequence_ids, Fasta_Heading = X1)
write.csv(aoa_sequence_ids, 'aoa_sequence_ids.csv', row.names = FALSE)

# read in files, save contig name, read and best match score

files <- list.files(pattern = 'fasta.txt')
contig_names <- c('11_contig_2', '13_21_contig_3', '13_21_contig_1', '5_21_3_contig_26')
time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')

top_hits <- data.frame(matrix(nrow = 0, ncol = 14))
names(top_hits) <- c('Transcript_ID', 'ID', 'Perc_ID', 'Length', 'Mismatch', 'GapOpen', 'QStart', 'QEnd', 'SStart', 'SEnd', 'EValue', 'Bitscore', 'Contig', 'Sample')

for(i in contig_names){
  for(k in time_points){
    for (j in files) {
      if (str_detect(j, i) == TRUE & str_detect(j, k) == TRUE){
        data <- read_delim(j, delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)
        names(data) <- c('Transcript_ID', 'ID', 'Perc_ID', 'Length', 'Mismatch', 'GapOpen', 'QStart', 'QEnd', 'SStart', 'SEnd', 'EValue', 'Bitscore')
        hits <- data %>% group_by(Transcript_ID) %>% arrange(desc(Perc_ID)) %>% filter(row_number() == 1)
        hits$Contig <- i
        hits$Sample <- k
        
        top_hits <- rbind(top_hits, hits)
      }
    }
  }
}

# bind top hits file to sequence IDs

top_hits <- left_join(top_hits, aoa_sequence_ids, by = 'ID')

non_hypothetical <- top_hits %>% filter(str_detect(Function, 'hypothetical protein') == FALSE)

non_hypothetical <- non_hypothetical %>% select(Transcript_ID, ID, Perc_ID, Length, Mismatch, EValue, Bitscore, Function, Contig, Sample)

write.csv(non_hypothetical, 'metatranscriptome_protein_matches.csv', row.names = FALSE)


# investigating quality of hits

matches <- read.csv('metatranscriptome_protein_matches.csv')

hist(matches$Bitscore)
