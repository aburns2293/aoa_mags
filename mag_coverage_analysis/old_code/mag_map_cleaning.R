
library(readr)
library(dplyr)
library(stringr)

# import aoa sequence ids

setwd('new_blast_mapped_reads')

aoa_sequence_ids <- read_csv("nitrososphaerota_sequence_info.csv") %>% select(-`...1`)

aoa_sequence_ids$description <- str_replace(aoa_sequence_ids$description, "\\bWP\\w*.1", '')
aoa_sequence_ids$description <- str_replace(aoa_sequence_ids$description, "MULTISPECIES: ", '')
write.csv(aoa_sequence_ids, 'nitro_sequence_info.csv', row.names = FALSE)

# read in files, save contig name, read and best match score

files <- list.files(pattern = 'mapped.txt')
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
        hits <- data %>% group_by(Transcript_ID) %>% arrange(desc(Bitscore)) %>% filter(row_number() == 1)
        hits$Contig <- i
        hits$Sample <- k
        
        top_hits <- rbind(top_hits, hits)
      }
    }
  }
}

# bind top hits file to sequence IDs

top_hits <- left_join(top_hits, aoa_sequence_ids, by = join_by('ID' == 'name'))

top_hits_short <- top_hits %>% select(Transcript_ID, ID, Perc_ID, Mismatch, EValue, Bitscore, Contig, Sample, description, sequence_length)

write.csv(top_hits_short, 'metatranscriptome_protein_matches.csv', row.names = FALSE)


## trying with full genome assignment

library(readr)
library(dplyr)
library(stringr)

blastx <- read_delim("mag_annotation/blastx_output/ec5.blastx.txt", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)

names(blastx) <- c('ORF_ID', 'ID', 'Perc_ID', 'Length', 'Mismatch', 'GapOpen', 'QStart', 'QEnd', 'SStart', 'SEnd', 'EValue', 'Bitscore')

seq_info <- read_csv('mag_annotation/blastx_output/nitrososphaerota_sequence_info.csv') %>% select(-`...1`)

names(seq_info) <- c('ID', 'BLAST_Description', 'Blast_Query_Length')

seq_info$BLAST_Description <- str_replace(seq_info$BLAST_Description, "\\bWP\\w*.1\\s", '')
seq_info$BLAST_Description <- str_sub(seq_info$BLAST_Description, end = -2)
seq_info$BLAST_Description <- sapply(strsplit(seq_info$BLAST_Description, "\\["), `[`, 1)

kegg <- read_delim("mag_annotation/blastx_output/ec5.kegg.txt", 
                   delim = "\t", escape_double = FALSE, 
                   col_names = FALSE, trim_ws = TRUE) %>% select(-X4, -X5, -X6)

names(kegg) <- c('ORF_ID', 'KEGG_Pathway', 'KEGG_Description')

blastx <- left_join(blastx, kegg, by = 'ORF_ID')
blastx <- left_join(blastx, seq_info, by = 'ID')

investigate_further <- c()

hits <- top_hits <- data.frame(matrix(nrow = 0, ncol = 16))
names(hits) <- names(blastx)

## save all hits that have a consensus of best: percent identity, bit score, e value and lowest mismatch

## otherwise, save identity names to investigate further

for(i in blastx$ORF_ID %>% unique()){
  transcript_matches <- blastx %>% filter(ORF_ID == i)
  best.bitscore <- transcript_matches %>% arrange(desc(Bitscore)) %>% filter(row_number() == 1)
  best.evalue <- transcript_matches %>% arrange(desc(EValue)) %>% filter(row_number() == 1)
  best.mismatch <- transcript_matches %>% arrange(Mismatch) %>% filter(row_number() == 1)
  best.percid <- transcript_matches %>% arrange(desc(Perc_ID)) %>% filter(row_number() == 1)
  
  matches <- c(best.bitscore$ID, best.evalue$ID, best.mismatch$ID, best.percid$ID)
  if (length(matches %>% unique()) == 1){
    hits <- rbind(hits, best.bitscore)
  } else{
    investigate_further <- c(investigate_further, i)
  }
    
}

## investigation of scores for unambiguous hits:
## culture 5: 1842 proteins (57 ambiguous)

poor_classification <- c()

## low bitscore (< 50) - visually inspect to see percent ID, if it has a KEGG hit
## if the low bitscore protein is already labeled 'hypothetical', leave it

summary(hits$Bitscore)

low.bitscore <- hits %>% filter(Bitscore < 50,
                                BLAST_Description != 'hypothetical protein')

## low percent ID (< 40%) for non-hypothetical proteins

summary(hits$Perc_ID)
low.ID <- hits %>% filter(Perc_ID < 40,
                          BLAST_Description != 'hypothetical protein')

poor_classification <- c(poor_classification, 'contig_26_1872', 'contig_26_625', 'contig_26_458')

## over 50% mismatch

hits <- hits %>% mutate(Perc_Mismatch = Mismatch/Length)
high.mismatch <- hits %>% filter(Perc_Mismatch > 0.5,
                                 BLAST_Description != 'hypothetical protein')

poor_classification <- c(poor_classification, 'contig_26_1843')

## investigation of ambiguous hits

ambig <- blastx[blastx$ORF_ID %in% investigate_further,]

ambig_cleaned <- data.frame(matrix(nrow = 0, ncol = 17))
names(ambig_cleaned) <- names(hits)


ambig_cleaned <- ambig[ambig$ORF_ID == 'contig_26_1089' & ambig$ID == 'WP_157927101.1',]

ambig_cleaned <- rbind(ambig_cleaned, ambig[ambig$ORF_ID %in% c('contig_26_112', 'contig_26_1126', 'contig_26_1130', 'contig_26_1132', 'contig_26_114', 'contig_26_1159', 'contig_26_1163', 'contig_26_1212', 
                                                                'contig_26_1332'),])

ambig_cleaned <- rbind(ambig_cleaned, ambig[ambig$ORF_ID == 'contig_26_1374' & ambig$ID == 'WP_157926875.1',])



## attemping manual analysis - keep hypotheticals, remove if percent id < 30% and bit score < 50,
## then pick top ten entries per protein by bitscore

blastx <- blastx %>% filter(BLAST_Description != 'hypothetical protein' & Perc_ID > 30 & Bitscore > 50) %>% group_by(ORF_ID) %>% top_n(10, Bitscore)
  #arrange(desc(Bitscore), .by_group = TRUE) #%>% filter(row_number %in% seq(1, 10, 1)) 

write.csv(blastx, 'manual_attempt.csv', row.names = FALSE)
