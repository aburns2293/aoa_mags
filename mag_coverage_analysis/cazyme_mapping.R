
setwd('~/aoa_mags/mag_mt_mapping/CAZymes/output')

files <- list.files()
contigs <- c('contig5', 'contig11', 'contig131', 'contig133')
time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')

top_hits <- data.frame(matrix(nrow = 0, ncol = 14))
names(top_hits) <- c('Transcript_ID', 'ID', 'Perc_ID', 'Length', 'Mismatch', 'GapOpen', 'QStart', 'QEnd', 'SStart', 'SEnd', 'EValue', 'Bitscore', 'Contig', 'Sample')

for(i in contigs){
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

# sort contig IDs
setwd('..')
out <- read_delim('133_cazyme_CGC/hmmer.out', delim = '\t')
out.short <- out %>% select(`HMM Profile`, `Profile Length`, `Gene ID`)
out.short$Contig <- 'contig133'

contig_ids_full <- rbind(contig_ids_full, out.short)
write.csv(contig_ids_full, 'dbcan_ids.csv', row.names = FALSE)

contig_ids_full <- rename(contig_ids_full, ID = `Gene ID`)
contig_ids_full <- contig_ids_full %>% select(-Contig)

# bind

top_hits_full <- left_join(top_hits, contig_ids_full, by = 'ID', relationship = 'many-to-many')

write.csv(top_hits_full, 'dbcan_output_detailed.csv', row.names = FALSE)

# tpm 

caz.sums <- top_hits_full %>% group_by(Contig, Sample, `HMM Profile`, `Profile Length`) %>% summarise(Read_Count = n())

caz.sums <- caz.sums %>% mutate(Profile_Length_Kb = `Profile Length`/1000,
                                    Reads_per_Kb = Read_Count/Profile_Length_Kb)

caz.sums$Per_Million_Factor <- NA

for(i in contigs){
  for(j in time_points){
    subset <- caz.sums[caz.sums$Contig == i & caz.sums$Sample == j,]
    total_reads_per_kb <- sum(subset$Reads_per_Kb)
    reads_per_kb_million <- total_reads_per_kb / 1000000
    caz.sums[caz.sums$Contig == i & caz.sums$Sample == j,]$Per_Million_Factor <- reads_per_kb_million
  }
}

caz.sums <- caz.sums %>% mutate(TPM = Reads_per_Kb/Per_Million_Factor)

write.csv(caz.sums, 'cazyme_hits_tpm.csv', row.names = FALSE)

## figure

caz.sums <- read.csv('cazyme_hits_tpm.csv')

caz.sums$`HMM.Profile` <- str_sub(caz.sums$`HMM.Profile`, end = -5)

ggplot(caz.sums, aes(x = Sample, y = `HMM.Profile`, fill = Read_Count)) + 
  geom_tile() + 
  facet_wrap(vars(Contig)) + 
  xlab('') + 
  ylab('') + 
  scale_fill_distiller(palette = 'BuGn', direction = 1, name = 'Reads') + 
  geom_vline(xintercept = 3.5) + 
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(labels = c('', 'April', '', '', 'August', '', '', 'October', ''), expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme_classic() + 
  ggtitle('CAZyme Expression in PD Metatranscriptomes') + 
  theme(panel.spacing = unit(0, 'lines'))
  