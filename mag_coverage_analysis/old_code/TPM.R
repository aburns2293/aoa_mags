
match <- read.csv('metatranscriptome_protein_matches.csv')
match$description <- sapply(strsplit(match$description, "\\["), `[`, 1)
match.sums <- match %>% group_by(Contig, Sample, description, ID, sequence_length) %>% summarise(Read_Count = n())

match.sums <- match.sums %>% mutate(ID_Length_Kb = sequence_length/1000,
                                    Reads_per_Kb = Read_Count/ID_Length_Kb)

contig_names <- c('11_contig_2', '13_21_contig_3', '13_21_contig_1', '5_21_3_contig_26')
time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')

match.sums$Per_Million_Factor <- NA

for(i in contig_names){
  for(j in time_points){
    subset <- match.sums[match.sums$Contig == i & match.sums$Sample == j,]
    total_reads_per_kb <- sum(subset$Reads_per_Kb)
    reads_per_kb_million <- total_reads_per_kb / 1000000
    match.sums[match.sums$Contig == i & match.sums$Sample == j,]$Per_Million_Factor <- reads_per_kb_million
  }
}

match.sums <- match.sums %>% mutate(TPM = Reads_per_Kb/Per_Million_Factor)

write.csv(match.sums, 'transcript_hits_tpm.csv', row.names = FALSE)

# heat map of ammonia oxidation functions 

Functions <- c(' urease subunit gamma ', ' urease subunit beta ', ' urease subunit alpha ', ' urease accessory protein UreG ', ' multicopper oxidase family protein ',
               ' multicopper oxidase domain-containing protein ', ' ammonium transporter ', ' ammonium transporter, partial ', ' methane monooxygenase/ammonia monooxygenase subunit C ',
               ' methane monooxygenase/ammonia monooxygenase subunit B ', ' ammonia monooxygenase ', ' methane monooxygenase/ammonia monooxygenase subunit C, partial ')

ao <- match.sums[match.sums$description %in% Functions,]

ao$description <- factor(ao$description, levels = c(' urease subunit gamma ', ' urease subunit beta ', ' urease subunit alpha ', ' urease accessory protein UreG ', ' multicopper oxidase family protein ',
                                              ' multicopper oxidase domain-containing protein ', ' ammonium transporter, partial ', ' ammonium transporter ', ' methane monooxygenase/ammonia monooxygenase subunit C ', ' methane monooxygenase/ammonia monooxygenase subunit C, partial ',
                                              ' methane monooxygenase/ammonia monooxygenase subunit B ', ' ammonia monooxygenase '))

ao.sums <- ao %>% group_by(Sample, Contig, description) %>% summarise(total_tpm = sum(TPM))

ggplot(ao.sums, aes(x = Sample, y = description, fill = total_tpm)) + 
  geom_tile() + 
  facet_wrap(vars(Contig)) + 
  xlab('') + 
  ylab('') + 
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = 'TPM') + 
  geom_vline(xintercept = 3.5) + 
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(labels = c('', 'April', '', '', 'August', '', '', 'October', ''), expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme_classic() + 
  ggtitle('Expression of mapped reads: Ammonia oxidation') + 
  theme(panel.spacing = unit(0, 'lines'))

## frequently expressed reads

freq <- match.sums %>% group_by(Sample, Contig, description) %>% summarise(total_tpm = sum(TPM))

freq.count <- freq %>% group_by(description) %>% summarise(count = n())

freq.count.short <- freq.count %>% filter(count >= 27 )
