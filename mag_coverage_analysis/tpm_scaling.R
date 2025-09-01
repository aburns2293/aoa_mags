
library(ggplot2)

setwd('~/aoa_mags/mag_coverage_analysis')

ec133 <- read.csv('ec133_mage_coverage.csv')
ec131 <- read.csv('ec131_mage_coverage.csv')
ec11 <- read.csv('ec11_mage_coverage.csv')
ec5 <- read.csv('ec5_mage_coverage.csv')

coverage <- rbind(ec5, ec11, ec131, ec133)

# tpm

time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
cultures <- c('Ca. Nitrosocosmicus', 'Nitrosotalea 5', 'Nitrosotalea 131', 'Nitrosotalea 133')

data.long <- pivot_longer(coverage, cols = -c('locus_tag', 'length', 'product', 'Total', 'Culture'), names_to = 'Sample', values_to = 'Reads')

data.long <- data.long %>% mutate(length_kb = length/1000,
                                  reads_per_kb_length = Reads/length_kb)

data.long$per_million_factor <- NA

for(i in cultures){
  for(j in time_points){
    subset <- data.long[data.long$Culture == i & data.long$Sample == j,]
    total_reads_per_kb <- sum(subset$reads_per_kb_length)
    reads_per_kb_million <- total_reads_per_kb / 1000000
    data.long[data.long$Culture == i & data.long$Sample == j,]$per_million_factor <- reads_per_kb_million
  }
}


data.long <- data.long %>% mutate(TPM = reads_per_kb_length/per_million_factor)

# scale by number of reads mapped - not necessary because TPM is already scaled

data.long$mrna_reads <- NA

ec131_cov <- read.csv('~/aoa_mags/featurecounts_mage/ec133_coverage_stats.csv') %>% select(Assigned)
ec131_cov$Sample <- time_points
for(i in time_points){
  val <- ec131_cov[ec131_cov$Sample == i,]$Assigned
  data.long[data.long$Culture == 'Nitrosotalea 133' & data.long$Sample == i,]$mrna_reads <- val
}

data.long <- data.long %>% mutate(TPM_Scaled_ReadsAssigned = TPM/mrna_reads)

data.long.pos <- data.long %>% filter(TPM_Scaled_ReadsAssigned > 0)

write.csv(data.long.pos, 'transcript_hits_tpm_scaled.csv', row.names = FALSE)
write.csv(data.long, 'transcript_hits_tpm_scaled_full.csv', row.names = FALSE)


#### removing scaling data

transcripts_hits <- transcript_hits_tpm_scaled_full %>% dplyr::select(-mrna_reads, -TPM_Scaled_ReadsAssigned, -GenomeLength, -mrna_reads_per_genome_length, -tpm_scaled_reads_genome)

write.csv(transcripts_hits, 'transcript_hits_tpm.csv', row.names = FALSE)
