
library(ggplot2)

setwd('~/aoa_mags/mag_coverage_analysis')

ec131 <- read.csv('ec131_mage_coverage.csv')
ec11 <- read.csv('ec11_mage_coverage.csv')

coverage <- rbind(ec11, ec131)

# tpm

time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
cultures <- c('Nitrosocosmicus', 'Nitrosotalea EC131')

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

data.long$product <- URLdecode(data.long$product)

write.csv(data.long, 'transcript_hits_tpm.csv', row.names = FALSE)

#### evaluating most highly transcribed genes

transcripts <- read.csv("transcript_hits_tpm.csv")

transcripts.avg <- transcripts %>% filter(Culture %in% c("Nitrosotalea EC131", "Nitrosocosmicus")) %>%
  group_by(Culture, locus_tag, product) %>%
  summarise(mean = mean(TPM), sd = sd(TPM))

transcripts.avg.nc <- transcripts.avg %>% filter(Culture == "Nitrosocosmicus") %>% arrange(desc(mean))

transcripts.avg.nt <- transcripts.avg %>% filter(Culture == "Nitrosotalea EC131") %>% arrange(desc(mean))
