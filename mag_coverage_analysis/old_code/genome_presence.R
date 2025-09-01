
key <- read_delim("~/aoa_mags/mag_annotation/ec11/EC11.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE) %>% select(-gene, -EC_number, -COG)

key$product <- sapply(strsplit(key$product, "\\["), `[`, 1)
key$product <- str_replace(key$product, 'MULTISPECIES: ', '')

key.sums <- key %>% group_by(product) %>% summarise(count = n())

key %>% group_by(ftype) %>% summarise(count = n())

key[key$ftype == 'rRNA',]

## coverage of genomes

library(dplyr)
library(ggplot2)
library(stringr)

files <- c('ec5_TPM.csv', 'ec11_TPM.csv', 'ec131_TPM.csv', 'ec133_TPM.csv')
total_tpm <- data.frame()

for(i in files){
  file_path <- paste0('mag_coverage_analysis/', i)
  temp <- read.csv(file_path) %>% filter(TPM > 0)
  total_tpm <- rbind(total_tpm, temp)
}

counts <- total_tpm %>% group_by(Culture, Sample) %>% summarise(Count = n())
counts$Sample <- factor(counts$Sample, levels = c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3'))
counts$Month <- str_sub(counts$Sample, 1, 3)

ggplot(counts, aes(x = Sample, y = Count, fill = Month)) + geom_col() + facet_wrap(vars(Culture)) + 
  ylab('# CDS Transcribed') + 
  
  theme_classic() + 
  theme(axis.text.x = element_blank())

ggplot(counts, aes(x = Culture, y = Count)) + geom_boxplot()
