
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(MetBrewer)
library(cowplot)

setwd('~/aoa_mags/mag_coverage_analysis')

exp <- read.csv('transcript_hits_tpm.csv')

#exp$product <- URLdecode(exp$product)

#write.csv(exp, 'transcript_hits_tpm.csv', row.names = FALSE)

##### merge with annotated list

expressed_genes_list <- read.csv("TableS2_DroughtResponseExpression.csv", skip = 1)

exp <- exp %>% filter(Total > 0)

eg_annotated <- left_join(exp, expressed_genes_list, by = join_by('product' == 'Name'))

setdiff(exp$product, expressed_genes_list$Name)
setdiff(expressed_genes_list$Name, exp$product)

eg_annotated$Date <- str_sub(eg_annotated$Sample, 1, 3)

write.csv(eg_annotated, 'annotated_transcript_hits_tpm.csv', row.names = FALSE)

## evaluating different groups

ggplot(eg_annotated, aes(x = Date, y = TPM, col = Culture)) + geom_point() + facet_wrap(vars(`Putative.Pathway`), scales = 'free') +
  theme_classic()

eg_annotated %>% select(locus_tag, Culture, `Putative.Pathway`) %>% unique() %>%
  group_by(Culture) %>% summarise(total=n())

eg_annotated %>% select(locus_tag, Culture, `Putative.Pathway`) %>% unique() %>%
  group_by(Culture, `Putative.Pathway`) %>% summarise(total=n())
