library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(stringr)
# process minimap2 and prokka outputs

coverage <- read.csv('~/aoa_mags/featurecounts_mage/ec131_coverage_details.csv')
names(coverage) <- c('locus_tag', 'Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')

key <- read_delim("~/aoa_mags/mag_annotation_mage/Nitrosotalea_sp._EC131.gff", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 3)
key <- key %>% mutate(length = X5-X4)
key.id <- key %>% filter(X3 == 'CDS') %>% select(X9, length)
key.id$locus_tag <- str_extract(key.id$X9, 'locus_tag=[^;]+')
key.id$locus_tag <- str_sub(key.id$locus_tag, start = 11)
key.id$product <- str_extract(key.id$X9, 'product=[^;]+')
key.id$product <- str_sub(key.id$product, start = 9)
key.id <- key.id %>% select(-X9)

coverage_full <- left_join(coverage, key.id, by = 'locus_tag')
coverage_full <- coverage_full %>% mutate(Total = rowSums(coverage_full[,2:10]))
coverage_full$Culture <- 'Nitrosotalea EC131'

write.csv(coverage_full, '~/aoa_mags/mag_coverage_analysis/ec131_mage_coverage.csv', row.names = FALSE)

coverage_long <- pivot_longer(coverage_full, cols = c(-locus_tag, -product, -Total, -Culture, -length), names_to = 'Sample', values_to = 'Reads')

coverage_long <- coverage_long %>% arrange(desc(Total))

ggplot(coverage_long[coverage_long$Total > 10,], aes(x = Sample, y = product)) + geom_tile(aes(fill = Reads)) + 
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = 'Reads', na.value = 'grey90', limits = c(0.1, 5000)) + 
  geom_vline(xintercept = 3.5) + 
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(labels = c('', 'April', '', '', 'August', '', '', 'October', ''), expand = c(0,0)) + 
  xlab('') + 
  ylab('')



# coverage_long$product <- sapply(strsplit(coverage_long$product, "\\["), `[`, 1)
# coverage_long$product <- str_replace(coverage_long$product, 'MULTISPECIES: ', '')
# coverage_long$product <- str_sub(coverage_long$product, end = -2)
# coverage_long[coverage_long$product == 'hypothetical protei',]$product <- 'hypothetical protein'
# 
# coverage_long <- coverage_long %>% mutate(ID_Length_Kb = length_bp/1000,
#                                     Reads_per_Kb = Reads/ID_Length_Kb)
# 
# time_points <- c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
# 
# coverage_long$Per_Million_Factor <- NA
# 
# for(j in time_points){
#   subset <- coverage_long[coverage_long$Sample == j,]
#   total_reads_per_kb <- sum(subset$Reads_per_Kb)
#   reads_per_kb_million <- total_reads_per_kb / 1000000
#   coverage_long[coverage_long$Sample == j,]$Per_Million_Factor <- reads_per_kb_million
# }
# 
# 
# coverage_long <- coverage_long %>% mutate(TPM = Reads_per_Kb/Per_Million_Factor,
#                                           Culture = 'Cul133')
# 
# write.csv(coverage_long, 'mag_coverage_analysis/ec133_TPM.csv', row.names = FALSE)
# 
# # number of transcripts mapped to each 
# 
# coverage_files <- list.files('featureCounts', pattern = 'coverage_details.csv')
# 
# transcript_coverage <- data.frame()
# #names(transcript_coverage) <-  c('Culture', 'Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
# 
# for(i in coverage_files){
#   cov <- read.csv(paste0('featureCounts/',i))
#   names(cov) <- c('locus_tag', 'Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
#   cov <- cov %>% remove_rownames() %>% column_to_rownames(var = 'locus_tag')
#   total_cov <- colSums(cov)
#   total_cov <- c(i, total_cov)
#   transcript_coverage <- rbind(transcript_coverage, total_cov)
#   
# }
# 
# names(transcript_coverage) <-  c('Culture', 'Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3')
# 
# transcript_coverage <- transcript_coverage %>% remove_rownames() %>% column_to_rownames(var = 'Culture')
# 
# transcript.cov.m <- as.matrix.data.frame(transcript_coverage)
# class(transcript.cov.m) <- "numeric"
# 
# pheatmap(transcript.cov.m, cluster_rows = FALSE, cluster_cols = FALSE)
# mean(transcript.cov.m)
