library(tidyverse)
library(dplyr)
library(Rsubread)
library(stringr)

setwd('~/aoa_mags/mag_mt_mapping/bowtie2_output_sensitive')

sorted_files <- list.files(pattern = 'sorted.bam$')

stats <- data.frame()
coverage <- data.frame()


for(i in sorted_files){
    
    fc <- featureCounts(files = i, 
                        #annot.ext = '~/aoa_mags/mag_annotation_mage/Candidatus_Nitrosocosmicus_sp._EC11.gff',
                        annot.ext = '~/aoa_mags/mag_annotation_mage/Nitrosotalea_sp._EC131.gff',
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType = 'CDS',
                        GTF.attrType = 'locus_tag')
    
    
    if(nrow(stats) == 0){
      stats <- fc$stat
    } else { stats <- left_join(stats, fc$stat, by = 'Status')}
    
    if(nrow(coverage) == 0){
      coverage <- fc$counts
    } else {
      coverage <- merge(coverage, fc$counts, by = 0)
      coverage <- coverage %>% column_to_rownames(var = 'Row.names')
    }
  
  
}


stats <- stats %>% column_to_rownames(var = 'Status')
stats <- stats %>% rbind(Total_Reads = colSums(stats))
stats <- as.data.frame(t(stats))
stats <- stats %>% mutate(Perc_Assigned = Assigned / Total_Reads)

write.csv(stats, '~/aoa_mags/featurecounts_mage/ec131_coverage_stats.csv', row.names = TRUE)
write.csv(coverage, '~/aoa_mags/featurecounts_mage/ec131_coverage_details.csv', row.names=TRUE)

library(ggplot2)
library(stringr)

stats$Sample <- sapply(str_split(row.names(stats), "[.]"), `[`, 2)

stats$Sample <- factor(stats$Sample, levels = c("Apr_PD1_archaea", "Apr_PD2_archaea", "Apr_PD2_co", "Apr_PD3_archaea",
                                                "Aug_PD1_archaea", "Aug_PD2_archaea", "Aug_PD3_archaea",
                                                "Oct_PD1_archaea", "Oct_PD2_archaea", "Oct_PD3_archaea"),
                       labels = c("Apr18_PD1", "Apr18_PD2", "Apr18_PD2**", "Apr18_PD3",
                                  "Aug18_PD1", "Aug18_PD2", "Aug18_PD3",
                                  "Oct18_PD1", "Oct18_PD2", "Oct18_PD3"))

ggplot(stats, aes(x = Sample, y = Perc_Assigned*100)) + geom_col() + xlab("") +
  ylab("% Coverage of MAG") +
  ggtitle("PD Metatranscriptome Coverage of Culture 13 MAG (Relative)") +
  theme_classic()

ggplot(stats, aes(x = Sample, y = Assigned)) + geom_col() + xlab("") +
  ylab("# Aligned Reads") +
  ggtitle("PD Metatranscriptome Coverage of Culture 13 MAG (Absolute)") +
  theme_classic()

coverage <- coverage[rowSums(coverage)>0,]

coverage <- coverage %>% mutate(Volume = rowSums(coverage)) %>% arrange(desc(Volume))

write.csv(coverage, '~/aoa_mags/featurecounts_mage/ec5_coverage.csv', row.names=TRUE)
