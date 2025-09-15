
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(readxl)

setwd('~/aoa_mags/mag_mt_mapping/CAZymes')

tpm <- read_csv('~/aoa_mags/mag_coverage_analysis/transcript_hits_tpm.csv') %>% filter(Culture %in% c('Ca. Nitrosocosmicus', 'Nitrosotalea 131'))

## Nitrosocosmicus (NC11)

nc11_hmm <- read_delim("NC11_cazyme/dbCAN_hmm_results.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

nc11_hmm <- nc11_hmm %>% filter(`i-Evalue` < 1e-15,
                          Coverage > 0.8)

########################### run in blast conda env #################################

#makeblastdb -in ../NC11_contig_2.fasta -dbtype nucl
#tblastn -query NC11_cazyme_hmm_positive_hits.fasta -db ../NC11_contig_2.fasta -outfmt 6 -out NC11_cazyme_genome_loc.txt

###################################################################################

nc11_blastlocs <- read.delim('NC11_cazyme/NC11_cazyme_genome_loc.txt', header = FALSE)

colnames(nc11_blastlocs) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
nc11_blastlocs <- nc11_blastlocs %>% filter(pident >= 99)

########################### in MaGe ###############################################
### manually locate locus tags based on sequence starts and ends in blast output ##
###################################################################################

nc11_locus_tags <- read_xlsx('NC11_cazyme/NC11_Cazyme_IDs.xlsx')
nc11_hmm_locus_tags <- left_join(nc11_hmm, nc11_locus_tags, by = 'Target Name')

nc11_tpm <- left_join(nc11_hmm_locus_tags, tpm, by = 'locus_tag', relationship = 'many-to-many')


ggplot(nc11_tpm, aes(y = product, x = Sample, fill = TPM, col = `HMM Name`)) +
  geom_tile() 

## Nitrosotalea (NT131)

nt131_hmm <- read_delim("NT131_cazyme/dbCAN_hmm_results.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

nt131_hmm <- nt131_hmm %>% filter(`i-Evalue` < 1e-15,
                                Coverage > 0.8)

########################### run in blast conda env #################################

#makeblastdb -in ../13_21_contig_1.fasta -dbtype nucl
#tblastn -query NT131_cazyme_hmm_positive_hits.fasta -db ../13_21_contig_1.fasta -outfmt 6 -out NT131_cazyme_genome_loc.txt

###################################################################################

nt131_blastlocs <- read.delim('NT131_cazyme/NT131_cazyme_genome_loc.txt', header = FALSE)

colnames(nt131_blastlocs) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
nt131_blastlocs <- nt131_blastlocs %>% filter(pident >= 99)

########################### in MaGe ###############################################
### manually locate locus tags based on sequence starts and ends in blast output ##
###################################################################################

nt131_locus_tags <- read_xlsx('NT131_cazyme/NT131_Cazyme_IDs.xlsx')
nt131_hmm_locus_tags <- left_join(nt131_hmm, nt131_locus_tags, by = 'Target Name')

nt131_tpm <- left_join(nt131_hmm_locus_tags, tpm, by = 'locus_tag', relationship = 'many-to-many')

# Join data
 
cazyme_tpm <- rbind(nc11_tpm, nt131_tpm)
 
write.csv(cazyme_tpm, '~/aoa_mags/mag_coverage_analysis/cazyme_tpm_coverage.csv', row.names = FALSE)


######################################################################################################

setwd('~/aoa_mags/mag_coverage_analysis/')
cazyme_tpm <- read_csv('cazyme_tpm_coverage.csv')

ggplot(cazyme_tpm, aes(y = `HMM Name`, x = Sample, fill = TPM)) +
  geom_tile() +
  facet_wrap(vars(Culture), scales = 'free') +
  scale_fill_continuous(limits = c(1, 3500), na.value = 'light grey') +
  theme_classic()

cazyme_summary <- cazyme_tpm %>% select(`Target Name`, `HMM Name`, `i-Evalue`, Coverage, locus_tag, length, product, Culture) %>% unique()

write.csv(cazyme_summary, 'cazyme_summary.csv', row.names = FALSE)
