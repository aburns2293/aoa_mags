library(dplyr)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(MetBrewer)
library(cowplot)

setwd('~/aoa_mags/mag_coverage_analysis')

## biofilm-related

eg_annotated <- read.csv('annotated_transcript_hits_tpm.csv')

biofilm <- eg_annotated %>% filter(Putative.Pathway %in% c('Adhesion', 'Motility', 'Quorum sensing', 'Signaling')) %>% arrange(`Putative.Pathway`)

biofilm_loci <- biofilm %>% select(product, locus_tag, Putative.Pathway) %>% unique()
biofilm$Sample <- factor(biofilm$Sample, levels = c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3'))
biofilm$product <- factor(biofilm$product, levels = c('3-hydroxy-2,4-pentadione 5-phosphate thiolase',
                                                      'Late endosomal/lysosomal adaptor and MAPK and MTOR activator 5',
                                                      'MEDS domain-containing protein',
                                                      'Adenylate and Guanylate cyclase catalytic domain protein',
                                                      'Adenylate cyclase 2',
                                                      'ArnR1-like winged helix-turn-helix domain-containing protein',
                                                      'Flagellin',
                                                      'Spondin_N',
                                                      'WxL domain-containing protein',
                                                      'HYR domain-containing protein',
                                                      'F5/8 type C domain protein',
                                                      'Archaeal Type IV pilin N-terminal domain-containing protein'))
biofilm$Type <- 'NA'

biofilm[biofilm$product %in% c('3-hydroxy-2,4-pentadione 5-phosphate thiolase'),]$Type <- 'Quorum sensing'
biofilm[biofilm$product %in% c('Late endosomal/lysosomal adaptor and MAPK and MTOR activator 5',
                               'MEDS domain-containing protein',
                               'Adenylate and Guanylate cyclase catalytic domain protein',
                               'Adenylate cyclase 2'),]$Type <- 'Signaling'
biofilm[biofilm$product %in% c('ArnR1-like winged helix-turn-helix domain-containing protein',
                               'Flagellin'),]$Type <- 'Motility'
biofilm[biofilm$product %in% c('Spondin_N',
                               'WxL domain-containing protein',
                               'HYR domain-containing protein',
                               'F5/8 type C domain protein',
                               'Archaeal Type IV pilin N-terminal domain-containing protein'),]$Type <- 'Adhesion'

biofilm[biofilm$Culture == 'Nitrosotalea 131',]$Culture <- 'Nitrosotalea'

biofilm.grouped <- biofilm %>% group_by(Culture, Date, Sample, Putative.Pathway) %>% summarise(total = sum(TPM))

colors <- met.brewer('Thomas')

biofilm_plot <- ggplot(biofilm.grouped, aes(x = Sample, y = Putative.Pathway, fill = total)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab('') +
  ylab('') +
  ggtitle('Biofilm-Related Transcripts') + 
  scale_fill_continuous(name = 'TPM', low = colors[5], high = colors[6], na.value="grey90", limits = c(0.1, 50000),
                        breaks=c(1e4, 2e4, 3e4, 4e4), labels=c("1e4", "2e4", "3e4", "4e4")) +
  facet_wrap(vars(Culture)) + 
  theme_classic() +
  theme(aspect.ratio = 1,
        panel.spacing=unit(0, "lines"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
biofilm_plot

biofilm.stats <- biofilm.grouped %>% ungroup()

biofilm %>% group_by(Culture, Putative.Pathway, product) %>% select(locus_tag) %>% unique() %>% summarise(n = n())

# NT biofilm statistics

ggplot(data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Adhesion",],
       aes(x = Date, y = total)) +
  geom_boxplot()

kruskal.test(total ~ Date, data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Adhesion",])
dunn_test(data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Adhesion",],
          formula = total ~ Date, p.adjust.method = "BH")


ggplot(data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Motility",],
       aes(x = Date, y = total)) +
  geom_boxplot()

kruskal.test(total ~ Date, data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Motility",])
dunn_test(data = biofilm.stats[biofilm.stats$Culture == "Nitrosotalea" & biofilm.stats$Putative.Pathway == "Motility",],
          formula = total ~ Date, p.adjust.method = "BH")