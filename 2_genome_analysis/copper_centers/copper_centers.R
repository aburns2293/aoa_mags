
library(stringr)
library(dplyr)
library(ggplot2)
library(MetBrewer)


setwd("mag_coverage_analysis")

annot <- read.csv('annotated_transcript_hits_tpm.csv')

copper <- annot[str_detect(annot$product, "copper") == TRUE |
                  str_detect(annot$product, "Copper") == TRUE |
                  str_detect(annot$product, "Nitrite") == TRUE |
                  str_detect(annot$product, "nitrite") == TRUE,] %>%
  filter(Culture != "Nitrosotalea 5")

ggplot(copper, aes(x = Sample, y = product, fill = TPM)) +
  geom_tile() +
  facet_wrap(vars(Culture)) +
  scale_fill_continuous(name = 'TPM', na.value="grey90", limits = c(0.1, 30000)) +
  facet_wrap(vars(Culture)) + 
  theme_classic() +
  theme(
        panel.spacing=unit(0, "lines"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

copper.counts <- copper %>% group_by(locus_tag) %>% filter(row_number()==1)

copper.counts <- copper.counts %>% group_by(Culture, product) %>% summarise(count = n())

## systematizing names

multicopper_oxidases <- c('Multicopper oxidase', 'Multicopper oxidase mco', 'putative Multicopper oxidase')
blue_coppers <- c('Blue (Type 1) copper domain protein', 'Blue (type 1) copper domain-containing protein', 'Blue (Type 1) copper domain protein (Modular protein)',
                  'Blue (Type1) copper domain-containing protein', 'Blue (type 1) copper domain-containing protein', 'putative Blue (Type 1) copper domain protein')
copper_binding_plastocyanin <- c('Copper binding protein, plastocyanin/azurin family', 'Copper-binding protein, plastocyanin/azurin family protein', 
                                 'putative Copper binding protein, plastocyanin/azurin family',
                                 'putative copper-binding protein, plastocyanin/azurin family protein')
copd <- c('Copper resistance protein D', 'Copper resistance protein D domain-containing protein')

copper$product2 <- copper$product

copper[copper$product2 %in% multicopper_oxidases,]$product2 <- 'Multicopper oxidase'
copper[copper$product2 %in% blue_coppers,]$product2 <- 'Blue (Type 1) copper domain protein'
copper[copper$product2 %in% copper_binding_plastocyanin,]$product2 <- 'Copper binding protein, plastocyanin/azurin family'
copper[copper$product2 %in% copd,]$product2 <- 'Copper resistance protein D'

colors <- met.brewer('Thomas')

ggplot(copper, aes(x = Sample, y = product2, fill = TPM)) +
  geom_tile() +
  facet_wrap(vars(Culture)) +
  scale_fill_continuous(name = 'TPM', high = colors[1], low = '#e0b4a6' ,na.value="grey90", limits = c(0.1, 30000)) +
  facet_wrap(vars(Culture)) + 
  theme_classic() +
  theme(
    panel.spacing=unit(0, "lines"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

copper.counts <- copper %>% group_by(locus_tag) %>% filter(row_number()==1)

copper.counts <- copper.counts %>% group_by(Culture, product2) %>% summarise(count = n())
