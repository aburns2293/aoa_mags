
library(stringr)
library(dplyr)


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
