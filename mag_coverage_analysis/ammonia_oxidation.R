
library(dplyr)
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(tidyverse)

# ammonia oxidation
setwd("mag_coverage_analysis")

exp <- read.csv("annotated_transcript_hits_tpm.csv")

## manual annotations

exp[exp$product == "methane/ammonia monooxygenase subunit C",]$product <- "ammonia monooxygenase subunit C"
exp[exp$product == "methane/ammonia monooxygenase subunit B",]$product <- "ammonia monooxygenase subunit B"
exp[exp$product == "methane/ammonia monooxygenase subunit A",]$product <- "ammonia monooxygenase subunit A"
exp[exp$product == "putative ammonia monooxygenase-assosciated protein",]$product <- "putative ammonia monooxygenase subunit X"

colors <- met.brewer('Thomas')

amm.list <- c(
  "Glutamine--fructose-6-phosphate aminotransferase [isomerizing]", "Nitrogen regulatory protein P-II", "putative nitrogen regulatory protein P-II", "Glutamate dehydrogenase",
  "Glutamine synthetase", "Glutamate synthase domain-containing protein","Glutamate--tRNA ligase", # ammonia assimilation 
  "PhnB protein", "Peroxyureidoacrylate/ureidoacrylate amidohydrolase RutB", # intracellular ammonia provisioning
  "Carbamoyl phosphate synthase small chain", "Carbamoylphosphate synthase large subunit", "carbamoyl phosphate synthetase subunit beta",
  "Urease accessory protein UreG" , "Urease subunit gamma" , "Urease subunit beta" , "Urease subunit alpha" ,
  "putative ammonium transporter MTH_663", "Ammonium transporter Amt3", 
  "Ammonia monooxygenase subunit Y", "putative ammonia monooxygenase subunit X", "ammonia monooxygenase subunit C", "ammonia monooxygenase subunit B", "ammonia monooxygenase subunit A")

amm <- exp[exp$product %in% amm.list,]

amm$product_short <- factor(amm$product, levels = c(
  "Glutamine--fructose-6-phosphate aminotransferase [isomerizing]", "Nitrogen regulatory protein P-II", "putative nitrogen regulatory protein P-II", "Glutamate dehydrogenase",
  "Glutamine synthetase", "Glutamate synthase domain-containing protein","Glutamate--tRNA ligase", # ammonia assimilation 
  "PhnB protein", "Peroxyureidoacrylate/ureidoacrylate amidohydrolase RutB", # intracellular ammonia provisioning
  "Carbamoyl phosphate synthase small chain", "Carbamoylphosphate synthase large subunit", "carbamoyl phosphate synthetase subunit beta",
  "Urease accessory protein UreG" , "Urease subunit gamma" , "Urease subunit beta" , "Urease subunit alpha" ,
  "putative ammonium transporter MTH_663", "Ammonium transporter Amt3", 
  "Ammonia monooxygenase subunit Y", "putative ammonia monooxygenase subunit X", "ammonia monooxygenase subunit C", "ammonia monooxygenase subunit B", "ammonia monooxygenase subunit A"),
  labels = c("GPI", "P-II", "P-II", "GDH", "GS", "GOGAT", "GTL", 
             "phnB", "rutB",
             "carA", "carB", "carB",
             "ureG", "ureC", "ureB", "ureA",
             "amt_MTH_663", "amt3",
             "amoY", "amoX", "amoC", "amoB", "amoA"))

amo.plot <- ggplot(amm[amm$product_short %in% c("amoA", "amoB", "amoC", "amoX", "amoY"),], aes(x = Sample, y = product_short, fill = TPM)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab('') +
  ylab('') +
  ggtitle('Ammonia Monooxygenase Transcripts') + 
  scale_fill_continuous(name = 'TPM', low = colors[3], high = colors[2], na.value = 'grey90', limits = c(0.1, 3e5)) +
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

amo.plot

other.amm.genes <- setdiff(amm, amm[amm$product_short %in% c("amoA", "amoB", "amoC", "amoX", "amoY"),])

other.amm <- ggplot(other.amm.genes, aes(x = Sample, y = product_short, fill = TPM)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab('') +
  ylab('') +
  ggtitle('Ammonia Oxidation Transcripts') + 
  scale_fill_continuous(name = 'TPM', low = colors[3], high = colors[2], na.value = 'grey90', limits=c(0.1,1e5)) +
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

other.amm

plot_grid(amo.plot, other.amm, align = "hv", nrow = 1)

# comparison with no3

nutrients <- read.csv('nutrients.pd.csv')

nitrogen <- nutrients %>% filter(depth == '05-10 cm', loc == 'PD', season2 %in% c('18-Apr', '18-Aug', '18-Oct'))
nitrogen$Sample <- NA
nitrogen[nitrogen$X == '18AprPD105',]$Sample <- 'Apr_PD1'
nitrogen[nitrogen$X == '18AprPD205',]$Sample <- 'Apr_PD2'
nitrogen[nitrogen$X == '18AprPD305',]$Sample <- 'Apr_PD3'
nitrogen[nitrogen$X == '18AugPD105',]$Sample <- 'Aug_PD1'
nitrogen[nitrogen$X == '18AugPD205',]$Sample <- 'Aug_PD2'
nitrogen[nitrogen$X == '18AugPD305',]$Sample <- 'Aug_PD3'
nitrogen[nitrogen$X == '18OctPD105',]$Sample <- 'Oct_PD1'
nitrogen[nitrogen$X == '18OctPD205',]$Sample <- 'Oct_PD2'
nitrogen[nitrogen$X == '18OctPD305',]$Sample <- 'Oct_PD3'

amo <- amm %>% filter(product_short %in% c('amoY', 'amoX', 'amoC', 'amoB', 'amoA'))
amo$season2 <- NA
amo[str_detect(amo$Sample, 'Apr') == TRUE,]$season2 <- '18-Apr'
amo[str_detect(amo$Sample, 'Aug') == TRUE,]$season2 <- '18-Aug'
amo[str_detect(amo$Sample, 'Oct') == TRUE,]$season2 <- '18-Oct'

amo.no3 <- left_join(amo, nitrogen, by = 'Sample')
amo.no3 <- amo.no3 %>% filter(product_short != 'amoX', product_short != 'amoY')

amo.no3$product_short <- factor(amo.no3$product_short, levels = c('amoA', 'amoB', 'amoC'))
amo.no3 <- amo.no3 %>% filter(locus_tag != 'EC11NC_v1_3459') # reads were mapped twice to two amoC in NC11

ggplot(amo.no3) +
  geom_boxplot(aes(x = season2.x, y = NO3*375), alpha = 0.3, col = 'gray') +
  geom_point(aes(x = season2.x, y = TPM, col = Culture), size = 2) + 
  stat_summary(aes(x = season2.x, y = TPM, col = Culture, group = Culture), fun = "mean", geom = "line", linetype = 'dashed') + 
  scale_y_continuous(
    name = "TPM",
    sec.axis = sec_axis(~./375, name = "NO3 (mg N/g DW soil)")) + 
  facet_wrap(vars(product_short)) + 
  theme_classic() + 
  xlab('') + 
  scale_color_manual(values = c(colors[2],colors[3],colors[8], colors[1])) + 
  theme(text = element_text(size = 14),
        legend.position = 'bottom')