
library(dplyr)
library(ggplot2)
library(MetBrewer)
setwd('~/aoa_mags')
colors <- met.brewer('Thomas')
hits <- read.csv('EnrichmentCultureMAGS_BlastHits.csv')

pos.id <- hits %>% filter(Perc_ID == 100)

knitr::kable(pos.id, 'rst')

aoa.ps <- readRDS('wetscapes_aoa.rds')

## when looking within PD and shallow samples R2 = 0.285 for sequencing prep...

## all ASVs

cultures <- subset_taxa(aoa.ps, taxa_names(aoa.ps) %in% c('ASV3', 'ASV29'))

cultures.melt <- psmelt(cultures)

cultures.melt[cultures.melt$season2 == '18Feburary',]$season2 <- '18-Feb'
cultures.melt[cultures.melt$season2 == '19Feburary',]$season2 <- '19-Feb'


## figures

cultures.melt.avgs <- cultures.melt %>% group_by(OTU, loc, season2, depth, Seq_Round) %>% summarise(Mean = mean(Abundance), sd = sd(Abundance))

cultures.melt.avgs[cultures.melt.avgs$OTU == 'ASV29',]$OTU <- 'Ca. Nitrosocosmicus'
cultures.melt.avgs[cultures.melt.avgs$OTU == 'ASV3',]$OTU <- 'Nitrosotalea'

culture.melt.old <- cultures.melt.avgs %>% filter(season2 %in% c('17-Apr', '17-May', '17-Jun', '17-July', '17-Aug', '17-Sep', '17-Oct', '17-Nov', '17-Dec',
                                                                '18-Jan', '18-Feb', '18-Mar', '18-Apr', '18-May', '18-Jun', '18-Jul', '18-Aug', '18-Sep', '18-Oct', '18-Nov', '18-Dec',
                                                                '19-Jan', '19-Feb', '19-Mar', '19-Apr', '19-May', '19-Jun', '19-Jul', '19-Aug', '19-Sep', '19-Oct'))

culture.melt.old$season2 <- factor(culture.melt.old$season2, levels = c('17-Apr', '17-May', '17-Jun', '17-July', '17-Aug', '17-Sep', '17-Oct', '17-Nov', '17-Dec',
                                                                        '18-Jan', '18-Feb', '18-Mar', '18-Apr', '18-May', '18-Jun', '18-Jul', '18-Aug', '18-Sep', '18-Oct', '18-Nov', '18-Dec',
                                                                        '19-Jan', '19-Feb', '19-Mar', '19-Apr', '19-May', '19-Jun', '19-Jul', '19-Aug', '19-Sep', '19-Oct'))
  
ggplot(culture.melt.old[culture.melt.old$loc == 'PD' & culture.melt.old$Seq_Round == 'Old' & culture.melt.old$depth == '25-30 cm',], 
       aes(x = season2, y = Mean*100, col = OTU, shape = OTU, group = OTU)) + 
  geom_point(size = 4) + 
  geom_line() + 
  facet_wrap(vars(depth), ncol = 1) + 
  geom_errorbar(aes(ymax = Mean*100 + sd*100, ymin = Mean*100, color = OTU), width = 0.1, linewidth = 1) +
  scale_x_discrete(drop = FALSE, labels = c('17-Apr', '', '', '', '17-Aug', '', '', '17-Nov', '', '', '18-Feb', '', '18-Apr', '', '18-Jun', '', '18-Aug', '', '18-Oct', '', '', '', '19-Feb',
                                            '', '', '', '19-Jun', '', '', '', '19-Oct')) + 
  xlab('') + 
  ylab('') + 
  scale_y_continuous(limits = c(0,0.5)) + 
  scale_color_manual(values = c(colors[2], colors[6]), name = 'Culture') + 
  scale_shape_manual(values = c(16,17), name = 'Culture') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5)

cultures.melt.new <- cultures.melt.avgs %>% filter(Seq_Round == 'New')

cultures.melt.new$season2 <- factor(cultures.melt.new$season2, levels = c('20-Sep', '20-Oct', '20-Nov', '20-Dec',
                                                                          '21-Jan', '21-Feb', '21-Mar', '21-Apr', '21-May', '21-Jun', '21-Jul', '21-Aug', '21-Sep', '21-Oct', '21-Nov', '21-Dec',
                                                                          '22-Jan', '22-Feb', '22-Mar', '22-Apr', '22-May', '22-Jun', '22-Jul', '22-Aug', '22-Sep', '22-Oct', '22-Nov', '22-Dec',
                                                                          '23-Jan', '23-Feb', '23-Mar', '23-Apr', '23-May', '23-Jun' , '23-Jul', '23-Aug', '23-Sep', '23-Oct', '23-Nov', '23-Dec',
                                                                          '24-Jan', '24-Feb', '24-Mar', '24-Apr', '24-May', '24-Jun' ))

ggplot(cultures.melt.new[cultures.melt.new$loc == 'PD',], 
       aes(x = season2, y = Mean*100, col = OTU, shape = OTU, group = OTU)) + 
  geom_point(size = 4) + 
  geom_line() + 
  facet_wrap(vars(depth), ncol = 1) + 
  geom_errorbar(aes(ymax = Mean*100 + sd*100, ymin = Mean*100, color = OTU), width = 0.1, linewidth = 1) +
  #scale_x_discrete(drop = FALSE, labels = c('20-Sep', '', '', '', '21-Jan', '', '', '', '', '21-Jun', '', '', '', '', '', '', '22-Jan', '', '', '', '', '22-Jun', '', '', '', '', '', '', '23-Jan',
  #                                          '', '', '', '', '23-Jun', '', '', '', '', '', '', '24-Jan', '', '', '', '', '24-Jun')) + 
  xlab('') + 
  ylab('') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_color_manual(values = c(colors[2], colors[6]), name = 'Culture') + 
  scale_shape_manual(values = c(16,17), name = 'Culture') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5)

