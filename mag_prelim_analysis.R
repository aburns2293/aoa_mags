
library(ggplot2)
library(dplyr)
library(knitr)

hits <- read.csv('EnrichmentCultureMAGS_BlastHits.csv')

ggplot(hits, aes(x = Perc_ID, fill = Subject)) + geom_histogram() + facet_wrap(vars(Query)) + theme_classic() +
  xlab('Percent Identity') + ylab('Count') + ggtitle('Enrichment Culture Hits vs. Wetscapes 16S rRNA Amplicon Sequences') + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        title = element_text(size = 22))

pos.id <- hits %>% filter(Perc_ID == 100)

knitr::kable(pos.id, 'rst')

# import phyloseq

aoa.ps <- readRDS('wetscapes_aoa.rds')

## culture 5 (ASV3)

cul5 <- subset_taxa(aoa.ps, taxa_names(aoa.ps) == 'ASV3')

cul5.melt <- psmelt(cul5)

cul5.melt[cul5.melt$season2 == '18Feburary',]$season2 <- '18-Feb'
cul5.melt[cul5.melt$season2 == '19Feburary',]$season2 <- '19-Feb'

ggplot(cul5.melt, aes(x = depth, y = Abundance*100, col = Seq_Round)) + geom_boxplot() + facet_wrap(vars(loc), scales = 'free') + theme_classic() + xlab('') + ylab('Relative Abundance (%)')

ggplot(cul5.melt[cul5.melt$loc == 'PD',], aes(x = season2, y = Abundance*100, col = Seq_Round)) + geom_point() + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'PD Culture 5 (ASV3, Ca. Nitrosotalea) Relative Abundance')

ggplot(cul5.melt[cul5.melt$loc == 'CD',], aes(x = season2, y = Abundance*100, col = Seq_Round)) + geom_point() + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'CD Culture 5 (ASV3, Ca. Nitrosotalea) Relative Abundance')

write.csv(cul5.melt, 'culture5_relabundance.csv')

## culture 11 (ASV29, ASV57854)

cul11 <- subset_taxa(aoa.ps, taxa_names(aoa.ps) == 'ASV29' | taxa_names(aoa.ps) == 'ASV57854')

cul11.melt <- psmelt(cul11)

cul11.melt[cul11.melt$season2 == '18Feburary',]$season2 <- '18-Feb'
cul11.melt[cul11.melt$season2 == '19Feburary',]$season2 <- '19-Feb'

ggplot(cul11.melt, aes(x = depth, y = Abundance*100, col = Seq_Round)) + geom_boxplot() + facet_wrap(vars(loc), scales = 'free') + theme_classic() + xlab('') + ylab('Relative Abundance (%)')

ggplot(cul11.melt[cul11.melt$loc == 'PD',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size = 3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'PD Culture 11 (Ca. Nitrosocosmicus) Relative Abundance')

ggplot(cul11.melt[cul11.melt$loc == 'CD',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size = 3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'CD Culture 11 (Ca. Nitrosocosmicus) Relative Abundance')

ggplot(cul11.melt[cul11.melt$loc == 'PW',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size = 3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'PW Culture 11 (Ca. Nitrosocosmicus) Relative Abundance')

write.csv(cul11.melt, 'culture11_relabundance.csv')

## culture 13 (ASV3, ASV92)

cul13 <- subset_taxa(aoa.ps, taxa_names(aoa.ps) == 'ASV3' | taxa_names(aoa.ps) == 'ASV92')

cul13.melt <- psmelt(cul13)

cul13.melt[cul13.melt$season2 == '18Feburary',]$season2 <- '18-Feb'
cul13.melt[cul13.melt$season2 == '19Feburary',]$season2 <- '19-Feb'

ggplot(cul13.melt, aes(x = depth, y = Abundance*100, col = Seq_Round)) + geom_boxplot() + facet_wrap(vars(loc), scales = 'free') + theme_classic() + xlab('') + ylab('Relative Abundance (%)')

ggplot(cul13.melt[cul13.melt$loc == 'PD',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size = 3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'PD Culture 13 (Ca. Nitrosotalea) Relative Abundance')

ggplot(cul13.melt[cul13.melt$loc == 'CD',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size = 3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'CD Culture 13 (Ca. Nitrosotalea) Relative Abundance')

ggplot(cul13.melt[cul13.melt$loc == 'PW',], aes(x = season2, y = Abundance*100, col = Seq_Round, shape = OTU)) + geom_point(size =3) + facet_wrap(vars(depth), ncol = 1) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12)) + 
  labs(x = '', y = 'Relative Abundance (%)', title = 'PW Culture 13 (Ca. Nitrosotalea) Relative Abundance')

write.csv(cul13.melt, 'culture13_relabundance.csv')
