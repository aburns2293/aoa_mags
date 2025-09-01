
hits <- read.csv('EnrichmentCultureMAGS_BlastHits.csv')

pos.id <- hits %>% filter(Perc_ID == 100)

knitr::kable(pos.id, 'rst')

aoa.ps <- readRDS('wetscapes_aoa.rds')

## influence of sequencing round on aoa generally

library(vegan)

aoa.otu <- otu_table(aoa.ps) %>% data.frame()
aoa.sample <- sample_data(aoa.ps) %>% data.frame()
aoa.dist <- vegdist(aoa.otu, method = 'bray') ## lots of zeros, some empty rows

# results in negative distances, removing samples with no aoa detected (54 samples of 520)

aoa.positive <- aoa.otu[rowSums(aoa.otu) > 0,]

# check within PD site

aoa.sample.pd <- aoa.sample %>% filter(loc == 'PD', depth == '05-10 cm')
pd.sample.names <- rownames(aoa.sample.pd)
aoa.otu.pd <- aoa.otu[rownames(aoa.otu) %in% pd.sample.names,]

aoa.dist <- vegdist(aoa.otu.pd, method = 'bray')

aoa.positive.samples <- rownames(aoa.positive)
aoa.sample.positive <- aoa.sample[rownames(aoa.sample) %in% aoa.positive.samples,]
adonis2(aoa.dist ~ aoa.sample.pd$Seq_Round, permutations = 1000)
## significant (p = 0.000999 ***) but R2 only 0.032

## when looking within PD and shallow samples R2 = 0.285 for sequencing prep...

## all ASVs

cultures <- subset_taxa(aoa.ps, taxa_names(aoa.ps) %in% c('ASV3', 'ASV29', 'ASV57854', 'ASV92'))

cultures.melt <- psmelt(cultures)

cultures.melt[cultures.melt$season2 == '18Feburary',]$season2 <- '18-Feb'
cultures.melt[cultures.melt$season2 == '19Feburary',]$season2 <- '19-Feb'

cultures.melt <- cultures.melt %>% filter(OTU != 'ASV57854')

cultures.melt$season2 <- factor(cultures.melt$season2, levels = c('17-Apr', '17-May', '17-Jun', '17-July', '17-Aug', '17-Sep', '17-Oct', '17-Nov', '17-Dec',
                                                                  '18-Jan', '18-Feb', '18-Mar', '18-Apr', '18-May', '18-Jun', '18-Jul', '18-Aug', '18-Sep', '18-Oct', '18-Nov', '18-Dec',
                                                                  '19-Jan', '19-Feb', '19-Mar', '19-Apr', '19-May', '19-Jun', '19-Jul', '19-Aug', '19-Sep', '19-Oct'))

## figures

cultures.melt.avgs <- cultures.melt %>% group_by(OTU, loc, season2, depth, Seq_Round) %>% summarise(Mean = mean(Abundance), sd = sd(Abundance))

ggplot(cultures.melt[cultures.melt$loc == 'PD' & cultures.melt$Seq_Round == 'Old' & cultures.melt$depth == '05-10 cm',], aes(x = season2, y = Abundance*100, col = OTU, shape = OTU)) + 
  geom_point(size = 3) + 
  facet_wrap(vars(depth), ncol = 1) + 
  scale_x_discrete(drop = FALSE) + 
  xlab('') + 
  ylab('') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_color_manual(values = c("#1CCAD8","#109648","#D68C45")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5)

ggplot(cultures.melt.avgs[cultures.melt.avgs$loc == 'PD' & cultures.melt.avgs$Seq_Round == 'Old' & cultures.melt.avgs$depth == '25-30 cm',], 
       aes(x = season2, y = Mean*100, col = OTU, shape = OTU, group = OTU)) + 
  geom_point(size = 4) + 
  geom_line() + 
  facet_wrap(vars(depth), ncol = 1) + 
  geom_errorbar(aes(ymax = Mean*100 + sd*100, ymin = Mean*100, color = OTU), width = 0.1, linewidth = 1) +
  scale_x_discrete(drop = FALSE) + 
  xlab('') + 
  ylab('') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_color_manual(values = c("#1CCAD8","#109648","#D68C45")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5)

cultures.melt.new <- cultures.melt.avgs %>% filter(Seq_Round == 'New')

cultures.melt.new$season2 <- factor(cultures.melt.new$season2, levels = c('20-Sep', '21-Jan', '21-Jun', '22-Jan', '22-Jun', '23-Jan', '23-Jun', '24-Jan', '24-Jun'))

ggplot(cultures.melt.new[cultures.melt.new$loc == 'PD',], 
       aes(x = season2, y = Mean*100, col = OTU, shape = OTU, group = OTU)) + 
  geom_point(size = 4) + 
  geom_line() + 
  facet_wrap(vars(depth), ncol = 1) + 
  geom_errorbar(aes(ymax = Mean*100 + sd*100, ymin = Mean*100, color = OTU), width = 0.1, linewidth = 1) +
  scale_x_discrete(drop = FALSE) + 
  xlab('') + 
  ylab('') + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_color_manual(values = c("#1CCAD8","#109648","#D68C45")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5)

