
library(dplyr)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(MetBrewer)
library(cowplot)

setwd('~/aoa_mags/mag_coverage_analysis')

exp <- read.csv('transcript_hits_tpm.csv')

exp$product <- URLdecode(exp$product)

write.csv(exp, 'transcript_hits_tpm.csv', row.names = FALSE)

##### merge with annotated list

expressed_genes_list <- read_excel("expressed_genes_list.xlsx")

eg_annotated <- left_join(exp, expressed_genes_list, by = join_by('product' == 'Name'))

eg_annotated %>% filter(is.na(`Putative Pathway`) == TRUE) %>% dplyr::select(product) %>% unique() %>% count()

## manually editing incorrect merges

eg_annotated[eg_annotated$product == 'Blue (Type1) copper domain-containing protein',]$`Putative Pathway` <- 'Copper oxidation'
eg_annotated[eg_annotated$product == 'Blue (Type 1) copper domain protein',]$`Putative Pathway` <- 'Copper oxidation'
eg_annotated[eg_annotated$product == 'Blue (Type 1) copper domain protein (Modular protein)',]$`Putative Pathway` <- 'Copper oxidation'
eg_annotated[eg_annotated$product == '4Fe-4S binding domain protein',]$`Putative Pathway` <- 'Iron-sulfur cluster biosynthesis'
eg_annotated[eg_annotated$product == 'Dentin sialophosphoprotein',]$`Putative Pathway` <- 'Secondary metabolite biosynthesis'
eg_annotated[eg_annotated$product == '2,5-furandicarboxylate decarboxylase 1',]$`Putative Pathway` <- 'Secondary metabolism'
eg_annotated[eg_annotated$product == '7,8-dihydro-6-hydroxymethylpterin dimethyltransferase',]$`Putative Pathway` <- 'Secondary metabolite biosynthesis'
eg_annotated[eg_annotated$product == '7,8-dihydro-6-hydroxymethylpterin dimethyltransferase',]$Notes <- 'Methanopterin , found in Methanobacteria'
eg_annotated[eg_annotated$product == 'chromosome segregation protein',]$`Putative Pathway` <- 'Cell division'
eg_annotated[eg_annotated$product == 'Isocitrate dehydrogenase (NAD(+))',]$`Putative Pathway` <- 'TCA cycle'
eg_annotated[eg_annotated$product == 'beta-1,2-mannobiose phosphorylase / 1,2-beta-oligomannan phosphorylase',]$`Putative Pathway` <- 'Secondary metabolism'
eg_annotated[eg_annotated$product == 'putative dihydroorotate dehydrogenase B (NAD(+)), electron transfer subunit',]$`Putative Pathway` <- 'Respiration'

eg_annotated$Date <- str_sub(eg_annotated$Sample, 1, 3)

write.csv(eg_annotated, 'annotated_transcript_hits_tpm.csv', row.names = FALSE)

## evaluating different groups

ggplot(eg_annotated[eg_annotated$Culture != 'Nitrosotalea 5',], aes(x = Date, y = TPM, col = Culture)) + geom_point() + facet_wrap(vars(`Putative Pathway`), scales = 'free') +
  theme_classic()

## biofilm-related

eg_annotated <- read.csv('annotated_transcript_hits_tpm.csv')

biofilm <- eg_annotated %>% filter(Putative.Pathway %in% c('Adhesion', 'Motility', 'Quorum sensing', 'Signaling'), Culture %in% c('Nitrosotalea 131', 'Ca. Nitrosocosmicus')) %>% arrange(`Putative.Pathway`)

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
  scale_fill_continuous(name = 'TPM', low = colors[5], high = colors[6], na.value="grey90", limits = c(0.1, 50000)) +
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

## stress

stress <- eg_annotated %>% filter(Putative.Pathway %in% c('Chaperone', 'Nutrient storage', 'Stress protein', 'Antitoxin'), Culture %in% c('Nitrosotalea 131', 'Ca. Nitrosocosmicus'))

stress[stress$Culture == 'Nitrosotalea 131',]$Culture <- 'Nitrosotalea'

stress$Sample <- factor(stress$Sample, levels = c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3'))

ggplot(stress, aes(x = Sample, y = product, fill = TPM)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab('') +
  ylab('') +
  ggtitle('Stress-Related Transcripts') + 
  scale_fill_continuous(name = 'TPM Scaled', low = colors[5], high = colors[6]) +
  facet_wrap(vars(Culture)) + 
  theme_classic() +
  theme(aspect.ratio = 1.5,
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

stress <- stress %>% filter(product != 'PAC2 family protein', product != 'prefoldin alpha subunit', product != 'Prefoldin subunit alpha',
                            product != 'Prefoldin subunit beta', product != 'peptidyl-prolyl cis-trans isomerase B (cyclophilin B)',
                            product != 'OmpH family outer membrane protein')

stress$Type <- 'NA'

stress[stress$product %in% c('Antitoxin SocA-like Panacea domain-containing protein',
                             'Type II toxin-antitoxin system PemK/MazF family toxin'),]$Type <- 'Antitoxin'
stress[stress$product %in% c('General stress protein CsbD',
                             'Universal stress protein', 'Universal stress protein family protein', 'UspA domain-containing protein'),]$Type <- 'Universal'
stress[stress$product %in% c('Mannosyl-3-phosphoglycerate synthase', 'putative RNA-binding protein containing KH domain', 'serine protease Do'),]$Type <- 'Heat'
stress[stress$product %in% c('membrane protein HdeD'),]$Type <- 'Acid'
stress[stress$product %in% c('ATP-dependent RNA helicase DeaD', 'Ice-binding protein'),]$Type <- 'Cold'
stress[stress$product %in% c('DNA protection during starvation protein', 'gamma-glutamyl hercynylcysteine S-oxide synthase',
                             'methionine sulfoxide reductase B', 'Non-heme bromoperoxidase BpoC',
                             'putative peroxiredoxin bcp', 'Superoxide dismutase [Mn]', 'Superoxide dismutase [Mn]', 'thioredoxin 1', 'Hypoxic response protein 1'),]$Type <- 'Oxidative'
stress[stress$product %in% c('Phasin domain-containing protein', 'Phasin family protein', 'Phasin protein', 'putative Patatin-like phospholipase'),]$Type <- 'Nutrient storage'
stress[stress$Putative.Pathway == 'Chaperone',]$Type <- 'Heat shock chaperones'

stress.grouped <- stress %>% group_by(Culture, Date, Sample, Type) %>% summarise(total = sum(TPM))

stress_loci <- stress %>% select(product, locus_tag) %>% unique()

stress_plot <- ggplot(stress.grouped, aes(x = Sample, y = Type, fill = total)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab('') +
  ylab('') +
  ggtitle('Stress-Related Transcripts') + 
  scale_fill_continuous(name = 'TPM Scaled', low = colors[8], high = colors[7], limits = c(1, 75000), na.value="grey90") +
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

stress_plot

### water content information

water_table <- read.csv('water_table_summary.csv') %>% filter(Site == 'PD')library(lubridate)
library(stringr)
water_table$Date2 <- str_replace(water_table$Date, '^00', '20')
water_table$Date2 <- ymd(water_table$Date2)

water_table <- water_table %>% filter(Date2 > '2018-03-15' & Date2 < '2018-11-15')

water_fig <- ggplot(water_table, aes(x = Date2, y = GW_level)) + 
  geom_rect(aes(xmin = ymd('2018-04-01'), xmax = ymd('2018-04-30'),
                ymin= -85, ymax = 10), fill = '#EDC9D5', alpha = 0.03) + 
  geom_rect(aes(xmin = ymd('2018-08-01'), xmax = ymd('2018-08-30'),
                ymin= -85, ymax = 10), fill = '#EDC9D5', alpha = 0.03) + 
  geom_rect(aes(xmin = ymd('2018-10-01'), xmax = ymd('2018-10-30'),
                ymin= -85, ymax = 10), fill = '#EDC9D5', alpha = 0.03) + 
  scale_y_continuous(limits = c(-85, 10), expand = c(0,0)) + 
  geom_line() + 
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  xlab('') +
  scale_x_date(breaks = c(ymd('2018-04-15'),  ymd('2018-08-15'), ymd('2018-10-15')), date_labels = '%b', expand = c(0,0)) + 
  ggtitle('Water table depth (cm)') + 
  ylab('') + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        title = element_text(size = 14))
water_fig

redox <- read_xlsx('Redox_potential_all.xlsx', sheet = 'Redox_potential_all')
redox.pd <- redox %>% filter(Site == 'PD', Date %in% c('2018-04-20', '2018-08-09', '2018-10-02'))
redox.pd$Season3 <- NA
redox.pd[redox.pd$Date == '2018-04-20',]$Season3 <- 'Apr'
redox.pd[redox.pd$Date == '2018-08-09',]$Season3 <- 'Aug'
redox.pd[redox.pd$Date == '2018-10-02',]$Season3 <- 'Oct'
redox.pd[redox.pd$Season3 == 'Jun',]$`0-30 cm` <- NA

redox.pd$Season3 <- factor(redox.pd$Season3, levels = c('Apr', 'Jun', 'Aug', 'Oct'))

redox_fig <- ggplot(redox.pd, aes(x = Season3, y = `0-30 cm`)) + 
  geom_point(size = 3) +
  geom_line(group = 1, linewidth = 1.2) + 
  geom_hline(yintercept = 300, linetype = 'dashed') + 
  theme_classic() +
  xlab('') +
  ylab('') + 
  ggtitle('Redox potential (mV)') +
  scale_x_discrete(breaks = c('Apr', 'Jun', 'Aug', 'Oct'), labels = c('Apr', '', 'Aug', 'Oct'), drop = FALSE) + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        title = element_text(size = 14))

soil_temperature <- read_excel("soil_temperature.xlsx")

soil_temp_pd <- soil_temperature %>% filter(site == 'PD')
soil_temp_pd$datetime_utc <- ymd_hms(soil_temp_pd$datetime_utc)
soil_temp_pd <- soil_temp_pd %>% filter(datetime_utc > '2018-03-15' & datetime_utc < '2018-11-15')
soil_temp_pd <- soil_temp_pd %>% mutate(day = date(datetime_utc)) %>% group_by(day) %>% summarise(daily_avg = mean(soil_t_5))

temp_fig <- ggplot(soil_temp_pd, aes(x = day, y = daily_avg)) + 
  scale_y_continuous(limits = c(0,24), breaks = c(0, 5, 10, 15, 20), expand = c(0,0)) + 
  geom_rect(aes(xmin = ymd('2018-04-01'), xmax = ymd('2018-04-30'),
                ymin= 0, ymax = 24), fill = '#EDC9D5', alpha = 0.03) + 
  geom_rect(aes(xmin = ymd('2018-08-01'), xmax = ymd('2018-08-30'),
                ymin= 0, ymax = 24), fill = '#EDC9D5', alpha = 0.03) + 
  geom_rect(aes(xmin = ymd('2018-10-01'), xmax = ymd('2018-10-30'),
                ymin= 0, ymax = 24), fill = '#EDC9D5', alpha = 0.03) + 
  scale_x_date(breaks = c(ymd('2018-04-15'),ymd('2018-08-15'), ymd('2018-10-15')), date_labels = '%b', expand = c(0,0)) + 
  geom_line() + 
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  xlab('') +
  ylab('') + 
  ggtitle('Soil temperature (°C)') + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        title = element_text(size = 14))

temp_fig
## assembling figure

library(cowplot)
grid1 <- plot_grid(water_fig, temp_fig, redox_fig, align = 'hv', ncol=1)

grid2 <- plot_grid(biofilm_plot, stress_plot, align = 'hv', ncol = 1)

grid1
grid2

plot_grid(biofilm_plot, stress_plot, grid1, ncol = 1)
