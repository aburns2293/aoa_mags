library(dplyr)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(MetBrewer)
library(cowplot)

setwd('~/aoa_mags/mag_coverage_analysis')

hits <- read.csv('transcript_hits_tpm.csv')

annotations <- read.csv("TableS2_DroughtResponseExpression.csv", skip=1)

annot.short <- annotations %>% select(Name, Putative.Pathway)

full_hits <- left_join(hits, annot.short, by=join_by('product' == 'Name'))

## stress

stress <- full_hits %>% filter(`Putative.Pathway` %in% c('Chaperone', 'Nutrient storage', 'Stress protein', 'Antitoxin'))

stress$Sample <- factor(stress$Sample, levels = c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3'))

# add phaC, phaE, phaP to nutrient storage

phaC <- c("EC11NC_v1_0009", "EC11NC_v1_0847", "EC11NC_v1_3324", "EC131NT_v1_1537", "EC131NT_v1_1757", "EC131NT_v1_2104")
phaE <- c("EC131NT_v1_2105", "EC11NC_v1_3300", "EC11NC_v1_1498", "EC11NC_v1_1498", "EC11NC_v1_0480", "EC11NC_v1_0445", "EC11NC_v1_0010")
phaP <- c("EC11NC_v1_0469", "EC11NC_v1_0549", "EC11NC_v1_1440", "EC11NC_v1_2801", "EC11NC_v1_2843", "EC11NC_v1_3027",
          "EC11NC_v1_3181", "EC11NC_v1_3248", "EC11NC_v1_3259", "EC11NC_v1_3363", "EC131NT_v1_2103", "EC131NT_v1_2142")

phasins <- c(phaC, phaE, phaP)

pha <- hits %>% filter(locus_tag %in% phasins)

pha$Putative.Pathway <- "Nutrient storage"

stress <- rbind(stress, pha)

stress[stress$Culture == 'Nitrosotalea EC131',]$Culture <- 'Nitrosotalea'

colors <- met.brewer('Thomas')

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

stress.full <- stress$product %>% unique()

stress.short <- c('dps', 'htrA', 'usp', 'phaP', 'hsp20', 'hrp1', 'thsB', 'mnSOD',
                  'RCF3', 'trx1', 'uspA','AHSA1',
                  'MPGS', 'pemK/mazF', 'dnaK', 'dnaJ', 'csbD', 'ths', 'PNPLA', 'phaP', 'shsp',
                  'thsA', 'ibp','phaC', 'phaE', 'phaE', "phaP")

stress.rename <- cbind(stress.full, stress.short) %>% as.data.frame() %>% rename(product = stress.full,
                                                                                 short_name = stress.short)

stress <- left_join(stress, stress.rename, by = 'product')

stress$Type <- 'NA'

stress[stress$product %in% c('Antitoxin SocA-like Panacea domain-containing protein',
                             'Type II toxin-antitoxin system PemK/MazF family toxin'),]$Type <- 'Antitoxin'
stress[stress$product %in% c('General stress protein CsbD',
                             'Universal stress protein', 'Universal stress protein family protein', 'UspA domain-containing protein'),]$Type <- 'Universal'
stress[stress$product %in% c('Mannosyl-3-phosphoglycerate synthase', 'putative RNA-binding protein containing KH domain', 'serine protease Do'),]$Type <- 'Heat'
stress[stress$product %in% c('ATP-dependent RNA helicase DeaD', 'Ice-binding protein'),]$Type <- 'Cold'
stress[stress$product %in% c('DNA protection during starvation protein', 'gamma-glutamyl hercynylcysteine S-oxide synthase',
                             'methionine sulfoxide reductase B', 'Non-heme bromoperoxidase BpoC',
                             'putative peroxiredoxin bcp', 'Superoxide dismutase [Mn]', 'Superoxide dismutase [Mn]', 'thioredoxin 1', 'Hypoxic response protein 1'),]$Type <- 'Oxidative'
stress[stress$product %in% c('Phasin domain-containing protein', 'Phasin family protein', 'Phasin protein', 'putative Patatin-like phospholipase', "polyhydroxyalkanoate synthase subunit PhaC", "polyhydroxyalkanoate synthase subunit PhaE"),]$Type <- 'Nutrient storage'
stress[stress$Putative.Pathway == 'Chaperone',]$Type <- 'Heat shock chaperones'
stress[stress$Putative.Pathway == 'Nutrient storage',]$Type <- 'Nutrient storage'

stress.grouped <- stress %>% group_by(Culture, Sample, Type) %>% summarise(total = sum(TPM))

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
  scale_fill_continuous(name = "TPM", low = colors[8], high = colors[7], limits = c(1, 150000), na.value="grey90",
                        breaks=c(3e4, 6e4, 9e4), labels=c("3e4", "6e4", "9e4")) +
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

# non-collapsed stress plot

stress.sum <- stress %>% group_by(Culture, Sample, Type, short_name) %>% summarise(total_TPM = sum(TPM))

stress.arranged <- rbind(stress.sum[stress.sum$Type == "Nutrient storage",], stress.sum[stress.sum$Type == "Universal",],
                         stress.sum[stress.sum$Type == "Oxidative",], stress.sum[stress.sum$Type == "Heat shock chaperones",],
                         stress.sum[stress.sum$Type == "Heat",], stress.sum[stress.sum$Type == "Cold",],
                         stress.sum[stress.sum$Type == "Antitoxin",])

stress.arranged.name <- stress.arranged$short_name %>% unique()

stress.sum$short_name <- factor(stress.sum$short_name, levels=stress.arranged.name)

stress_plot_2 <- ggplot(stress.sum, aes(x = Sample, y = short_name, fill = total_TPM)) + 
  geom_tile() + 
  geom_vline(xintercept = 3.5, linetype = 'dashed') + 
  geom_vline(xintercept = 6.5, linetype = 'dashed') + 
  scale_x_discrete(labels = c('', 'Apr', '', '', 'Aug', '', '', 'Oct', ''), expand = c(0,0), drop = FALSE) +
  scale_y_discrete(expand = c(0,0), limits=rev) + 
  xlab('') +
  ylab('') +
  ggtitle('Stress-Related Transcripts') + 
  scale_fill_continuous(name = "TPM", low = colors[8], high = colors[7], limits=c(1,100000),na.value="grey90",
                        breaks=c(8e4,6e4,4e4,2e4), labels=c("8e4", "6e4", "4e4", "2e4")) +
  facet_wrap(vars(Culture)) + 
  theme_classic() +
  theme(#aspect.ratio = 1,
        panel.spacing=unit(0, "lines"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

stress_plot_2

#### stress statistics - only on genes with meaningful presence across transcriptomes

stress.stats <- stress.sum %>% ungroup()


### water content information

water_table <- read.csv('water_table_summary.csv') %>% filter(Site == 'PD')
library(lubridate)
library(stringr)
library(readxl)
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

water_table %>% mutate(month = month(Date2), year = year(Date2)) %>% group_by(month, year) %>%
  summarise(mean = mean(GW_level), sd = sd(GW_level))

redox <- read_xlsx('Redox_potential_all.xlsx', sheet = 'Redox_potential_all')
redox$Date <- as.Date(redox$Date)
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

redox_fig

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

soil_temp_pd %>% mutate(month = month(day), year = year(day)) %>% group_by(month, year) %>% summarise(mean = mean(daily_avg), sd = sd(daily_avg))
## assembling figure

library(cowplot)
grid1 <- plot_grid(water_fig, temp_fig, redox_fig, align = 'hv', ncol=1)

grid2 <- plot_grid(biofilm_plot, stress_plot, align = 'hv', ncol = 1)

grid1
grid2

plot_grid(biofilm_plot, stress_plot, grid1, nrow = 1)
