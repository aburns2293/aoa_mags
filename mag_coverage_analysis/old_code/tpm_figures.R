
library(dplyr)
library(ggplot2)

cul <- read.csv('ec133_TPM.csv')

cul.expressed <- cul %>% filter(TPM > 0)
cul.expressed$ID <- str_c(cul.expressed$locus_tag, cul.expressed$product, sep=' ')

cul.expressed <- cul.expressed %>% group_by(ID) %>% arrange(desc(TPM), .by_group = TRUE)

cul.expressed$Sample <- factor(cul.expressed$Sample, levels = c('Apr_PD1', 'Apr_PD2', 'Apr_PD3', 'Aug_PD1', 'Aug_PD2', 'Aug_PD3', 'Oct_PD1', 'Oct_PD2', 'Oct_PD3'))

cul.expressed <- cul.expressed %>% group_by(ID) %>% mutate(Frequency = n())

ggplot(cul.expressed[cul.expressed$Frequency >= 5,], aes(x = Sample, y = ID, fill = TPM)) + 
  geom_tile() + 
  xlab('') + 
  ylab('') + 
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = 'TPM') + 
  geom_vline(xintercept = 3.5) + 
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(labels = c('', 'April', '', '', 'August', '', '', 'October', ''), expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme_classic() + 
  ggtitle('EC13.3 PD Transcripts') + 
  theme(panel.spacing = unit(0, 'lines'))

### ammonia oxidiation

cul <- read.csv('ec133_TPM.csv')

ao_genes <- c('ammonium transporter', 'ammonia monooxygenase', 'methane monooxygenase/ammonia monooxygenase subunit B', 'methane monooxygenase/ammonia monooxygenase subunit C', 
              'urea transporter', 'urease accessory protein UreD', 'urease accessory protein UreE', 'urease accessory protein UreF', 'urease accessory protein UreG', 'urease subunit alpha',
              'urease subunit beta', 'urease subunit gamma', 'cupredoxin domain-containing protein', 'multicopper oxidase domain-containing protein', 'multicopper oxidase family protein',
              'plastocyanin/azurin family copper-binding protein', 'copper resistance CopC family protein', 'copper resistance D family protein', 'copper resistance protein CopC')

#cul.ao <- cul[cul$product %in% ao_genes,]
cul.ao <- rbind(cul.ao, cul[cul$product %in% ao_genes,])

cul.ao.sums <- cul.ao %>% group_by(Culture, Sample, product) %>% summarise(Total_TPM = sum(TPM))

## three copies of amoC all had the exact same transcript counts mapped to them, this wasn't an issue with other multiple copy functional genes 
cul.ao.sums[cul.ao.sums$Culture == 'Cul11' & cul.ao.sums$product == 'methane monooxygenase/ammonia monooxygenase subunit C',]$Total_TPM <- cul.ao.sums[cul.ao.sums$Culture == 'Cul11' & cul.ao.sums$product == 'methane monooxygenase/ammonia monooxygenase subunit C',]$Total_TPM/3

cul.ao.sums$product <- factor(cul.ao.sums$product, levels = c('copper resistance D family protein', 'copper resistance CopC family protein', 'copper resistance protein CopC', 
                                                              'plastocyanin/azurin family copper-binding protein', 'multicopper oxidase domain-containing protein', 'multicopper oxidase family protein', 'cupredoxin domain-containing protein',
                                                              'urease subunit gamma', 'urease subunit beta', 'urease subunit alpha', 'urease accessory protein UreG', 'urease accessory protein UreF', 'urease accessory protein UreE',
                                                              'urease accessory protein UreD', 'urea transporter', 'ammonium transporter', 'methane monooxygenase/ammonia monooxygenase subunit C', 'methane monooxygenase/ammonia monooxygenase subunit B', 'ammonia monooxygenase'),
                              labels= c('CopD', 'CopC family protein', 'CopC', 'Plastocyanin/azurin family copper-binding protein', 'Multicopper oxidase domain-containing protein', 'Multicopper oxidase family protein', 'Cupredoxin domain-containing protein',
                                        'Urease subunit gamma', 'Urease subunit beta', 'Urease subunit alpha', 'UreG', 'UreF', 'UreE', 'UreD', 'Urea transporter', 'Amt', 'AmoC', 'AmoB', 'AmoA'))

#cul.ao.sums <- cul.ao.sums %>% filter(Total_TPM > 0)

ggplot(cul.ao.sums, aes(x = Sample, y = product, fill = Total_TPM)) + 
  geom_tile() + 
  facet_wrap(vars(Culture)) + 
  xlab('') + 
  ylab('') + 
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = 'TPM', na.value = 'grey90', limits = c(0.1, 3e5)) + 
  geom_vline(xintercept = 3.5) + 
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 9.5) +
  scale_x_discrete(labels = c('', 'April', '', '', 'August', '', '', 'October', ''), expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0), drop = FALSE) + 
  theme_classic() + 
  ggtitle('Ammonia Oxidation PD Transcripts', subtitle = 'Gray = present in genome but not expressed; white = absent in genome') + 
  theme(panel.spacing = unit(0, 'lines'))
