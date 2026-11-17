
library(dplyr)
library(ggplot2)
library(tidyverse)

nutrients <- read.csv('mag_coverage_analysis/nutrients.pd.csv')

nitrogen <- nutrients %>% select(season2, depth, NH4, NO3)
nitrogen <- pivot_longer(nitrogen, cols = c(-season2, -depth), names_to = 'Nutrient', values_to = 'Volume')
nitrogen <- nitrogen %>% group_by(season2, depth, Nutrient) %>% summarise(avg_volume = mean(Volume), sd_volume = sd(Volume))

ggplot(nitrogen, aes(x = season2, y = avg_volume, col = Nutrient, shape = Nutrient)) + 
  geom_point(size = 3) + 
  geom_line(aes(group = Nutrient), linewidth = 1) + 
  geom_errorbar(aes(min = avg_volume, max = avg_volume + sd_volume, width = 0), linewidth = 1) + 
  facet_wrap(vars(depth), scales = "free") +
  xlab('') + 
  ylab('') + 
  theme_classic() 

nitrogen.drought <- nitrogen %>% filter(season2 %in% c("18-Apr", "18-Jun", "18-Aug", "18-Oct"))

ggplot(nitrogen.drought, aes(x = season2, y=avg_volume, col = Nutrient, shape = Nutrient)) +
  geom_point() +
  facet_wrap(vars(depth), scales = "free")
## relationship between amoC expression and NO3 variation

ao.sums$season2 <- NA
ao.sums[str_detect(ao.sums$Sample, 'Apr') == TRUE,]$season2 <- '18-Apr'
ao.sums[str_detect(ao.sums$Sample, 'Aug') == TRUE,]$season2 <- '18-Aug'
ao.sums[str_detect(ao.sums$Sample, 'Oct') == TRUE,]$season2 <- '18-Oct'

amo <- c('AmoC', 'AmoB', 'AmoA')

amoc.no3 <- left_join(ao.sums, nutrients, by = 'season2', relationship = 'many-to-many')
amoc.no3 <- amoc.no3 %>% filter(product %in% amo, depth == '05-10 cm')

amoc.no3$product <- factor(amoc.no3$product, levels = c('AmoA', 'AmoB', 'AmoC'), labels = c('amoA', 'amoB', 'amoC'))


ggplot(amoc.no3) +
  geom_boxplot(aes(x = season2, y = NO3*4e02), alpha = 0.3, col = 'gray') +
  geom_point(aes(x = season2, y = Total_TPM, col = Culture), size = 2) + 
  stat_summary(aes(x = season2, y = Total_TPM, col = Culture, group = Culture), fun = "mean", geom = "line", linetype = 'dashed') + 
  scale_y_continuous(
    name = "TPM",
    sec.axis = sec_axis(~ . / 4e02, name = "NO3 (mg N/g DW soil)")) + 
  facet_wrap(vars(product)) + 
  theme_classic() + 
  xlab('') + 
  scale_color_manual(values = c("#1CCAD8","#109648","#feb24c", '#f03b20')) + 
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        legend.position = 'bottom')
