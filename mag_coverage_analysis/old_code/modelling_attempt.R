
## amo3 modelling and variable importance

library(Hmisc)
library(corrplot)
library(pls)
library(caret)
library(vip)
library(tidymodels)
library(dplyr)

amo.no3.model <- amo.no3 %>% dplyr::select(locus_tag, Culture, Sample, TPM_Scaled_ReadsAssigned, product_short, season, water, TOM, NO3, NH4) %>% filter(product_short == 'amoC', locus_tag != 'EC11NC_v1_1858')


ggplot(amo.no3.model, aes(x = NO3, y = log(TPM_Scaled_ReadsAssigned))) + geom_point() + facet_wrap(vars(Culture), scales = 'free') + geom_smooth(method = 'lm')

amo.no3.model <- amo.no3.model %>% dplyr::select(-product_short, -locus_tag, -Sample)

hist(amo.no3.model$TPM_Scaled_ReadsAssigned)
hist(log(amo.no3.model$TPM_Scaled_ReadsAssigned))

model <- lm(data = amo.no3.model, log(TPM_Scaled_ReadsAssigned) ~ .,)
plot(model)
summary(model)

# partial least squares linear model

library(caret)
#BiocManager::install('mixOmics')
library(mixOmics)
#install.packages('plsmod')
library(plsmod)

set.seed(124)

amo.no3.model2 <- amo.no3.model %>% pivot_wider(names_from = Culture, values_from = TPM_Scaled_ReadsAssigned)

amo_split <- initial_split(amo.no3.model, prop = 0.90)
train <- training(amo_split)
test <- testing(amo_split)
folds <- vfold_cv(train, v=5)
grid <- expand.grid(num_comp = seq(1,7,1))

pls_model <- parsnip::pls(num_comp = tune()) %>% set_mode('regression') %>% set_engine('mixOmics')

pls_recipe <- 
  recipe(TPM_Scaled_ReadsAssigned ~ .,
         data = train) %>%
  step_log(TPM_Scaled_ReadsAssigned, base = 10) %>% 
  step_relevel(Culture, ref_level = 'Nitrosotalea 131') %>%
  step_normalize(everything()) %>%
  step_dummy(all_nominal_predictors())

pls_wflow <- workflow() %>% add_model(pls_model) %>% add_recipe(pls_recipe)

pls_tune <- pls_wflow %>%
  tune_grid(resamples = folds,
            grid = grid)
pls_tune %>% 
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  ggplot() + 
  geom_line(aes(x =num_comp, y = mean))

pls_tune %>% 
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  ggplot() + 
  geom_line(aes(x =num_comp, y = mean))

pls_tune_result <- parsnip::pls(num_comp = 6) %>% set_mode('regression') %>% set_engine('mixOmics')

pls_final <-
  workflow() %>%
  add_model(pls_tune_result) %>%
  add_recipe(pls_recipe) %>%
  fit(data = train)

pls_final



vip(pls_final)

amo.no3.model2 <- amo.no3.model %>% mutate(TPM_Scaled_ReadsAssigned = log(TPM_Scaled_ReadsAssigned)) %>% pivot_wider(names_from = Culture, values_from = TPM_Scaled_ReadsAssigned)
amo.no3.model2$season <- factor(amo.no3.model2$season)
cor.m <- rcorr(as.matrix(amo.no3.model2.m))
corrplot(cor.m$r, p.mat = cor.m$P, type = 'upper', tl.col = "black", sig.level = 0.05, insig = 'blank')


# partial least squares model with full dataset

annot <- read.csv('annotated_transcript_hits_tpm_scaled.csv')

annot[annot$product == "methane/ammonia monooxygenase subunit C",]$product <- "ammonia monooxygenase subunit C"
annot[annot$product == "methane/ammonia monooxygenase subunit B",]$product <- "ammonia monooxygenase subunit B"
annot[annot$product == "methane/ammonia monooxygenase subunit A",]$product <- "ammonia monooxygenase subunit A"
annot[annot$product == "putative ammonia monooxygenase-assosciated protein",]$product <- "putative ammonia monooxygenase subunit X"

ggplot(annot[annot$TPM_Scaled_ReadsAssigned > 1,], aes(x = TPM_Scaled_ReadsAssigned)) + geom_histogram()

annot.wide <- annot %>% pivot_wider(id_cols = c('Sample', 'Putative.Pathway', 'Date', 'product'), names_from = Culture, values_from = TPM_Scaled_ReadsAssigned)

NT11.list <- annot.wide$`Ca. Nitrosocosmicus`
new.NT11 <- c()
NT131.list <- annot.wide$`Nitrosotalea 131`
new.NT131 <- c()
NT133.list <- annot.wide$`Nitrosotalea 133`
new.NT133 <- c()
for(i in NT133.list) {
  if(length(i) > 1 & length(unique(i))==1) {
    new <- i[1]
  }
  else if(length(i) > 1 & length(unique(i)>1)) { new <- sum(i)}
  else if(length(i) == 0){ new <- NA}
  else {new <- i}
  new.NT133 <- c(new.NT133, new)
}
new.NT133

annot.new <- cbind(Sample = annot.wide$Sample, PutativePathway = annot.wide$Putative.Pathway, Product = annot.wide$product, Date = annot.wide$Date, 
                   Nitrosocosmicus = new.NT11, Nitrosotalea131 = new.NT131, Nitrosotalea133 = new.NT133) %>% as.data.frame()

annot.new$Nitrosocosmicus <- as.numeric(annot.new$Nitrosocosmicus)
annot.new$Nitrosotalea131 <- as.numeric(annot.new$Nitrosotalea131)
annot.new$Nitrosotalea133 <- as.numeric(annot.new$Nitrosotalea133)

annot.new$total=rowSums(cbind(annot.new$Nitrosocosmicus, annot.new$Nitrosotalea131, annot.new$Nitrosotalea133),na.rm=TRUE)
annot.new <- annot.new %>% filter(total > 0) %>% dplyr::select(-total)

annot.long <- annot.new %>% pivot_longer(cols = c(-'Sample', -'PutativePathway', -'Product', -'Date'), names_to = 'Culture', values_to = 'TPM_Scaled')

annot.long[annot.long$PutativePathway == 'Secondary metaboilte biosynthesis',]$PutativePathway <- 'Secondary metabolite biosynthesis'

annot.pathways <- annot.long %>% group_by(Sample, PutativePathway, Culture) %>% summarise(Total_TPM_Scaled = sum(TPM_Scaled, na.rm = TRUE))

annot.pathways.wide <- annot.pathways %>% pivot_wider(id_cols = c('Sample', 'Culture'), names_from = PutativePathway, values_from = Total_TPM_Scaled)

annot.pathways.wide[is.na(annot.pathways.wide)] <- 0

library(stringr)
library(vegan)

annot.pathways.wide$SampleID <- str_c(annot.pathways.wide$Sample, '_', annot.pathways.wide$Culture)

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

nitrogen <- nitrogen %>% dplyr::select(Sample, NO3, TOM, water)

annot.pathways.nutrients <- left_join(annot.pathways.wide, nitrogen, by = 'Sample')

### nmds



##### nitrogen model

annot.pathways.no3 <- annot.pathways.nutrients %>% remove_rownames() %>% column_to_rownames(var = 'SampleID') %>%  dplyr::select(-water, -TOM, -Culture, -Sample)

no3_folds <- vfold_cv(no3_train, v=5, repeats=5)
no3_grid <- expand.grid(num_comp = seq(1,15,1))

no3_pls_model <- parsnip::pls(num_comp = tune()) %>% set_mode('regression') %>% set_engine('mixOmics')

no3_pls_recipe <- 
  recipe(NO3 ~ .,
         data = annot.pathways.no3) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors())

no3_pls_recipe

no3_pls_wflow <- workflow() %>% add_model(no3_pls_model) %>% add_recipe(no3_pls_recipe)

no3_pls_tune <- no3_pls_wflow %>%
  tune_grid(resamples = no3_folds,
            grid = no3_grid)
no3_pls_tune %>% 
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  ggplot() + 
  geom_line(aes(x =num_comp, y = mean))

no3_pls_tune %>% 
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  ggplot() + 
  geom_line(aes(x =num_comp, y = mean))

pls_tune_result <- parsnip::pls(num_comp = 3) %>% set_mode('regression') %>% set_engine('mixOmics')

pls_final <-
  workflow() %>%
  add_model(pls_tune_result) %>%
  add_recipe(no3_pls_recipe) %>%
  fit(data = annot.pathways.no3)

pls_final

pls_fit <- extract_fit_parsnip(pls_final)$fit
vip <- vip(pls_fit)

