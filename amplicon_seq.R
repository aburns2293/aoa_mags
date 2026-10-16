
library(ggplot2)

tax <- readRDS("aoa_tax_final.rds")
seq <- readRDS("aoa_seqtab_final.rds")

library(phyloseq)
library(tidyverse)

tax <- tax %>% as.data.frame() %>% rownames_to_column("Sequence")

tax$ASV <- NA
for(i in 1:nrow(tax)){
  tax[i,]$ASV <- paste0("ASV", i)
}

seq <- seq %>% t() %>% as.data.frame() %>% rownames_to_column("Sequence")

tax_asv <- tax %>% select(ASV, Sequence)

seq <- left_join(seq, tax_asv, by = "Sequence")
seq <- seq %>% select(-Sequence)

seq <- seq %>% column_to_rownames(var = "ASV")

tax <- tax %>% column_to_rownames(var = "ASV")

aoa.ps <- phyloseq(otu_table(seq, taxa_are_rows = TRUE), tax_table(as.matrix(tax)))

aoa.ps.rel <- transform_sample_counts(aoa.ps, function(otu) otu/sum(otu))

aoa.ps.clean <- psmelt(aoa.ps.rel)

nt <- aoa.ps.clean %>% filter(Sample %in% c("S5", "S6", "S11", "S12"))
nt_max <- nt %>% filter(Sample == "S5")
nt_max <- nt_max %>% filter(Abundance > 0) %>% mutate(Abundance = Abundance*100)
nt_contam <- nt_max %>% filter(Class == "Nitrososphaeria")

ggplot(nc, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col()

nc <- aoa.ps.clean %>% filter(Sample %in% c("S1", "S2", "S7", "S8"))

nc_max <- nc %>% filter(Sample == "S2")
nc_max <- nc_max %>% filter(Abundance > 0) %>% mutate(Abundance = Abundance*100)
nc_contam <- nc_max %>% filter(Class == "Nitrososphaeria")

plot(aoa.ps.clean)
write.fasta(as.list(tax$Sequence), tax$ASV, "ampseq.fasta")
