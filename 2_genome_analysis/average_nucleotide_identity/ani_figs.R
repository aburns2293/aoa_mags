
library(pheatmap)

setwd("mag_coverage_analysis")

talea <- read.delim("talea_ani/ANIm_identity.tsv", row.names = 1)

pheatmap(talea)

cosmicus <- read.delim("cosmicus_ani/ANIm_identity.tsv", row.names=1)

pheatmap(cosmicus)
