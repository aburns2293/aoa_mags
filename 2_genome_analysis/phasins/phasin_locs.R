
library(ape)
library(dplyr)
library(ggbio)
library(stringr)

nt <- read.gff("mag_annotation_mage/Nitrosotalea_sp._EC131.gff")

seqinfo(nt)

nt <- nt %>% mutate(length=end-start)

nt.phasins <- c("EC131NT_v1_1537", "EC131NT_v1_1757", "EC131NT_v1_2103", "EC131NT_v1_2104", "EC131NT_v1_2105", "EC131NT_v1_2142")
for(i in 1:length(nt.phasins)){
  nt.phasins[i] <- paste0("locus_tag=", nt.phasins[i])
}
nt.phasins.data <- nt[nt$attributes %in% nt.phasins,]

write.csv(nt.phasins.data, "nt.phasins.data.csv", row.names=FALSE)

nt.phasins.ggbio <- data.frame(locus_tag=c("phaC", "phaC", "phaP", "phaC", "phaE", "phaP"),
                               start=nt.phasins.data$start[1:6],
                               end=nt.phasins.data$end[1:6],
                               strand=nt.phasins.data$strand[1:6])

gr.nt <- GRanges(seqnames="contig_1",
                 ranges=IRanges(start = nt.phasins.data$start, end=nt.phasins.data$end),
                 strand=nt.phasins.data$strand,
                 locus_tag =c("phaC", "phaC", "phaP", "phaC", "phaE", "phaP") )

autoplot(gr.nt, arrow=TRUE) + 
  theme_minimal()

