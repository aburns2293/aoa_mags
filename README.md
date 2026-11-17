# Code and data repository for "Genomic insights into the physiology, niche differentiation, and stress responses of two novel peatland ammonia-oxidizing archaea" (Burns et al., submitted)

### Author: Anna Burns

### Last updated: 17.11.2026

### Contact: [anna.burns\@uni-greifswald.de](mailto:anna.burns@uni-greifswald.de){.email}

This repository contains the data and code for the genomic and metatranscriptomic analysis of two novel ammonia-oxidizing archaea (AOA) species, *Candidatus* Nitrosocosmicus percolatus Mir11 and *Candidatus* Nitrosotalea inundata Mir13. If you use or refer to this repository, please cite the following manuscript:

Burns A, Pozaric F, Kerou M, Melcher M, Schlueter R, Zak D, Hagenau L, Kuss AW, Schleper C, Urich T & Wang H. Genomic insights into the physiology, niche differentiation, and stress responses of two novel peatland ammonia-oxidizing archaea. Submitted 2026.

The repository is organized as follows:

- 0_raw_genome_data:

  - Consensus metagenome-assembled genomes for both species (.fasta)

  - 16S rRNA genes from both species (BothGenomes_16SrRNAgenes.fasta)

  - Annotated genomes from MaGe MicroScope (.gff)

  - Available genomes from other *Nitrosotalea* and *Ca.* Nitrosocosmicus species, as accessed for the manuscript from MaGe MicroScope (.zip)

- 1a_seq_metatranscriptomes: all archaeal mRNA transcripts from 9 samples extracted from the Tribsees percolation fen peatland, sampled in triplicates in April, August and October 2018 (.fasta)

- 1b_seq_16Samplicons: amplicon sequencing data of the enrichment cultures, including:

  - aoa_seqtab_final.rds: reads assigned to ASVs

  - aoa_tax_final.rds: reads assigned to ASVs and taxonomy, based on SILVA 138.2 database

- 2_genome_analysis: all code and results from the genomic analysis, including:

  - average_nucleotide_identity

  - barrnap (identification of rRNA in genomes)

  - CAZYmes

  - copper_centers

  - phasins

- 3_mapping: scripts for figures related to ammonia_oxidation.R (Fig. 2) and biofilm_stress_analysis.R (Fig. 4), along with the code/data for:

  - TPM transformation (annotated_transcript_hits_tpm.csv, transcript_hits_tpm.csv, tpm_scaling.R)

  - Expressed gene pathways (expressed_genes_list.csv, expressed_genes_list.xlsx)

  - bowtie2 mapping output and code (mapping_output)

  - the results from RSubRead FeatureCounts (featurecounts)

- 4_environmental_data: nutrient, redox, soil temperature and water table summary data from the Tribsees percolation fen site
