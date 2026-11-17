#!/bin/bash

cd /home/burnsa/aoa_mags/mag_mt_mapping

# loop through fasta files of MAGs and create minimap2 dictionary
for f in *.fasta
do
  # get the base name of the file
  base=$(basename $f .fasta)
  # run the mapping script
  minimap2 -d $base.mmi $base.fasta
done

# loop through mag dictionaries and map reads to them
for f in *.mmi 
do
  # get the base name of the file
  dict_base=$(basename $f .mmi)
  # loop through each mrna sequence and map to the mag
  for r in metatranscriptomes/*.fasta
    do
      # get the base name of the file
      read_base=$(basename $r .fasta)
      # run the mapping script
      minimap2 -ax map-ont $dict_base.mmi $r > $dict_base.$read_base.sam
    done
done

# loop through sam files and convert to bam
for f in *.sam
do
  base=$(basename $f .sam)
  samtools view -Sb $base.sam > $base.bam
done

# sort and index the bam files
for f in *.bam
do
  base=$(basename $f .bam)
  samtools sort $base.bam -o $base.sorted.bam
  samtools index $base.sorted.bam
done

# run R script to use featureCounts
# run Python script to extract sequences for covered reads

cd aoa_mags/mag_mt_mapping/output

conda activate blast_env
makeblastdb -in aoa_sequence.fasta -dbtype prot -out aoa_db
for f in *.fasta
do
  blastx -query $f -db aoa_db -outfmt 6 -evalue 1e-5 -out aoa_blast_result.$f.txt
done