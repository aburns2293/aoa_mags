!#/bin/bash

makeblastdb -in phasin_NC.fasta -dbtype prot -out phasin_Nc_percolatus_db

blastp -query A4241.1-Proteins.fasta -db phasin_Nc_percolatus_db -out phasin_Nc_exaquare.out -evalue 1e-5 -outfmt 6

blastp -query MY3.1-Proteins.fasta -db phasin_Nc_percolatus_db -out phasin_Nc_oleophilus.out -evalue 1e-5 -outfmt 6

