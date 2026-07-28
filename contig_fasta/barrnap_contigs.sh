
conda activate barrnap

barrnap 11_contig_2.fasta --kingdom arc > nc_mag_barrnap_out.txt

conda activate hmmer

## adapted from https://github.com/aziele/pfam_scan

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

mkdir pfamdb
gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

hmmpress pfamdb/Pfam-A.hmm

hmmscan --domtblout nt131.domtblout -E 1e-5 pfamdb/Pfam-A.hmm EC131NT.1-Proteins.fasta 
hmmscan --domtblout nc11.domtblout -E 1e-5 pfamdb/Pfam-A.hmm EC11NC.1-Proteins.fasta 