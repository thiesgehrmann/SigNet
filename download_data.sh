#!/bin/sh

mkdir data;
mkdir analysis;

# Install DTcut in current working directory

git clone git@github.com:thiesgehrmann/DTcut.git

# Download data
wget http://string-db.org/download/protein.links.v10/10090.protein.links.v10.txt.gz        -O data/10090.protein.links.v10.txt.gz
wget http://string-db.org/download/protein.aliases.v10/10090.protein.aliases.v10.txt.gz    -O data/10090.protein.aliases.v10.txt.gz
wget http://mutapedia.nki.nl/uploads/downloads/NCBIBuildM37_mm9_All_Genotypes_Tissue.csv   -O data/NCBIBuildM37_mm9_All_Genotypes_Tissue.csv
wget http://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz -O data/Mus_musculus.NCBIM37.67.gtf

# Unzip data
gunzip data/10090.protein.links.v10.txt.gz
gunzip data/10090.protein.aliases.v10.txt.gz

# Convert GTF to GFF data
gffread -E data/Mus_musculus.NCBIM37.67.gtf -o- > data/Mus_musculus.NCBIM37.67.gff
