#!/bin/bash

#### Chaitanya Srinivasan #####
#### This script installs dependencies for the motif discovery pipeline #####


echo "Installing python packages..."

if ! [ -x "$(command -v conda)" ]
then
	echo "Error : conda is not executable."
	exit 0
else
	#Setup channels and install python packages
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda install cython
	conda install pandas
	conda install matplotlib
	conda install argparse
	conda install seqlogo
fi
if ! [ -x "$(command -v wget)" ]
then
	echo "Warning! wget is not executable. hg38 reference files were not downloaded."
	echo "Download the following files and place in scripts/"
	echo "GENCODE release 33: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz"
	echo "hg38 FASTA: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz "
	echo "hg38 chromosome sizes: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
else
	echo "Downloading hg38 reference files..."
	wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
	wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
	echo "Done"
fi
echo "Cythonizing scripts..."
python setup.py build_ext --inplace
echo "Installation complete."
exit 1
