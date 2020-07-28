#!/bin/bash

echo "Installing python packages..."
#Setup channels
if ! [ -x "$(command -v conda)" ]
then
	echo "Error : conda is not executable."
	exit 0
else
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda install -c anaconda cython
	conda install -c anaconda pandas
	conda install -c conda-forge matplotlib
	conda install -c conda-forge argparse
	conda install -c bioconda seqlogo
fi
echo "Done"
echo "Downloading hg38 reference files..."
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
echo "Done"
echo "Cythonizing scripts..."
python setup.py build_ext --inplace
echo "Done"
echo "Installation complete."
exit 1
