#!/bin/bash

echo "Installing python packages..."
conda install -c anaconda cython
conda install -c anaconda pandas
conda install -c conda-forge matplotlib
conda install -c conda-forge argparse
echo "Done"
echo "Downloading hg38 GENCODE and fasta files..."
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
echo "Done"
echo "Cythonizing scripts..."
python setup.py build_ext --inplace
echo "Done"
echo "Installing seqlogo"
cd ../
git clone https://github.com/betteridiot/seqlogo.git
cd seqLogo
python setup.py install
conda env create -f environment.yml
echo "Installation complete."
exit 1
