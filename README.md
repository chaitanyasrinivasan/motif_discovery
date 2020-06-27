# Installation

1. Download anaconda from https://docs.anaconda.com/anaconda/ to ensure your system has the compatible version of Python along with the necessary packages.

2. Install Cython with the following command:

```shell
   conda install -c anaconda cython
```
3. Clone this repository:

```shell
git clone https://github.com/chaitanyasrinivasan/motif_discovery.git
```

4. Within the motif_discovery folder, install seqlogo using the following commands:
```shell
  git clone https://github.com/betteridiot/seqlogo.git
  cd seqLogo
  python setup.py install
  conda env create -f environment.yml
```

# **Gibbs Sampler Motif Finder**

## Parameters

-w --width : Integer length of motif
-i --input : Processed .fasta file that contains only lower case DNA nucleotide characters. The file should contain at least 2 sequences, and all sequences should have length greater than or equal to the width. Raw .fasta files can be processed into a suitable format by placing them in the data/ folder and running the following command:

```shell
cd scripts
./preprocess.sh
```

