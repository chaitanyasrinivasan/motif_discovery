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

4. Within the `motif_discovery` folder, install seqlogo using the following commands:
```shell
  git clone https://github.com/betteridiot/seqlogo.git
  cd seqLogo
  python setup.py install
  conda env create -f environment.yml
```
# Dependencies


# **Gibbs Sampler Motif Finder**

## Parameters

### Required Arguments

-w --width : Integer length of motif
-i --input : Preprocessed .fasta file that contains only lower case DNA nucleotide characters. The file should contain at least two sequences, and the sequences should have length greater than or equal to the specified width.
-t --type : Input file type [FASTA/BED/GENES]. The GENES file type is a single column of official human gene names crosslisted with GENCODE Human Release 33.


### Optional Arguments

-s --sequential : Run the motif discovery algorithm sequentially.
-p --parallel : Run the motif discovery algorithm in parallel.
-v --verbose
-h --help

If the -s or -p flags are not specified, the program will automatically choose the mode that gives better performance.

The preprocessing step will discard sequences that contain 'X' or 'N' and convert all characters to lowercase.

## Usage

```shell
cd scripts
python setup.py build_ext --inplace
./gibbs_parallel.sh -i [/path/to/input] -w [integer size]
```
