# Installation


1. Clone this repository:

```shell
git clone https://github.com/chaitanyasrinivasan/motif_discovery.git
```

2. Download anaconda from https://docs.anaconda.com/anaconda/ to ensure your system has the compatible version of Python along with the necessary packages.

3. Use the following commands to install all dependencies for this pipeline.

```shell
  cd motif_discovery/scripts
  chmod 755 *.sh
  ./install.sh
```
# Dependencies

- Python 3
	- argparse
	- cython
	- matplotlib
	- numpy
	- pandas
	- pickle
	- seqlogo
	- sys
- bedtools
- Slurm Workload Manager (for submitting parallel jobs)

All dependencies excluding Slurm will be installed in the third step above.

# **Motif Discovery**

This tool can perform de novo regulatory motif discovery on hg38 fasta sequences, BED coordinates, or marker genes.

![Image of tool](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/motif_discovery.jpg)

## Parameters

### Required Arguments

```
-w --width : Integer length of motif
-i --input : Preprocessed .fasta file that contains only lower case DNA nucleotide characters. The file should contain at least two sequences, and the sequences should have length greater than or equal to the specified width.  
-t --type : Input file type [FASTA/BED/GENES]. BED coordinates must be mapped to the hg38 reference genome. The GENES file type is a single column of official human gene names crosslisted with GENCODE Human Release 33.
```

### Optional Arguments

```
-s --sequential : Run the motif discovery algorithm sequentially.
-p --parallel : Run the motif discovery algorithm in parallel. Requires Slurm Workload Manager.  
-h --help
```

If the `-s` or `-p` flags are not specified, the program will automatically choose the mode that gives better performance.

## Usage

```shell
cd motif_discovery/scripts
./gibbs_parallel.sh -i [/path/to/input] -w [integer size] -t [FASTA/BED/GENES]
```
## Output

The program will output the motif logo and a plot of the loss function.