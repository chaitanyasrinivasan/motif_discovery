# Installation


1. Clone this repository:

```shell
git clone https://github.com/chaitanyasrinivasan/motif_discovery.git
```

2. Download anaconda from https://docs.anaconda.com/anaconda/ to ensure your system has the compatible version of Python.

3. Use the following commands to install all dependencies for this pipeline.

```shell
  cd motif_discovery/scripts
  chmod 755 *.sh
  ./install.sh
```
# Dependencies

- [Python 3](https://docs.anaconda.com/anaconda/install/)
	- argparse
	- cython
	- matplotlib
	- numpy
	- pandas
	- pickle
	- seqlogo
	- sys
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [wget](https://www.gnu.org/software/wget/)
- [Slurm Workload Manager (for submitting parallel jobs)](https://slurm.schedmd.com/download.html)
- 2 GB of disk space

All dependencies excluding Slurm and wget are installed in steps above.

# **Motif Discovery**

This tool can perform de novo regulatory motif discovery on hg38 fasta sequences, BED coordinates, or marker genes.

![Image of tool](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/motif_discovery.jpg)

## Parameters

### Required Arguments

```
-i --input : .fasta file. The file should contain at least two sequences, and the sequences should have length greater than or equal to a specified width.  
-t --type : Input file type [FASTA/BED/GENES]. BED coordinates must be mapped to the hg38 reference genome. The GENES file type is a single column of official human gene names crosslisted with GENCODE Human Release 33.
```

The file extensions for the input passed to `-i` must be specified as follows:
TYPE | File Extension
------------ | -------------
FASTA | .fa
BED | .bed
GENES | .txt

### Optional Arguments

```
-w --width : Integer length of motif. If this flag is not provided, the motif width will be inferred.
-s --sequential : Run the motif discovery algorithm sequentially.
-p --parallel : Run the motif discovery algorithm in parallel. Requires Slurm Workload Manager.  
-h --help
```

If the `-s` or `-p` flags are not specified, the program will automatically choose the mode that gives better performance.

## Usage

```shell
cd motif_discovery/scripts
./gibbs_parallel.sh -i [/path/to/input] -t [FASTA/BED/GENES]
```
## Output

The program will output the motif logo and a plot of the motif entropy over sampling iterations. The motif logo represents the alignment with maximimum entropy.

![Image of motif](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/example_motif.png)

![Image of entropy](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/example_entropy.png)
