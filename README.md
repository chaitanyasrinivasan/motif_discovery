# Installation


1. Clone this repository:

```shell
git clone https://github.com/chaitanyasrinivasan/motif_discovery.git
```

2. Download anaconda from https://docs.anaconda.com/anaconda/ to ensure your system has the compatible version of Python.

3. Use the following commands to install dependencies for this pipeline.

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
- [**wget**](https://www.gnu.org/software/wget/)
- [**Slurm Workload Manager**](https://slurm.schedmd.com/download.html)
- 2 GB of disk space

The dependencies in bold are not installed in the steps above. `wget` is necessary to download hg38 reference files. Slurm is not necessary for running the sequential motif discovery tool.

# **Motif Discovery**

This tool can perform de novo regulatory motif discovery on hg38 fasta sequences, BED coordinates, or marker genes.

![Image of tool](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/motif_discovery.jpg)

## Parameters

### Required Argument

```
-i --input : Input data [TYPE : FASTA/BED/GENES]
```

The file extensions for the input passed to `-i` must be specified as follows:

TYPE | File Extension
------------ | -------------
FASTA | .fa
BED | .bed
GENES | .txt

The input data must correspond to at least two sequences. If the width is specified, it cannot be greater than the length of any sequence. BED coordinates must be mapped to the hg38 reference genome. The GENES file type is a single column of official human gene names crosslisted with GENCODE Human Release 33.

### Optional Arguments

```
-w --width : Integer length of motif.
-s --sequential : Run the motif discovery algorithm sequentially.
-p --parallel : Run the motif discovery algorithm in parallel. Requires Slurm Workload Manager.  
-h --help
```

If the `-w` flag is not provided, the program will infer the motif width.

If the `-s` or `-p` flags are not specified, the program will automatically choose the mode that gives better performance.

## Usage

```shell
cd motif_discovery/scripts
./find_motif.sh -i [/path/to/input]
```
## Output

The program will output the motif logo and a plot of the motif entropy over sampling iterations. The motif logo represents the alignment with maximimum entropy.

![Image of motif](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/example_motif.png)

![Image of entropy](https://github.com/chaitanyasrinivasan/motif_discovery/blob/master/images/example_entropy.png)
