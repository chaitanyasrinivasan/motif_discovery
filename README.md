# Installation

1. Download anaconda from https://docs.anaconda.com/anaconda/ to ensure your system has the compatible version of Python along with the necessary packages.

2. Install Cython with the following command:

```shell
   conda install -c anaconda cython
```

3. Install seqlogo with the following commands:
```shell
  git clone https://github.com/betteridiot/seqlogo.git
  cd seqLogo
  python setup.py install
  conda env create -f environment.yml
```

# **Gibbs Sampler Motif Finder**



