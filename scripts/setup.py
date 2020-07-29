from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

'''
Chaitanya Srinivasan

This script is called during installation to produce Cython executable scripts
usage: python setup.py build_ext --inplace
'''

#cythonize cygibbs.pyx and align.pyx 
setup(
    ext_modules=cythonize("*.pyx"),
    include_dirs=[numpy.get_include()]
)
