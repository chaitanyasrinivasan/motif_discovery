from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=[
        Extension("parallel_cygibbs", ["parallel_cygibbs.c"],
                  include_dirs=[numpy.get_include()]),
    ],
)

# Or, if you use cythonize() to make the ext_modules list,
# include_dirs can be passed to setup()

setup(
    ext_modules=cythonize("parallel_cygibbs.pyx"),
    include_dirs=[numpy.get_include()]
)    
