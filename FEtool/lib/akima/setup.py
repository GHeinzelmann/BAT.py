# setup.py
# Usage: ``python setup.py build_ext --inplace``
from distutils.core import setup, Extension
import numpy
setup(name='_akima', ext_modules=[Extension('_akima', ['akima.c'], include_dirs=[numpy.get_include()])])
