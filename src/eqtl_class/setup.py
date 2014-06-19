from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "cpp_eqtl.pyx",                 # our Cython source
           sources=["eqtl_m.cpp"],  # additional source file(s)
           language="c++",             # generate C++ code
      ))