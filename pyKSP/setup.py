from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

"""
Build ksp and import from anywhere:
    cd pyKSP
	python setup.py build
    python setup.py install

Build ksp and import from this directory
    cd pyKSP
	python setup.py build_ext --inplace
"""

setup(
    name = 'pyKSP',
    ext_modules=cythonize(Extension("ksp", 
        sources=["kgraph.pyx", "ksp_graph.cpp","ksp_computer.cpp"], # Note, you can link against a c++ library instead of including the source
        language="c++",
    )),
)