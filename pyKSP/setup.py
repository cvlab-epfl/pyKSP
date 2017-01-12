from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
 
setup(
  name = 'Demos',
  ext_modules=[ 
    Extension("ksp", 
              sources=["kgraph.pyx", "ksp_graph.cpp","ksp_computer.cpp"], # Note, you can link against a c++ library instead of including the source
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
 
)