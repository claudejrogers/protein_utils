from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name='protein_utils',
    version='0.1',
    description='Utilities for protein coordinate files',
    author='Claude Rogers',
    author_email='cjrogers@caltech.edu',
    url='https://github.com/claudejrogers/protein_utils',
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("protutils.cealign.cealign",
                  ["protutils/cealign/cealign.pyx"],
                  include_dirs=[numpy.get_include()],
                  libraries=["m"])
    ]
)
