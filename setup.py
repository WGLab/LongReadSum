"""
setup.py:
Compile the module and its dependencies.
"""

import os
import glob
import setuptools
# from distutils.core import setup, Extension
from setuptools import setup, Extension

print("Running setup.py...")

# Get the project dependencies
src_files = []
project_dir = 'src/'
project_src_files = []
project_src_files.extend(glob.glob(project_dir + '*.cpp'))
project_src_files.extend(glob.glob(project_dir + '*.cxx'))

project_headers = glob.glob('include/*.h')
print("header files: ")
print(project_headers)

# Set up the extension
# include_dirs = ['include/', 'lib/hdf/HDF_Group/HDF5/1.12.2/include/']
# include_dirs = ['/home/perdomoj/github/LongReadSum/include/', 'include/hdf5-1_12_1/']
include_dirs = ['include', 'include/hdf5-1_12_1/']
print("CWDS: ", os.getcwd())
print("INCLUDE DIRS: ")
print(include_dirs)
lrst_mod = Extension("_lrst",
                     sources=project_src_files,
                     language='c++',
                     extra_compile_args=['-std=c++11'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts", "hdf5_cpp", "hdf5", "hdf5_hl_cpp", "hdf5_hl"],
                     library_dirs=['lib/hdf5-1_12_1/'],
                     include_dirs=include_dirs,
                     depends=project_headers)

#library_dirs=['lib/hdf/HDF_Group/HDF5/1.12.2/lib/']

# Set up the module
setup(name = "longreadsum",
      version = '1.0.1',
      author      = "WGLab",
      description = """A fast and flexible QC tool for long read sequencing data""",
      ext_modules=[lrst_mod],
      py_modules=['lrst'],
      packages=setuptools.find_packages(),
      headers=project_headers,
      include_dirs=include_dirs,
      test_suite='tests',
      entry_points={
          'console_scripts': [
              'longreadsum = src.__main__:main'
          ]
      },
      )


# python3 -m build
# pip install --editable .
# longreadsum
# python3 -m twine upload --repository testpypi dist/*

# docker pull quay.io/pypa/manylinux1_x86_64
# docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/repair_wheel.sh /io/dist/longreadsum-1.0.1-cp39-cp39-linux_x86_64.whl

# TODO: Try the above steps using the Python 3.6 environment to compile the wheel for Python 3.6

# Auditwheel error: cannot repair "dist/longreadsum-1.0.1-cp39-cp39-linux_x86_64.whl" to "manylinux_2_5_x86_64" ABI because of the presence of too-recent versioned symbols.
# https://github.com/pypa/auditwheel
