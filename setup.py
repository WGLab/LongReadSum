"""
setup.py:
Compile the module and its dependencies.
"""

import os
import glob
from distutils.core import setup, Extension

print("Running setup.py...")

# Get the project dependencies
src_files = []
project_dir = 'src/'
project_src_files = []
project_src_files.extend(glob.glob(project_dir + '*.cpp'))
project_src_files.extend(glob.glob(project_dir + '*.cxx'))

# Set up the extension
# include_dirs = ['include/', 'lib/hdf/HDF_Group/HDF5/1.12.2/include/']
include_dirs = ['include/', 'include/hdf5-1_12_1/']
lrst_mod = Extension("_lrst",
                     sources=project_src_files,
                     language='c++',
                     extra_compile_args=['-std=c++14'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts", "hdf5_cpp", "hdf5", "hdf5_hl_cpp", "hdf5_hl"],
                     library_dirs=['lib/hdf5-1_12_1/'],
                     include_dirs=include_dirs)

#library_dirs=['lib/hdf/HDF_Group/HDF5/1.12.2/lib/']

# Set up the module
setup(name = "lrst",
      version = '0.1',
      author      = "Qian Liu",
      description = """Long-read statistics""",
      ext_modules=[lrst_mod],
      py_modules=['lrst'],
      include_dirs=include_dirs,
      test_suite='tests')
