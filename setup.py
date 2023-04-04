"""
setup.py:
Compile the module and its dependencies.
"""

import os
import glob
import setuptools
from setuptools import setup, Extension

print("Running setup.py...")

# Get the project dependencies
src_files = []
project_dir = 'longreadsum/'
project_src_files = []
project_src_files.extend(glob.glob(project_dir + '*.cpp'))
project_src_files.extend(glob.glob(project_dir + '*.cxx'))

project_headers = glob.glob('include/*.h')
print("header files: ")
print(project_headers)

# Set up the extension
include_dirs = ['include']
print("CWDS: ", os.getcwd())
print("INCLUDE DIRS: ")
print(include_dirs)
lrst_mod = Extension("_lrst",
                     sources=project_src_files,
                     language='c++',
                     extra_compile_args=['-std=c++11'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts", "hdf5_cpp", "hdf5", "hdf5_hl_cpp", "hdf5_hl"],
                     include_dirs=include_dirs,
                     depends=project_headers, )

# Set up the module
setup(name="longreadsum",
      version='1.0.1',
      author="WGLab",
      description="""A fast and flexible QC tool for long read sequencing data""",
      ext_modules=[lrst_mod],
      script_args=['build_ext', '--inplace', '--build-lib', 'lib'],
      py_modules=['lrst'],
      packages=setuptools.find_packages(),
      headers=project_headers,
      include_dirs=include_dirs,
      test_suite='tests',
      entry_points={
          'console_scripts': [
              'longreadsum = longreadsum.__main__:main'
          ]
      },
      )
