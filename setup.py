"""
setup.py:
Compile the module and its dependencies.
"""

import os
import glob
import subprocess
import setuptools
from setuptools import setup, Extension

print("Running setup.py...")

def get_git_version():
    """Get version from git tag, fallback to default if not available."""
    try:
        # Get the latest git tag
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, check=True)
        version = result.stdout.strip()
        # Remove 'v' prefix if present
        if version.startswith('v'):
            version = version[1:]
        return version
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Fallback to default version if git is not available or no tags exist
        print("Warning: Could not get version from git tag, using default version")
        return '1.5.0'
    
# Set the version
version = get_git_version()
print(f"Setting version to: {version}")

# Get the project dependencies
src_files = []
project_dir = 'src/'
project_src_files = []
project_src_files.extend(glob.glob(project_dir + '*.cpp'))
project_headers = glob.glob('include/*.h')

# Set up the extension
include_dirs = ['include']
lrst_mod = Extension("_lrst",
                     sources=project_src_files,
                     language='c++',
                     extra_compile_args=['-std=c++11'],
                     libraries=["hts", "hdf5_cpp", "hdf5", "hdf5_hl_cpp", "hdf5_hl"],
                     include_dirs=include_dirs,
                     depends=project_headers, )

# Set up the module
setup(name="longreadsum",
      version=version,
      author="WGLab",
      description="""A fast and flexible QC tool for long read sequencing data""",
      ext_modules=[lrst_mod],
      script_args=['build_ext', '--build-lib', 'lib'],
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
