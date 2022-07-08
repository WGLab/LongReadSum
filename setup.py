"""
setup.py:
Compile the module and its dependencies.
"""

import os
from distutils.core import setup, Extension

print("Running setup.py...")

# Get the C++ dependencies
project_dir = 'src/longreadsum'
project_src_files = []
for file in os.listdir(project_dir):
    if file.endswith(".cxx") or file.endswith(".cpp"):
        src_filepath = os.path.join(project_dir, file)
        project_src_files.append(src_filepath)

lrst_mod = Extension("_lrst",
                     sources=project_src_files,
                     language='c++',
                     extra_compile_args=['-std=c++14'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts"],
                     include_dirs=[project_dir])

setup(name = "lrst",
      version = '0.1',
      author      = "Qian Liu",
      description = """Long-read statistics""",
      ext_modules=[lrst_mod],
      py_modules=['lrst'],
      test_suite='tests')
