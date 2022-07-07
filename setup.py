from distutils.core import setup, Extension
print("Running setup.py...")
lrst_mod = Extension("_lrst",
                     sources=["src/lrst_wrap.cxx", "src/Python_to_CXX_class.cpp", "src/CXX_to_python_class.cpp", "src/LRST_function.cpp","src/BAM_module.cpp", "src/BamReader.cpp", "src/F5_module.cpp", "src/kseq.cpp", "src/FASTQ_module.cpp", "src/ComFunction.cpp", "src/FASTA_module.cpp", "src/ComStruct.cpp" ],
                     language='c++',
                     extra_compile_args=['-std=c++14'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts"],
                     include_dirs=['src'])

setup(name = "lrst",
      version = '0.1',
      author      = "Qian Liu",
      description = """Long-read statistics""",
      ext_modules=[lrst_mod],
      py_modules=['lrst'])
