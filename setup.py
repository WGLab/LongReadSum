from distutils.core import setup, Extension

lrst_mod = Extension("_lrst",
                     sources=["lrst_wrap.cxx", "Python_to_CXX_class.cpp", "CXX_to_python_class.cpp", "LRST_function.cpp","BAM_module.cpp", "BamReader.cpp", "ComFunction.cpp", "ComStruct.cpp" ],
                     language='c++',
                     extra_compile_args=['-std=c++14'],
                     libraries=["rt", "pthread", "z", "dl", "m", "hts"])

setup(name = "lrst",
      version = '0.1',
      author      = "Qian Liu",
      description = """Long-read statistics""",
      ext_modules=[lrst_mod],
      py_modules=['lrst'])

