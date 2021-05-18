from distutils.core import setup, Extension

lrst_mod = Extension("_lrst",
                     sources=["lrst_wrap.cxx", "BAM_module.cpp", "BamReader.cpp", "ComFunction.cpp", "ComStruct.cpp", "CXX_to_python_class.cpp", "LRST_function.cpp", "Python_to_CXX_class.cpp" ], 
                     language='c++',
                     extra_compile_args=['-std=c++14'],
                     libraries=["hts","rt", "pthread", "z", "dl", "m"])

setup(name = "lrst", 
      version = '0.1',
      author      = "Qian Liu",
      description = """Long-read statistics""",
      ext_modules=[lrst_mod], 
      py_modules=['lrst'])
