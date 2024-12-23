INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
LIB_DIR := $(CURDIR)/lib

# Set the library paths for the compiler
#LIBRARY_PATHS := -L$(LIB_DIR) -L/usr/share/miniconda/envs/longreadsum/lib
#INCLUDE_PATHS := -I$(INCL_DIR) -I/usr/share/miniconda/envs/longreadsum/include
CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)
LIBRARY_PATHS := -L$(LIB_DIR) -L$(CONDA_PREFIX)/lib
INCLUDE_PATHS := -I$(INCL_DIR) -I$(CONDA_PREFIX)/include

# All targets
all: swig_build compile

# Generate the SWIG Python/C++ wrappers
swig_build:
	swig -c++ -python -outdir $(LIB_DIR) -I$(INCL_DIR) -o $(SRC_DIR)/lrst_wrap.cpp $(SRC_DIR)/lrst.i

# Compile the C++ shared libraries into lib/
compile:
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(CONDA_PREFIX)/lib \
	CXXFLAGS="$(INCLUDE_PATHS)" LDFLAGS="$(LIBRARY_PATHS)" python3 setup.py build_ext --build-lib $(LIB_DIR)

# LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/usr/share/miniconda/envs/longreadsum/lib \
