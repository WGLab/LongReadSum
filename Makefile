INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
LIB_DIR := $(CURDIR)/lib

# All targets
all: swig_build compile

# Generate the SWIG Python/C++ wrappers
swig_build:
	swig -c++ -python -outdir $(LIB_DIR) -I$(INCL_DIR) -o $(SRC_DIR)/lrst_wrap.cpp $(SRC_DIR)/lrst.i

# Compile the C++ shared libraries into lib/
compile:
	python3 setup.py build_ext --build-lib $(LIB_DIR)
