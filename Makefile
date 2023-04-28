INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
LIB_DIR := $(CURDIR)/lib

all:
    # Generate the SWIG Python/C++ wrappers
	swig -c++ -python -outdir $(LIB_DIR) -I$(INCL_DIR) -o $(SRC_DIR)/lrst_wrap.cpp $(SRC_DIR)/lrst.i

	# Compile the C++ shared libraries into lib/
	python setup.py build_ext --build-lib $(LIB_DIR)
