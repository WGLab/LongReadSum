INCL_DIR := $(CURDIR)/include/longreadsum
SRC_DIR := $(CURDIR)/src/longreadsum
LIB_DIR := $(CURDIR)/lib

all:
	# Create the output directory
	mkdir -p $(LIB_DIR)

    # Generate the SWIG python <-> C++ wrappers (Both ways)
    # The Python wrapper is stored in the main directory, while the C++ wrapper code is stored in /src
	swig -c++ -python -outdir $(LIB_DIR) -I$(INCL_DIR) -o $(SRC_DIR)/lrst_wrap.cxx $(SRC_DIR)/lrst.i

	# Compile the C++ shared libraries
	# The C++ wrapper binary _lrst is generated in the /src directory
	python setup.py build_ext --build-lib $(LIB_DIR)
