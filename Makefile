SRC_DIR := $(CURDIR)/src/longreadsum
#SRC_DIR := src
# $(info $(SRC_DIR))

all:
    # Generate the SWIG python <-> C++ wrappers (Both ways)
    # The Python wrapper is stored in the main directory, while the C++ wrapper code is stored in /src
	swig -c++ -python -outdir $(SRC_DIR) -I$(SRC_DIR) -o $(SRC_DIR)/lrst_wrap.cxx lrst.i

	# Compile the C++ shared libraries
	# The C++ wrapper binary _lrst is generated in the main directory
	python setup.py build_ext --inplace
