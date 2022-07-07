SRC_DIR := src

all:
    # Generate the SWIG python <-> C++ wrappers (Both ways)
	swig -c++ -python -I$(SRC_DIR) lrst.i

	# Compile the C++ shared libraries
	python setup.py build_ext --inplace
