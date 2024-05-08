#!/bin/bash

# Add the library path to the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib

# Generate the SWIG files
swig -c++ -python -outdir "${SRC_DIR}"/lib -I"${SRC_DIR}"/include -I"${PREFIX}"/include -o "${SRC_DIR}"/src/lrst_wrap.cpp "${SRC_DIR}"/src/lrst.i

# Generate the shared library
$PYTHON setup.py -I"${PREFIX}"/include -L"${PREFIX}"/lib install

# Create the src directory
mkdir -p "${PREFIX}"/src

# Copy python source files to the bin directory
cp -r "${SRC_DIR}"/src/*.py "${PREFIX}"/bin

# Copy the SWIG generated library to the lib directory
cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/lib
cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib
# cp -r "${SRC_DIR}"/lib/lrst.py "${PREFIX}"/lib
# cp -r "${SRC_DIR}"/lib/_lrst.*.so "${PREFIX}"/lib

# Download the pod5 library for linux x64
# wget https://github.com/nanoporetech/pod5-file-format/releases/download/0.3.10/lib_pod5-0.3.10-linux-x64.tar.gz

# Extract the pod5 library
# tar -xzvf lib_pod5-0.3.10-linux-x64.tar.gz

# echo "Showing the contents of the extracted pod5 library"
# ls lib

# Copy the pod5 library to the lib directory
# cp lib/libpod5_format.so "${PREFIX}"/lib

# Copy the pod5 header files to the include directory under pod5_format
# cp -r include/pod5_format "${PREFIX}"/include

# Check that the pod5_format folder is in the include directory
# echo "Checking that the pod5_format folder is in the include directory (ls ${PREFIX}/include/pod5_format)"
# ls "${PREFIX}"/include/pod5_format
