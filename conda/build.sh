#!/bin/bash

# Add the library path to the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib

# Generate the SWIG files
swig -c++ -python -outdir "${SRC_DIR}"/lib -I"${SRC_DIR}"/include -I"${PREFIX}"/include -o "${SRC_DIR}"/src/lrst_wrap.cpp "${SRC_DIR}"/src/lrst.i

# Generate the shared library
$PYTHON setup.py -I"${PREFIX}"/include -L"${PREFIX}"/lib install

# Copy the entry point to the bin directory
cp "${SRC_DIR}"/longreadsum "${PREFIX}"/bin

# Create the src directory
mkdir -p "${PREFIX}"/src

# Copy the lib files to the src directory
cp -r "${SRC_DIR}"/src/*.py "${PREFIX}"/src

# Copy the SWIG generated library to the lib directory
cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/src
cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib
