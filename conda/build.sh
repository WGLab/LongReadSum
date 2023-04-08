#!/bin/bash

#export CFLAGS="$CFLAGS -I$PREFIX/include"
#export LDFLAGS=${LDFLAGS}" -I${PREFIX}/include -L${PREFIX}/lib"
#export CPATH=${CPATH}" -I${PREFIX}/include -L${PREFIX}/lib"
#export C_INCLUDE_PATH=${C_INCLUDE_PATH}:${PREFIX}/include
#export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:${PREFIX}/include
#export LIBRARY_PATH=${LIBRARY_PATH}:${PREFIX}/lib

echo "Calling entry_point.py and printing outputs"
#export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

echo "Lib contents: "
ls -l "${PREFIX}"/lib

# Add the library path to the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib

#export LDFLAGS="$LDFLAGS -Wl,-rpath,$PREFIX/lib"

# Generate the SWIG files
swig -c++ -python -outdir "${SRC_DIR}"/lib -I"${SRC_DIR}"/include -I"${PREFIX}"/include -o "${SRC_DIR}"/src/lrst_wrap.cpp "${SRC_DIR}"/src/lrst.i

#echo "HTS location: "
#ldconfig -p | grep libhts
#
#echo "All library files in PREFIX/lib containing the word hdf5: "
#ls -l "${PREFIX}"/lib | grep hdf5

# Generate the SWIG library
#$PYTHON setup.py install
# Run the setup.py script linking all libraries
$PYTHON setup.py -I"${PREFIX}"/include -L"${PREFIX}"/lib install

# Copy the entry point to the bin directory
cp "${SRC_DIR}"/longreadsum "${PREFIX}"/bin
#chmod +x "$PREFIX"/bin/longreadsum

# Create the src directory
mkdir -p "${PREFIX}"/src  # Create the longreadsum src directory

# Copy the lib files to the src directory
cp -r "${SRC_DIR}"/src/*.py "${PREFIX}"/src
#mkdir -p "${PREFIX}"/lib/longreadsum  # Create the longreadsum lib directory
#cp -r "${SRC_DIR}"/src/*.py "${PREFIX}"/lib/longreadsum

# Copy the SWIG generated library to the lib directory
#cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib/longreadsum
#cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/lib
cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/src
#cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib
cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib
#cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/lib

#echo "Bin contents: "
#ls -l "${PREFIX}"/bin
#
#echo "Lib contents: "
#ls -l "${PREFIX}"/lib
#
#echo "Src contents: "
#ls -l "${PREFIX}"/src
#ls -l "${PREFIX}"/lib/longreadsum
