#!/bin/bash

export CFLAGS="$CFLAGS -I$PREFIX/include"
export LDFLAGS=${LDFLAGS}" -I${PREFIX}/include -L${PREFIX}/lib"
export CPATH=${CPATH}" -I${PREFIX}/include -L${PREFIX}/lib"
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:${PREFIX}/include
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:${PREFIX}/include
export LIBRARY_PATH=${LIBRARY_PATH}:${PREFIX}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib
export LDFLAGS="$LDFLAGS -Wl,-rpath,$PREFIX/lib"

# Generate the SWIG files
swig -c++ -python -outdir "${SRC_DIR}"/lib -I"${SRC_DIR}"/include -o "${SRC_DIR}"/src/lrst_wrap.cpp "${SRC_DIR}"/src/lrst.i

# Generate the SWIG library
python setup.py install --single-version-externally-managed --record=record.txt  # SO should generate in $PREFIX/lib

# Copy the entry point to the bin directory
cp "${SRC_DIR}"/longreadsum "${PREFIX}"/bin
chmod +x "$PREFIX"/bin/longreadsum

# Copy the lib files to the lib directory
mkdir -p "${PREFIX}"/lib/longreadsum  # Create the longreadsum lib directory
cp -r "${SRC_DIR}"/src/*.py "${PREFIX}"/lib/longreadsum

# Copy the SWIG generated library to the lib directory
#cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib/longreadsum
cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/lib/longreadsum
cp -r "${SRC_DIR}"/lib/*.so "${PREFIX}"/lib
#cp -r "${SRC_DIR}"/lib/*.py "${PREFIX}"/lib

echo "Lib contents: "
ls -l "${PREFIX}"/lib

echo "Longreadsum ib contents: "
ls -l "${PREFIX}"/lib/longreadsum
