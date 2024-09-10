#!/bin/bash

set -e

if [ -u $1 ]; then
    echo "Usage: _build_s2 <install_prefix>"
    exit 1
fi

PREFIX=$(readlink -f $1)

pushd .

[[ -z "$PYTHON_ROOT" ]] && PYTHON_ROOT=`python3 -c "import sys; print(sys.exec_prefix)"`

test ! -e build && mkdir build
cd build

ABSEIL_TAG=20220623.1
if [ ! -e $PREFIX/lib/libabsl_base.a ]; then
    test ! -e ${ABSEIL_TAG}.tar.gz && wget https://github.com/abseil/abseil-cpp/archive/refs/tags/${ABSEIL_TAG}.tar.gz
    test ! -e abseil-cpp-${ABSEIL_TAG} && tar xzf ${ABSEIL_TAG}.tar.gz

    cd abseil-cpp-${ABSEIL_TAG}
    mkdir -p build
    cd build
    cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX=$PREFIX -DABSL_ENABLE_INSTALL=ON -DABSL_PROPAGATE_CXX_STD=ON ..
    make -j $(nproc)
    make install
    cd ../..
fi

S2_TAG=v0.10.0
if [ ! -e $PREFIX/lib/libs2.so ]; then
    test ! -e ${S2_TAG}.tar.gz && wget https://github.com/google/s2geometry/archive/refs/tags/${S2_TAG}.tar.gz
    test ! -e s2geometry-${S2_TAG} && tar xzf ${S2_TAG}.tar.gz

    export LDFLAGS="-L$PREFIX/lib -Wl,-rpath=$PREFIX/lib"
    export CXXFLAGS="-I$PREFIX/include"
    cd s2geometry-${S2_TAG:1}
    mkdir -p build
    cd build
    cmake -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=${PREFIX} -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX=${PREFIX} -Wno-dev -DPython3_FIND_STRATEGY=LOCATION -DPython3_ROOT_DIR=${PYTHON_ROOT} ..
    make -j $(nproc)
    make install
fi

popd