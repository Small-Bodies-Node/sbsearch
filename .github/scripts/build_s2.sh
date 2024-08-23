#!/bin/bash

# build and install s2 in a Github runner

set -e

S2PREFIX=/usr/local
ABSEIL_TAG=20220623.1
S2_TAG=v0.10.0

mkdir build
cd build

wget https://github.com/abseil/abseil-cpp/archive/refs/tags/${ABSEIL_TAG}.tar.gz
tar xzf ${ABSEIL_TAG}.tar.gz

cd abseil-cpp-${ABSEIL_TAG}
mkdir -p build
cd build
cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=$S2PREFIX -DABSL_ENABLE_INSTALL=ON -DABSL_PROPAGATE_CXX_STD=ON ..
make -j $(nproc)
sudo make install
cd ../..

wget https://github.com/google/s2geometry/archive/refs/tags/${S2_TAG}.tar.gz
tar xzf ${S2_TAG}.tar.gz

export LDFLAGS="-L$S2PREFIX/lib -Wl,-rpath=$S2PREFIX/lib"
export CXXFLAGS="-I$S2PREFIX/include"
cd s2geometry-${S2_TAG:1}
mkdir -p build
cd build
cmake -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=$S2PREFIX -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=$S2PREFIX -Wno-dev ..
make -j $(nproc)
sudo make install
