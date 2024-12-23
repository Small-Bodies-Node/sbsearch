#!/bin/bash

set -e
pushd .

[[ -z $S2PREFIX ]] && echo "Requires env variable S2PREFIX set to desired installation prefix" && exit 1
[[ -z "$PYTHON_ROOT" ]] && PYTHON_ROOT=`python3 -c "import sys; print(sys.exec_prefix)"`

# checkout code, S2 and Abseil
test ! -e build && mkdir build
cd build

# install abseil from source, as per s2geometry's readme
# (it must be configured with -DCMAKE_POSITION_INDEPENDENT_CODE=ON)
# We are using C++11, and the last version to support it is LTS 20220623.1
if [ ! -e $S2PREFIX/lib/libabsl_base.a ]; then
    ABSEIL_TAG=20220623.1
    test ! -e ${ABSEIL_TAG}.tar.gz && wget https://github.com/abseil/abseil-cpp/archive/refs/tags/${ABSEIL_TAG}.tar.gz
    test ! -e abseil-cpp-${ABSEIL_TAG} && tar xzf ${ABSEIL_TAG}.tar.gz

    cd abseil-cpp-${ABSEIL_TAG}
    mkdir -p build
    cd build
    cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=$S2PREFIX -DABSL_ENABLE_INSTALL=ON -DABSL_PROPAGATE_CXX_STD=ON ..
    make -j $(nproc)
    make install
    cd ../..
fi

if [ ! -e $S2PREFIX/lib/libs2.so ]; then
    S2_TAG=v0.10.0
    test ! -e ${S2_TAG}.tar.gz && wget https://github.com/google/s2geometry/archive/refs/tags/${S2_TAG}.tar.gz
    test ! -e s2geometry-${S2_TAG} && tar xzf ${S2_TAG}.tar.gz

    export LDFLAGS="-L$S2PREFIX/lib -Wl,-rpath=$S2PREFIX/lib"
    export CXXFLAGS="-I$S2PREFIX/include"
    cd s2geometry-${S2_TAG:1}
    mkdir -p build
    cd build
    cmake -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=${S2PREFIX} -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=${S2PREFIX} -Wno-dev -DPython3_FIND_STRATEGY=LOCATION -DPython3_ROOT_DIR=${PYTHON_ROOT} ..
    make -j $(nproc)
    make install
fi

popd