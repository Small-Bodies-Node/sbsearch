#!/bin/bash -e

set -e

PYTHON=python3
DEV=0

help() { echo "Usage: $0 [--python=<python3>] [--dev]" 1>&2; exit 1; }

PARSED=$(getopt -o "hp:d" --long "help,python:,dev" -- "$@")
if [ $? -ne 0 ]; then
        echo 'Terminating...' >&2
        exit 1
fi
eval set -- "$PARSED"
unset PARSED

while true; do
    case "$1" in
        "-h"|"--help")
            help
            exit
        ;;
        "-p"|"--python")
            PYTHON=$2
            shift 2
            continue
        ;;
        "-d"|"--dev")
            DEV=1
            shift
            continue
        ;;
        "--")
            shift
            break
        ;;
        *)
            echo "Internal error" >&2
            exit 1
        ;;
    esac
done

test ! -e .venv && $PYTHON -m venv .venv --prompt='sbsearch-v3'
source .venv/bin/activate
python3 -m pip install -q -U pip setuptools wheel

bash build_s2.sh ${VIRTUAL_ENV}

export LDFLAGS="-L${VIRTUAL_ENV}/lib -Wl,-rpath=${VIRTUAL_ENV}/lib"
export CXXFLAGS="-I${VIRTUAL_ENV}/include"
if [ $DEV -eq 0 ]; then
    python3 -m pip install -e .[recommended]
else
    python3 -m pip install -e .[recommended,test,docs]
fi

pushd .

cd build
test ! -e sbsearch && mkdir sbsearch
cd sbsearch

cmake ../.. -DCMAKE_PREFIX_PATH=${VIRTUAL_ENV} -DCMAKE_INSTALL_PREFIX=${VIRTUAL_ENV}
make -j$(nproc)
make install

popd

