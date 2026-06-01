#!/usr/bin/env bash
# Linux (Ubuntu/Debian) build. Requires:
#   sudo apt install build-essential gfortran cmake mpi-default-dev \
#       libopenblas-dev liblapacke-dev libhdf5-dev
set -euo pipefail
cd "$(dirname "$0")"

cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="/usr/lib/x86_64-linux-gnu/openblas-openmp/cmake/openblas" \
  -DBLA_VENDOR=OpenBLAS

cmake --build build -j"$(nproc)"
