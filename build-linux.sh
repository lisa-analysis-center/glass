#!/usr/bin/env bash
# Linux (Ubuntu/Debian) build. Requires:
#   sudo apt install build-essential gfortran cmake mpi-default-dev \
#       libopenblas-dev liblapacke-dev libhdf5-dev
set -euo pipefail
cd "$(dirname "$0")"

# Locate the OpenBLAS CMake config. Debian/Ubuntu ship it under a variant
# directory (openblas-openmp, openblas-pthread, openblas-serial, ...) depending
# on which libopenblas*-dev package is installed, so don't hardcode the variant.
openblas_config="$(find /usr/lib /usr/lib64 -path '*/cmake/openblas/OpenBLASConfig.cmake' 2>/dev/null | head -n1)"
if [ -z "$openblas_config" ]; then
  echo "error: could not find OpenBLASConfig.cmake; install libopenblas-dev" >&2
  exit 1
fi

cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$(dirname "$openblas_config")" \
  -DBLA_VENDOR=OpenBLAS

cmake --build build -j"$(nproc)"
