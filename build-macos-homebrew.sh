#!/usr/bin/env bash
# macOS (Homebrew) build. Requires:
#   brew install cmake llvm open-mpi openblas libomp hdf5
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v brew >/dev/null 2>&1; then
  echo "error: Homebrew not found; install it from https://brew.sh" >&2
  exit 1
fi

for formula in llvm open-mpi openblas libomp hdf5; do
  if ! brew --prefix --installed "$formula" >/dev/null 2>&1; then
    echo "error: missing Homebrew formula '$formula'; run: brew install $formula" >&2
    exit 1
  fi
done

llvm_prefix="$(brew --prefix llvm)"
openblas_prefix="$(brew --prefix openblas)"
libomp_prefix="$(brew --prefix libomp)"
hdf5_prefix="$(brew --prefix hdf5)"
openmpi_prefix="$(brew --prefix open-mpi)"

# Apple clang has no OpenMP support, so point Open MPI's mpicc wrapper at
# Homebrew's llvm clang and hand CMake the libomp location explicitly.
export OMPI_CC="$llvm_prefix/bin/clang"

CFLAGS="-I$openblas_prefix/include -I$libomp_prefix/include -L$libomp_prefix/lib" \
cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
  -DCMAKE_BUILD_TYPE=Release\
  -DCMAKE_PREFIX_PATH="$openblas_prefix;$hdf5_prefix;$openmpi_prefix" \
  -DBLA_VENDOR=OpenBLAS \
  -DOpenMP_C_FLAGS="-fopenmp -L$libomp_prefix/lib" \
  -DOpenMP_C_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY="$libomp_prefix/lib/libomp.dylib"

cmake --build build -j"$(sysctl -n hw.ncpu)"
