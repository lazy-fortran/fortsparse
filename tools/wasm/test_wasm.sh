#!/usr/bin/env bash
# Builds fortsparse and its test suite for WebAssembly with LFortran and
# emscripten and runs the tests under node, the way a browser runs them:
#   tools/wasm/test_wasm.sh [BUILD_DIR] [CTEST_ARGUMENT ...]
# Needs LFortran built with WITH_TARGET_WASM (LFORTRAN, default lfortran on
# PATH), an emscripten SDK (EMSDK_PATH) and node (NODE, default node on PATH).
# Without CTEST_ARGUMENT it runs the whole suite; -R PATTERN selects tests.
# The browser has no processes, threads or C solver libraries: the build has
# the native solvers only, without SuperLU, the UMFPACK helper, METIS, AMD
# and OpenMP.
set -euo pipefail
top=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out=${1:-$top/build-wasm}
shift || true
: "${EMSDK_PATH:?set EMSDK_PATH to the emscripten SDK}"
cmake -S "$top" -B "$out" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE="$top/tools/wasm/emscripten-lfortran.cmake" \
    -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
    -DFORTSPARSE_ENABLE_METIS=OFF -DFORTSPARSE_ENABLE_AMD=OFF \
    -DFORTSPARSE_ENABLE_OPENMP=OFF
cmake --build "$out"
ctest --test-dir "$out" --output-on-failure "$@"
