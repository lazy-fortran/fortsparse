#!/usr/bin/env bash
set -euo pipefail
sym_source=${1:?Usage: regenerate_stability.sh FORTSYM_SOURCE [FORTSYM_BUILD [GEN_BUILD]]}
sym_build=${2:-$sym_source/build}
gen_build=${3:-build/stability-codegen}
test "$(git -C "$sym_source" rev-parse HEAD)" = a24e64dbd0b8dba0f7d7e686b1f71c0331b5c336
mkdir -p "$gen_build"
"${FC:-gfortran}" -I "$sym_build/include" gen/gen_stability_oracle.f90 \
    -L "$sym_build/lib" -Wl,-rpath,"$sym_build/lib" -lfortsym \
    -o "$gen_build/gen_stability_oracle"
"$gen_build/gen_stability_oracle"
