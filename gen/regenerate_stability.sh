#!/usr/bin/env bash
set -euo pipefail
sym_source=${1:?Usage: regenerate_stability.sh FORTSYM_SOURCE [FORTSYM_BUILD [GEN_BUILD]]}
sym_build=${2:-$sym_source/build}
gen_build=${3:-build/stability-codegen}
test "$(git -C "$sym_source" rev-parse HEAD)" = cf7efc3864f6de26024c940bcf2ec96b3698b7af
mkdir -p "$gen_build"
"${FC:-gfortran}" -I "$sym_build/include" gen/gen_stability_oracle.f90 \
    -L "$sym_build/lib" -Wl,-rpath,"$sym_build/lib" -lfortsym \
    -o "$gen_build/gen_stability_oracle"
"$gen_build/gen_stability_oracle"
