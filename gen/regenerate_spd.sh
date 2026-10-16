#!/usr/bin/env bash
set -euo pipefail
fortsym_source=${1:?Usage: regenerate_spd.sh FORTSYM_SOURCE [GEN_BUILD [FORTSYM_BUILD]]}
generator_build=${2:-build/codegen}
fortsym_build=${3:-"$fortsym_source/build"}
expected_revision=ded561c5e4119aa176fec1541afac7c106f30fa4
test "$(git -C "$fortsym_source" rev-parse HEAD)" = "$expected_revision"
mkdir -p "$generator_build" src/generated
"${FC:-gfortran}" -I "$fortsym_build/include" \
    gen/gen_spd_leaves.f90 -L "$fortsym_build/lib" \
    -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
    -o "$generator_build/gen_spd_leaves"
"$generator_build/gen_spd_leaves"
