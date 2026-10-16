#!/usr/bin/env bash
set -euo pipefail
fortsym_source=${1:?Usage: regenerate_hpd.sh FORTSYM_SOURCE [GEN_BUILD [FORTSYM_BUILD]]}
generator_build=${2:-build/codegen-hpd}
fortsym_build=${3:-"$fortsym_source/build"}
expected_revision=cf7efc3864f6de26024c940bcf2ec96b3698b7af
test "$(git -C "$fortsym_source" rev-parse HEAD)" = "$expected_revision"
mkdir -p "$generator_build" src/generated
"${FC:-gfortran}" -I "$fortsym_build/include" \
    gen/gen_hpd_leaves.f90 -L "$fortsym_build/lib" \
    -Wl,-rpath,"$fortsym_build/lib" -lfortsym \
    -o "$generator_build/gen_hpd_leaves"
"$generator_build/gen_hpd_leaves"
