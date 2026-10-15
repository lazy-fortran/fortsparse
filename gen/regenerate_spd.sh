#!/usr/bin/env bash
set -euo pipefail
fortsym_source=${1:?Usage: gen/regenerate_spd.sh /path/to/fortsym}
generator_build=${2:-build/codegen}
expected_revision=283fbb8e3b6b14e269fa404284d8a0c47a9e6227
test "$(git -C "$fortsym_source" rev-parse HEAD)" = "$expected_revision"
mkdir -p "$generator_build" src/generated
"${FC:-gfortran}" -I "$fortsym_source/build/include" \
    gen/gen_spd_leaves.f90 -L "$fortsym_source/build/lib" \
    -Wl,-rpath,"$fortsym_source/build/lib" -lfortsym \
    -o "$generator_build/gen_spd_leaves"
"$generator_build/gen_spd_leaves"
