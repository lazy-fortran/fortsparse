#!/usr/bin/env bash
# Optional local comparator provision; no downloads, source edits or system install.
set -euo pipefail
retained_source=${1:?Usage: provision_openblas.sh RETAINED_SUITESPARSE TASK_ROOT [OPENBLAS_LIB [OPENBLAS_INCLUDE]]}
task_root=${2:?Usage: provision_openblas.sh RETAINED_SUITESPARSE TASK_ROOT [OPENBLAS_LIB [OPENBLAS_INCLUDE]]}
openblas_library=${3:-/usr/lib/libopenblas.so}
openblas_include=${4:-/usr/include/openblas}
source_revision=34db4dfbe21c53cea7b0e3c5ee2da37bfaef369a
mkdir -p "$task_root"
task_root=$(cd "$task_root" && pwd)
retained_source=$(cd "$retained_source" && pwd)
test -f "$openblas_library"
test -f "$openblas_include/openblas_config.h"
git -C "$retained_source" cat-file -e "$source_revision^{commit}"
if test -d "$task_root/source"; then
    test "$(cat "$task_root/source-revision.txt")" = "$source_revision"
else
    git -C "$retained_source" archive --format=tar --output="$task_root/source.tar" "$source_revision"
    mkdir "$task_root/source"
    tar -xf "$task_root/source.tar" -C "$task_root/source"
    echo "$source_revision" > "$task_root/source-revision.txt"
fi
# Exact frozen committed source: uncommitted retained-clone edits are excluded.
cmake -S "$task_root/source" -B "$task_root/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$task_root/prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_RPATH="$task_root/prefix/lib" \
    -DSUITESPARSE_ENABLE_PROJECTS=umfpack -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF -DBUILD_TESTING=OFF -DSUITESPARSE_DEMOS=OFF \
    -DSUITESPARSE_USE_CUDA=OFF -DSUITESPARSE_USE_64BIT_BLAS=OFF \
    -DBLA_VENDOR=OpenBLAS -DBLA_SIZEOF_INTEGER=4 \
    -DBLAS_LIBRARIES="$openblas_library" -DLAPACK_LIBRARIES="$openblas_library" \
    -DBLAS_INCLUDE_DIRS="$openblas_include" -DLAPACK_INCLUDE_DIRS="$openblas_include"
cmake --build "$task_root/build" --parallel 2
cmake --install "$task_root/build"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 - "$retained_source" "$task_root" "$openblas_library" "$openblas_include" <<'PY'
import ctypes
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
retained, root, blas, include = map(Path, sys.argv[1:])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
config = ctypes.CDLL(str(blas))
config.openblas_get_config.restype = ctypes.c_char_p
blas_config = config.openblas_get_config().decode()
if "USE64BITINT" in blas_config:
    raise SystemExit("LP64 OpenBLAS required; runtime reports USE64BITINT")
cache = {}
for line in (root/"build/CMakeCache.txt").read_text().splitlines():
    if not line or line.startswith(("#", "//")) or "=" not in line:
        continue
    name, value = line.split("=", 1)
    cache[name.split(":", 1)[0]] = value
probe_source = root/"lp64_probe.c"
probe_source.write_text(r'''/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <cblas.h>
#include <umfpack.h>
#include <stdio.h>
#include <string.h>
_Static_assert(sizeof(blasint) == 4, "OpenBLAS header must be LP64");
_Static_assert(sizeof(SUITESPARSE_BLAS_INT) == 4, "SuiteSparse BLAS must be LP64");
_Static_assert(sizeof(int32_t) == 4, "UMFPACK di/zi indices must be I4");
int main(void) {
    int v[3];
    const double a = 3.0, b = 7.0;
    double c = 0.0;
    umfpack_version(v);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        1, 1, 1, 1.0, &a, 1, &b, 1, 0.0, &c, 1);
    if (v[0] != UMFPACK_MAIN_VERSION || v[1] != UMFPACK_SUB_VERSION ||
        v[2] != UMFPACK_SUBSUB_VERSION || c != 21.0 ||
        strstr(openblas_get_config(), "USE64BITINT")) return 1;
    printf("UMFPACK %d.%d.%d; blasint=%zu; SuiteSparse BLAS=%zu; DGEMM=%.0f\n",
        v[0], v[1], v[2], sizeof(blasint), sizeof(SUITESPARSE_BLAS_INT), c);
    return 0;
}
''')
probe_binary = root/"lp64_probe"
probe_command = [cache["CMAKE_C_COMPILER"], "-std=c11", "-Wall", "-Wextra", "-Werror",
    "-I"+str(include), "-I"+str(root/"prefix/include/suitesparse"), str(probe_source),
    str(root/"prefix/lib/libumfpack.so"), str(blas),
    "-Wl,-rpath,"+str(root/"prefix/lib"), "-o", str(probe_binary)]
subprocess.run(probe_command, check=True)
probe_output = subprocess.check_output([str(probe_binary)], text=True).strip()
paths = [root/"prefix/lib/libumfpack.so", root/"prefix/lib/libcholmod.so",
         root/"prefix/lib/libsuitesparseconfig.so"]
record = dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    source_revision=(root/"source-revision.txt").read_text().strip(),
    source_tree=subprocess.check_output(["git", "-C", str(retained), "rev-parse",
        "34db4dfbe21c53cea7b0e3c5ee2da37bfaef369a^{tree}"], text=True).strip(),
    source_archive_sha256=digest(root/"source.tar"),
    retained_clone_status_ignored=subprocess.check_output(
        ["git", "-C", str(retained), "status", "--porcelain"], text=True),
    source_policy="git archive of pinned commit; no retained-clone edits copied",
    package_version="UMFPACK 6.3.8; same version as inspected installed header",
    license="UMFPACK GPL-2.0-or-later; exact source Doc/License.txt retained",
    umfpack_header_sha256=digest(root/"prefix/include/suitesparse/umfpack.h"),
    umfpack_license_sha256=digest(root/"source/UMFPACK/Doc/License.txt"),
    source_license_paths=[str(p.relative_to(root)) for p in sorted((root/"source").glob("*/Doc/License.txt"))],
    openblas=dict(path=str(blas), resolved=str(blas.resolve()), sha256=digest(blas.resolve()),
        config=blas_config, header=str(include/"openblas_config.h"),
        header_sha256=digest(include/"openblas_config.h"), integer_abi="LP64/I4",
        license_path="/usr/share/licenses/openblas/LICENSE",
        license_sha256=digest(Path("/usr/share/licenses/openblas/LICENSE"))),
    abi_probe=dict(source_path=str(probe_source), source_sha256=digest(probe_source),
        command=probe_command, binary_sha256=digest(probe_binary), output=probe_output,
        ldd=subprocess.check_output(["ldd", str(probe_binary)], text=True)),
    compiler_version=subprocess.check_output([cache["CMAKE_C_COMPILER"], "--version"], text=True),
    prefix=str(root/"prefix"), cache={k: v for k, v in cache.items()
        if k.startswith(("BLA", "LAPACK", "SUITESPARSE_", "UMFPACK_", "CHOLMOD_",
                         "BUILD_", "CMAKE_INSTALL_", "CMAKE_BUILD_TYPE", "CMAKE_C_COMPILER", "CMAKE_C_FLAGS"))},
    libraries=[dict(path=str(p), resolved=str(p.resolve()), sha256=digest(p.resolve()),
        ldd=subprocess.check_output(["ldd", str(p)], text=True)) for p in paths],
    status="built; consumer correctness checks still required; no performance admission",
    baseline_scope="different dense runtime from original system/reference-BLAS checkpoint; old binaries/libraries unchanged")
(root/"provenance.json").write_text(json.dumps(record, indent=2)+"\n")
print("Provisioned", record["package_version"], "with", blas_config)
print("Prefix:", root/"prefix")
PY
