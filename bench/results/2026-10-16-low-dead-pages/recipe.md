# Reproduction scope

The archived `freeze.py`/`run.py` hashes identify the recipes actually used.
The parameterized copies here preserve the same two cases, alternating order,
watchdogs, observer and separate counter diagnostics; directory arguments make
archive relocation possible. They were syntax-checked and reviewed, without a
second measurement campaign. The original measurements remain identified by
`provenance.json`, not by a subsequent integration or rerun.

For numeric regeneration only, with no compilation or solver execution:

```sh
python regenerate.py "$FS_AB_ARCHIVE" --output "$FS_AB_TABLES"
```

This validates raw SHA256 and source identities, then regenerates all20 ordinary
rows, four min/median/max groups and both excluded diagnostic rows. Committed
`trials.csv`, `summary.csv` and `diagnostics.csv` contain the portable numeric
evidence; full stdout/observer receipts remain in the archive.

For a new measurement, first restore baseline8d7f634 and candidatefda5cbd in
separate source trees. The tested candidate patch was subsequently integrated
as4fd03ef; that integration is not a newly measured source identity. Configure
both trees identically, using different build directories:

```sh
cmake -S "$FS_AB_SOURCE" -B "$FS_AB_BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_Fortran_FLAGS_RELEASE=-O3 -march=znver3 -DNDEBUG' \
  '-DCMAKE_C_FLAGS_RELEASE=-O3 -march=znver3 -DNDEBUG' \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_METIS=OFF \
  -DFORTSPARSE_ENABLE_AMD=ON -DFORTSPARSE_ENABLE_NATIVE_CUDA=OFF \
  -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DFORTSPARSE_BUILD_TESTING=OFF \
  -DFORTSPARSE_BUILD_BENCHMARKS=ON
cmake --build "$FS_AB_BUILD" --target bench_native_families -j 2
```

Use the qualified C observer ELF and its adjacent `measure.c`. The committed
qualification receipt identifies the exact observed ELF and the heavy-parent
RSS/watchdog oracle; recompiling on another toolchain needs its own observer
qualification. `measure.c` is retained as the source recipe, not a portable
promise about another machine's inherited RSS behavior.

Freeze into a fresh archive, supplying both source/build trees:

```sh
python freeze_checkpoint.py --archive "$FS_AB_NEW_ARCHIVE" \
  --baseline-source "$FS_AB_BASE_SOURCE" --candidate-source "$FS_AB_PORT_SOURCE" \
  --baseline-build "$FS_AB_BASE_BUILD" --candidate-build "$FS_AB_PORT_BUILD" \
  --observer "$FS_AB_OBSERVER" --qualification "$FS_AB_QUALIFICATION"
```

This checks matched flags and resolved shared-library hashes, pins fixture/
driver sources, and links the distinct diagnostic ELF from the exact candidate
core archive/driver objects plus `advice_counter.c`. The wrapper proves helper
origin before counting every actual advised span and syscall failure.

Acquire the controller's quiet window before running:

```sh
python run_checkpoint.py --archive "$FS_AB_NEW_ARCHIVE" \
  --fixture-source "$FS_AB_PORT_SOURCE" --go
```

This runs exactly20 ordinary children on CPU30, one OpenMP/BLAS thread, five
alternating pairs for drift-real128 training RHS8 and refactor, each with a0.9s
wall watchdog. It then runs two instrumented children under a separate3s wall
watchdog. No diagnostic time/CPU/RSS enters the ordinary summary. There are no
retries or size/permutation/tolerance changes; failures stay in raw receipts.
The exclusive-create raw outputs deliberately reject accidental overwrite.
