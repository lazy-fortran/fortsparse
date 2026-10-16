# Single-node native-family comparison foundation

These benchmark-only fixtures cover real LU, complex LU, complex Hermitian
positive-definite Cholesky and real symmetric-indefinite LDLT. They complement
the existing SPD corpus; no existing SPD fixture, driver or scientific source
is changed. Performance admission is **unknown**. Small correctness examples
do not establish throughput, scaling, representativeness of a production NEO-2
solve, or membership in a time/memory frontier.

The frozen definitions in `fixtures.f90` are checked by `contract.json` SHA256.
Three mathematical families provide five cases:

- `neo2-mini5`: actual retained NEO-2 integer matrix and rational oracle from
  `../fixtures/neo2-mini5/`; training normal RHS1 is the original all-ones RHS.
  Held-out RHS uses different independently known values; the matrix stays real.
- `drift-real` and `drift-complex`: strictly diagonally dominant 2D nonsymmetric
  stencils, with unequal left/right and up/down coefficients. They explicitly
  use GENERAL symmetry hints. Held-out coefficients change drift and phase.
- `hpd-chain` and `ldlt-chain`: conjugate-paired positive-dominant Hermitian
  chain, or coupled indefinite zero-diagonal 2x2 blocks requiring paired pivots.
  Held-out coefficients change the coupling. These are deliberately bounded
  algorithmic fixtures; larger/real consumer classes still need admission.

All RHS use independently known x. Checks evaluate the original full CSC
matrix, including genuine transpose or adjoint actions. The fixture test also
checks an independent dense action and the HPD dominance/indefinite structure.
A refactor scales A and b by 1.125 and reuses analysis where the backend does.

Build with benchmarks enabled and, when available, sequential OpenMP MUMPS:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFORTSPARSE_BUILD_BENCHMARKS=ON -DFORTSPARSE_USE_SYSTEM_DEPS=ON \
  -DFORTSPARSE_ENABLE_MUMPS=ON \
  -DMUMPS_DIR=/mnt/storage/fs1-envs/mumps-local/cmake \
  -DFORTSPARSE_BENCH_MKL_ROOT=/opt/intel/oneapi/mkl/2026.0
cmake --build build --target bench_native_families \
  bench_native_families_pardiso test_native_family_fixtures -j 2
ctest --test-dir build -R '^native_family_' --output-on-failure -j 2
```

The general driver takes backend, case, level, threads, workload, operation,
variant and mini5 matrix path. The Intel-only driver omits the backend argument:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/bin/bench_native_families \
  native drift-complex 3 1 rhs8 H heldout bench/fixtures/neo2-mini5/matrix.mtx
```

Workloads are `rhs1`, `rhs8`, `refactor`, `warm`; operations are `N`, `T`, `H`.
`warm` primes the numerical runtime with one factorization then repeats it;
`prime_factor_s` is separate and `total_s` still includes the complete process
workflow, including priming, fixture construction, checking and destruction.
`factor_s` includes public factor setup/analysis; `analysis_s`/`numeric_s`
separate reported phases. Refactor/solve are timed calls; total includes their
surrounding copies. The runner adds subprocess wall time including startup.

The CSV keeps parent CPU and terminated-helper CPU/RSS separate, plus actual
reported active factor teams, requested threads, factor family, reuse, fill,
owned payload when available, backward error and known-solution error. RSS
comes from existing `process_metrics.c` (`VmHWM` on Linux). Parent/helper RSS
must not be summed. Negative diagnostics mean unavailable. PARDISO memory is
vendor-reported, not exact owned payload; its actual active team is unknown.
A serial native factor family is flagged even when more threads were requested.
External FortSparse bindings at base revision a3bd4a2 lack retained T/H solve methods and
return backend unavailable there; those missing comparisons remain explicit.
The underlying external libraries may offer capabilities absent from the binding.

`bench_native_families_pardiso` links fixtures plus Intel oneMKL only; it never
links FortSparse or OpenBLAS. Its LP64 C ABI was checked against installed
`mkl_pardiso.h` and `mkl_service.h`; their hashes are in the contract. It uses
mtypes 11, 13, 4 and -2, verified in the retained official Intel Fortran manual,
not Panua examples. True CSR is constructed; HPD/LDLT use the upper triangle.
Intel `iparm(12)` maps N/H/T to 0/1/2. Hermitian transpose additionally uses
conjugation around a normal solve. C service symbols use the SDK's uppercase
`MKL_Set_Num_Threads_Local` and `MKL_Set_Dynamic`; lowercase aliases have the
Fortran by-reference ABI. This driver is Intel PARDISO, not licensed Panua
PARDISO. Missing SDK, library, runtime or license access remains a recorded gap.

`run.py` defaults to five fresh processes per point and physical-core sweeps
1/2/4/8/16 with fixed affinity, explicit thread settings and a 0.9 s process
limit. Timeout kills the entire process group, including IPC helpers. It
records commands, environment, binary hashes/ldd, raw stdout/stderr, failures,
subprocess wall time and median/ranges; it never treats RSS sums as memory.
Runs need a coordinated quiet measurement window. First calibrate external
backends only, without timing native execution to choose sizes:

```sh
python bench/native_families/run.py --build build --mode calibrate \
  --output /home/ert/code/sparse-dev/benchmarks/DATE-native-families-calibration
```

The generated `calibration_selection.json` admits a candidate size only when
all requested external configurations/trials succeeded below 0.5 s for cold
factorization plus normal RHS8. This initial calibration uses one provisional
level per case (500 attempts with four backends, five thread budgets and five
trials), not the full workload cross product. Peers are interleaved at each
matched point with rotating backend order, fixed affinity and recorded host
load. Missing
comparators, failures and timeouts remain raw receipts and prevent silently
claiming a complete selection. Freeze that file before tuning, then pass it
with `--levels-from PATH` in measurement mode. No performance levels or
calibration success are asserted in the initial source contract. Measurements
can be restricted with `--cases`, `--backends`, `--threads`; exclusions must be
reported. `--workloads`, `--operations` and `--variants` select bounded follow-up
subsets; record those exclusions. `--window-note` retains the controller
measurement-window qualification. `--dry-run --mode calibrate` prints commands
without execution.
