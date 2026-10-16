# Single-node native-family comparison foundation

Hashed private real consumer matrices use the bounded `real-mm`/`rhs4`
[adapter and independent audit](PRIVATE_REAL_FIXTURES.md). Numeric inputs and
candidate solutions stay outside Git; the original frozen corpus below is unchanged.

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
  bench_native_families_pardiso bench_native_families_umfpack \
  test_native_family_fixtures -j 2
ctest --test-dir build -R '^native_family_' --output-on-failure -j 2
```

The general driver takes backend, case, level, threads, workload, operation,
variant and mini5 matrix path. The Intel-only and direct UMFPACK drivers omit
the backend argument:

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

The CSV now records cumulative parent user+system CPU through the final
`RUSAGE_SELF` sample as `parent_cpu_s`, including loader/startup CPU. The old
`CPU_TIME` interval from Fortran main entry is retained as `parent_main_cpu_s`.
Receipts produced before this instrumentation change have only that interval
under the old `parent_cpu_s` name; never relabel them as full-process CPU.
Final CSV printing/termination is after the cumulative sample. Helper CPU is
separate full-lifetime `RUSAGE_CHILDREN` CPU after reaping. `workflow_cpu_s`
prefers the runner’s completed-process `wait4` user+system CPU, including
reaped descendants; helpers must not be added again. The guarded pre-CSV sum
is only a fallback when `wait4` is unavailable. `process_tree_cpu_s` is -1
on timeout: killed or unreaped descendants make whole-tree CPU incomplete.
`wait4_reaped_cpu_s` retains the partial launch usage separately. This CPU
counter is not an aggregate RSS measurement. Linux propagates waited-child
usage through ancestors (getrusage(2), RUSAGE_CHILDREN; wait4(2)).

The CSV keeps parent CPU and terminated-helper CPU/RSS separate, plus actual
reported active factor teams, requested threads, factor family, reuse, fill,
owned payload when available, backward error and known-solution error. RSS
comes from existing `process_metrics.c` (`VmHWM` on Linux). Parent/helper RSS
must not be summed. Helper RSS is a `RUSAGE_CHILDREN` launch-lifecycle peak
and can retain pre-exec/inherited high-water; it is not an exact post-exec
helper peak or a safe aggregate lower bound. UMFPACK aggregate memory remains
unknown; the measured parent `VmHWM` alone is a safe lower bound. Negative
diagnostics mean unavailable. PARDISO memory is
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

The frozen scan uses drift levels 16/32/64/128 (n=256/1024/4096/16384) and
chain levels 128/512/2048/8192 (n=256/1024/4096/16384). Drift levels are bounded
by 128; chain levels separately by 8192. Mini5 remains latency/correctness only.
Fixture assembly, CSR conversion and original-matrix checking require O(n+nnz)
storage and work (bounded stencil degree), never a dense large-case oracle.
The fixture test checks the maximum dimensions with an independent tridiagonal
action. Dense fixture checks remain confined to the small smoke examples.

The generated `calibration_selection.json` admits a candidate only after every
declared external/backend/thread/trial point succeeds below 0.5 s process wall
for cold factorization plus normal RHS8. Peers are interleaved with rotating
backend order, fixed affinity and recorded host load. Scan levels ascend;
the first failed, unavailable, timeout or slow attempt ends that family's scan.
Every attempted receipt is retained, the rejected candidate is classified, and
all larger candidates remain explicitly unmeasured. Thus selection means the
largest passing **scanned prefix**, not a claim of monotonic solver performance
or an extrapolation beyond the first failure. Missing comparators are recorded
separately from size/time or numerical failures.

With six external configurations, five thread budgets and five trials, a completely
passing four-level scan needs at most 2550 processes across these five cases.
Each failure ends its family immediately; the process watchdog is 0.9 s.
Native execution never participates in calibration or size selection. No
ordinary measurements accompany this source change. Freeze the raw receipts
and resulting selection before tuning or timing larger native cases.

`--levels-from PATH` verifies the raw receipt hash, selected family/candidate,
complete passing coverage, and the declared SuperLU/UMFPACK-IPC/UMFPACK-direct/MUMPS/Intel-PARDISO
baselines at all requested 1/2/4/8/16 thread budgets and at least five trials.
A selection from fewer declared baselines is rejected by default. Explicit
`--diagnostic` permits an incomplete selection, labels its provenance as a
diagnostic subset, and cannot establish performance admission. It still needs
complete passing raw observations for its declared subset. No performance
levels or calibration success are asserted by the source contract. Measurements
can be restricted with `--cases`, `--backends`, `--threads`; exclusions must be
reported. `--workloads`, `--operations` and `--variants` select bounded follow-up
subsets; record those exclusions. `--window-note` retains the controller
measurement-window qualification. `--dry-run --mode calibrate` prints commands
without execution.

`native-amd` explicitly selects native LU with AMD ordering while retaining its
own CSV backend label. It does not change the default native ordering. Non-LU
uses may honestly report unsupported. All native aliases are removed before
external-only calibration and rejected as size-selection comparators.


The initial-checkpoint report recipe preserves 500 matched groups, including
all timeout counts and successful-trial minimum/median/maximum values:

```sh
python3 bench/native_families/report.py /path/to/initial/measurement \
  --output /external/initial-analysis --plots \
  --cpu-scope 'parent main-entry CPU_TIME interval; helper lifetime CPU separately'
```

The default `--cpu-model initial` preserves the initial checkpoint's CPU
interpretation. Use `--cpu-model completed` only with corrected wait4 receipts;
the report rejects successful receipts lacking that counter. Both checkpoints'
reviewable numeric evidence is retained in
`bench/results/native-families-2026-10-16`. PNGs stay outside Git.
The recipe records the raw attempts hash and filtering. Error bars span every
successful trial; timeout crosses denote the configured censoring deadline.
Per-panel axes adapt to the fixture's regime. The compact CSV retains solve,
refactor, main-workflow and subprocess timings separately, actual factor-team
reports (unknown externally) and memory diagnostics without aggregating RSS.

The current comparator set includes `umfpack` (the FortSparse IPC backend),
`umfpack-direct` (standalone GPL general LU with native defaults) and
`umfpack-direct-noir` (the same isolated executable with only IRSTEP=0).
The direct target is optional when UMFPACK is found and never links FortSparse.
It reuses these same fixtures, records, process metrics and original-CSC oracles.
Its packed complex I4 API maps ordinary T to `UMFPACK_Aat` and adjoint H to
`UMFPACK_At`. Symbolic analysis is retained across numerical refactorization;
the selected controls persist across refactorization and warm-context recreation. HPD/LDLT cases are compared using
UMFPACK's general LU, not mislabelled Cholesky/LDL execution.

The standalone process has no helper or IPC: its Linux `VmHWM` is the complete
postexec process peak, and is an available memory comparator. Vendor factor
memory separately multiplies `UMFPACK_PEAK_MEMORY` Units by
`UMFPACK_SIZE_OF_UNIT` bytes/Unit, then converts to decimal MB. It is not an
exact owned-allocation count; owned fields stay unavailable. Default thread
budgets/environment and actual loaded libraries remain in runner provenance;
actual active UMFPACK teams are not reported. The installed library inspected
for this addition links reference BLAS/LAPACK, so these source checks establish
no best-tuned UMFPACK performance claim. See `umfpack_direct/README.md`.

Historical contract schema 3 added the fifth required external baseline. Existing receipts
and four-comparator selections keep their original historical scope; they do
not silently acquire the new baseline. Those old selections retain their original scope.
No ordinary measurements or calibration were performed for this source change.

The scalar/indexed diagnostic has a separate recipe that keeps the revision
label in every group and pairs matching trial IDs before forming speedups:

```sh
python3 bench/native_families/report_ab.py \
  /home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/indexed-ab \
  --output /external/indexed-ab-analysis --plots
```

This recipe verifies frozen binary hashes, records complete wait4 CPU and retains
all medians/ranges without treating this native-only comparison as admission.

`report_front.py` combines the completed five-peer drift/chain receipts, validates
all five external peer/thread configurations and separates robust observed-range
nondominance, known domination and indeterminate comparisons. It reads actual
matrix dimensions from the numeric records. Its retained checkpoint and recipe
are in `bench/results/native-families-2026-10-16/five-peer-openblas`.

`report_private_real.py` reads the immutable completed private-real measurement
receipts without executing a solver or rereading numeric candidates. It
normalizes wrapper, child and independent-audit statuses conservatively before
using the shared summaries/range assessments. Its explicit scalar whitelist
publishes only aggregate benchmark metrics, status counts and hashes; numeric
inputs, candidates, exact per-RHS fractions and the full private fixture
manifest remain private. The qualified fixture/policy guard prevents reporting
a different private case with this recipe. See the retained checkpoint in
`bench/results/private-neo2-2026-10-16`.

## Bounded AMD memory selections

The benchmark driver accepts `native-amd-low` (AMD with `MEMORY_LOW`, optional
host cap, default unlimited) and `native-amd-capped` (AMD with `MEMORY_AUTO` and an explicit positive
host cap). These are benchmark labels for existing solver options. They add no
production solver option and do not select synthetic sizes.

For `run.py` or `run_private_real.py`, capped execution requires
`--native-host-limit-bytes BYTES`, an exact positive decimal integer through
9223372036854775807. The runner passes `FORTSPARSE_BENCH_HOST_LIMIT_BYTES` to
the capped child and to LOW when the positive cap is explicitly requested. Direct capped-driver execution requires that environment
variable; missing, zero, negative, malformed or overflowing input fails before
factorization. LOW accepts the same optional positive cap; without the option it remains
unlimited. The runner clears an inherited cap variable so omitted options do
not silently change that default. Direct LOW execution can use the environment
variable; invalid supplied values fail before factorization.

CSV records include actual `memory_policy` (0 AUTO, 1 LOW) and exact integer
`host_memory_limit_bytes` (0 unlimited; -1 unavailable for comparators that do
not report these options). Runners preserve the integer in raw
records, summary configuration and requested provenance, confirm that successful
memory-alias records match the requested configuration, and retain mismatches as
failures. The private runner uses the same controls after its unchanged external
qualification gates; numeric fixtures and candidates stay in the private archive.

All four native labels (`native`, `native-amd`, `native-amd-low`,
`native-amd-capped`) are excluded from external calibration. Existing frozen
fixture contracts, source hashes and comparator IDs are unchanged. Plotting
includes both memory aliases with distinct styles; unknown styles and mixed
policy/cap records within one group are rejected. Caps remain explicit in numeric
CSV/provenance; do not merge separate cap campaigns under one backend label.
Owned allocation payload and process RSS remain separate, and unavailable runs
are retained. These driver checks do not establish a LOW/capped performance win.


Contract schema 4 requires six external configurations, retaining direct UMFPACK
both with its default `UMFPACK_IRSTEP=2` and with only that control set to zero.
The IPC helper uses IRSTEP=0 when `solver%refine` is false, as in these native
family runs. This adds a matched no-refinement alternative without discarding
the default comparison. No ordering, scaling, pivot, scientific fixture or
accuracy gate changes. Both variants must independently pass the unchanged
original-matrix backward/known-solution gates.

The direct driver's optional eighth argument is `default` (also the omitted
argument) or `no-refinement`. The runner supplies it only for
`umfpack-direct-noir`; CSV backend labels stay distinct. A bounded
`umfpack_controls,` JSON line records the runtime version, zero-based IRSTEP
index, complete native-default vector and complete effective vector after
execution. Receipts retain this snapshot and reject a successful result whose
label/controls differ from the request, including any override outside IRSTEP.
The arrays passed to vendor calls are exactly those recorded; the official
UMFPACK API declares Control const. Warm and refactor paths retain this choice.

Existing four/five-comparator selections and reports remain immutable historical
evidence. New matched admission needs a complete schema-4 six-configuration
calibration; older selections require explicit `--diagnostic` use.
`report_front.py --coverage historical-five` remains the default historical
verifier. `--coverage contract-v4-six` requires six-comparator calibration,
schema-4 provenance and confirmed controls, and compares all 30 external
configuration/thread combinations. The private runner likewise requires all
six qualified configurations; its old fixed-scope report remains historical.

Best-tuned external baseline coverage remains open until calibration and matched
measurements are complete. Testing this justified IRSTEP variation does not
establish exhaustive tuning of ordering, pivoting, scaling, dense runtime or
other documented vendor choices. No new ordinary timings accompany the source
change; existing default-mode results retain their recorded controls and scope.
