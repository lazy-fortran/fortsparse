# Native SPD performance evidence

Build the optional benchmark targets with the owning project's benchmark
option. Set `FORTSPARSE_BENCH_MKL_ROOT` to an installed oneMKL prefix to build
`bench_spd_pardiso`. It links only oneMKL and the fixture library; it never
loads the FortSparse/OpenBLAS backend stack. This measures **Intel oneMKL
PARDISO**, not Panua PARDISO. Python is optional orchestration tooling, not
a FortSparse runtime dependency. The runner requires Linux `taskset` and
sysfs CPU topology.

The frozen small corpus has three banded SPD sizes (128/512/2048 unknowns,
half-bandwidth four), anisotropic 2D grids (16/24/32 squared) and 3D grids
(5/7/9 cubed). Every matrix is a weighted graph Laplacian plus a 0.25 identity
shift, stored as full canonical CSC. Deterministic held-out edge coefficients
are selected with `--heldout`. Right-hand sides come from known alternating
solutions. Residual checks use the original full matrix independently of the
solver. The fixture test also checks symmetry and positive definiteness by
independent dense Cholesky.

```sh
python3 bench/run_spd.py --bin-dir build/bin --output /var/tmp/spd-small.json
```

Defaults are five fresh processes per sample group, one/eight RHS,
same-pattern changed-value refactorization, derivatives and warm-runtime
refactor-from-scratch, and 1/2/4/8/16 physical cores.
The CPU affinity is fixed, OpenBLAS uses one thread, OpenMP and oneMKL use
the requested thread count, and dynamic threading is disabled. SuperLU and
UMFPACK remain serial factorization baselines; sweeping their affinity does
not make them parallel implementations. Unbuilt backends, failed residuals
and timeouts remain explicit records. Narrow runs use `--backends`,
`--families`, `--levels`, `--threads`, `--workloads` and `--trials`.

Each process has a one-second timeout, including startup; a timeout kills
its entire process group, including the helper. The full benchmark collection
contains many short invocations and can take longer than a second. Calibrate
larger thread-scaling cases exclusively against external implementations:

```sh
python3 bench/run_spd.py --bin-dir build/bin --output /var/tmp/spd-calibration.json \
  --backends superlu,umfpack,mumps,mkl-pardiso --calibrate-ladder
python3 bench/run_spd.py --bin-dir build/bin --output /var/tmp/spd-scaling.json \
  --scaling-corpus /var/tmp/spd-calibration.json
```

The candidate ladder is fixed at 2D 48/64/96/128 and 3D 12/16/20/24. For each
family, calibration admits the largest candidate for which **every selected
external backend and thread count** completes in under 0.5 seconds. Preserve
that calibration record unchanged before tuning native kernels. Missing
external backends prevent admission; choose the available comparison set
explicitly and record the omission. Native failures never shrink the corpus.

Raw JSON retains commands, binary hashes/dynamic dependencies, Git revision,
CPU affinity, every trial, and median/minimum/maximum summaries. Times include
ordering, factorization and solves separately where the API exposes them,
plus total workflow and externally observed process wall time. `parent_cpu_s`
is user+system process CPU time, including all parent threads;
`helper_cpu_s` is reported separately after helper termination. Refactor time
reports whether the backend actually reused analysis. Missing backend
statistics remain `-1`; they are not inferred from total timings.

Memory reports Linux post-exec `VmHWM` for the parent and the reaped helper's
peak RSS separately. Do not sum the two: shared mappings may be counted in
both. These are process peaks, not attributed solver allocations. Backend
factor memory and fill are separate columns when available. CPU workflow
totals include fixture assembly, independent residual checking and cleanup;
`factor_solve_s` excludes fixture/residual/cleanup costs. The warm workload
frees the first factors then repeats analysis and factorization in a loaded
runtime; `cold_factor_s` preserves its first factorization. Its total workflow
and CPU totals include both factorizations, so use its factor-plus-solve
column when comparing amortized runtime costs. Changed-value refactorization
is a different workload and records actual symbolic reuse.

Derivative workloads time implicit JVP and VJP, including sparse matrix
cotangent construction, with the same retained SPD factors. They check a
known tangent and the adjoint identity independently. The isolated PARDISO
driver implements the equivalent implicit-equation wrapper around its solves.
No GPU performance claim is made by these drivers yet. Summaries include
speedup and efficiency against each backend's one-thread result; missing
statistics are JSON `null`. Compare both cold and warm cases, five-trial
spread, and best-tuned external thread counts before claiming competitiveness.
UMFPACK's separate parent/helper memory peaks do not establish a unique
aggregate footprint; a combined memory claim needs additional evidence.

## Match the MUMPS ordering provider

A prebuilt MUMPS archive must link the METIS provider whose headers it was
compiled against. Matching the METIS version number alone is insufficient.
The local MUMPS 5.9.1 package used a bundled METIS 5.2.1 header with
`METIS_OPTION_NUMBERING=18`; the system METIS 5.2.1 header/library use index 19.
Combining those providers made MUMPS pass a one-based graph to a zero-based
METIS call and crash during nested dissection. The captured graph was valid;
calling the same library on the same graph with the correct numbering option
succeeded. This is a dependency ABI mismatch, not a sparse-factorization result.

For a MUMPS installation that provides its own matched METIS, select both
explicitly when configuring the benchmark build:

```sh
cmake -S . -B build \
  -DMUMPS_DIR=/path/to/mumps-prefix/cmake \
  -DMETIS_INCLUDE_DIR=/path/to/mumps-prefix/include \
  -DMETIS_LIBRARY=/path/to/mumps-prefix/lib/libmetis.a
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  timeout 1s build/bin/bench_spd mumps grid3 6 1 rhs1 training
```

The static METIS provider also needs its matching GKlib dependency. Its package
target should supply that dependency; a dynamic system METIS must not acquire
an unrelated static GKlib on the executable's link line. Preserve failed
provider configurations as failed evidence and rerun calibration after changing
the linked provider.

## Admission and profiling records

Backend names `native-reference` and `native-cuda` select the simplicial
reference and native device request respectively. Requested threads and the
native factor's actual active team are separate columns. A device request
does not establish GPU execution: `actual_device=-1` remains explicit until
the backend supplies verified device statistics.

Every run writes `.trials.csv` and `.summaries.csv` alongside its JSON. Trials
include failures/timeouts; summaries include CPU time, wall time, memory,
scaling and missing statistics. JSON retains compiler versions, CMake cache
hash/options, source-content hash, executable hashes and dynamic dependencies.

Admission uses factor-plus-solve for ordinary/warm RHS workloads,
refactor-plus-solve for changed-value reuse, and factor-plus-primal-solve plus
JVP/VJP for derivatives. Derivative reuse alone is also reported. Each case
compares the native path with the *same* external Pareto-frontier configuration
on both time and memory, considering all measured external thread counts.
A clear observed-range pass requires the native maximum to be no greater
than the comparator minimum on both axes. Overlapping observed ranges are
indeterminate; no percentage tolerance or confidence interval is invented.

UMFPACK's unknown aggregate peak is retained as unknown. It cannot dominate a
known frontier candidate when it is clearly slower, or when its parent RSS
lower bound already exceeds the candidate's peak. Otherwise that candidate's
frontier membership is indeterminate. Missing comparators also prevent a
complete admission claim. The reference path is labeled separately from the
production path, and this performance analysis does not itself establish all
MVP requirements.

The default comparison includes `native` (public AUTO ordering) and `native-nd`
(the same production backend with public nested-dissection ordering). Raw CSV
retains `ordering_used`, actual thread/device counts and attributed storage.
`admission.best_native_by_case` selects an explicitly measured CPU ordering and
thread configuration only when one strict comparison passes on both elapsed
time and parent peak RSS against the **same** protected external frontier
configuration. Its companion `best_native_overall` passes only if every
evaluated case/workload/variant group passes. Failed alternatives remain in
the per-backend admission and raw trials; missing comparators and uncertain
aggregate UMFPACK memory retain indeterminate status. Reference, CUDA and
external-fallback aliases cannot supply this production CPU union.

This union requires case-specific ordering and thread choices. It is not a
universal claim about AUTO/default performance. Plots show `Native CPU (ND)`
as a separate marker in the CPU family color, including CPU/CUDA comparisons;
they do not silently combine ordering variants.

`report_spd.py` combines raw run JSONs into a reviewable admission JSON, a
per-case CSV and a configuration/scaling CSV. It reuses the strict admission
logic; selected native and comparator rows include elapsed phase, parent RSS,
whole-workflow CPU ranges, requested/actual teams and source/binary identities.
Supply files oldest first. A newer executable or changed declared environment
replaces the backend's entire thread series for that case/workload/variant.
The same executable/settings may supply disjoint thread counts, while a newer
overlapping configuration replaces its previous trials instead of pooling
them. Failed newer configurations stay authoritative. Input hashes, original
status counts and source/build provenance retain the superseded history.

```sh
python3 bench/report_spd.py /var/tmp/fortsparse-mvp/external-final.json \
  /var/tmp/fortsparse-mvp/bcsstk09-external.json \
  /var/tmp/fortsparse-mvp/native-final.json \
  /var/tmp/fortsparse-mvp/native-heldout.json \
  --output /var/tmp/fortsparse-mvp/final-report.json
```

Defaults declare all 11 frozen synthetic cases plus HB/bcsstk09, five workloads,
both native CPU ordering configurations, four external backends, five trials
and 1/2/4/8/16 threads. Held-out variants apply to synthetic cases only. State
a narrower intended corpus explicitly with `--cases`, `--workloads`,
`--variants`, `--backends`, `--threads` and `--trials`; missing cases, thread
counts and trial counts remain visible in coverage. Incomplete declared
evidence makes `scope_admission_status` indeterminate, even if the measured
subset passes. Failures are recorded independently of coverage. The report
does not turn empirical performance admission into a complete MVP claim.

`admission.front_membership` is a separate nondominance diagnostic. The strict
match/beat gate above is unchanged: a native time/memory tradeoff can be
nondominated without matching or beating any single external frontier point
on both axes. For each native configuration, this diagnostic tests **every**
measured external configuration. It reports a robust observed-range front
member only if native maximum time is below external minimum time, or native
maximum parent RSS is below external minimum parent RSS, for every comparison.
UMFPACK parent RSS is used only as a lower bound on its unknown aggregate peak.
Missing/failed external evidence blocks a positive membership result. Overlap
remains indeterminate; a known comparator that dominates on both separated
axes is reported explicitly. Report CSVs retain separate front-membership
status and native configuration; these are not relabeled strict-gate passes.

For CPU stack profiling, replay complete solver operations for about 0.2 s:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 perf record -g -- \
  build/bin/bench_spd native grid2 3 4 profile-factor training
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 perf record -g -- \
  build/bin/bench_spd native grid2 3 4 profile-refactor training
```

`profile-factor` frees/rebuilds analysis and factors then solves each iteration.
`profile-refactor` retains analysis and alternates matrix/RHS scales, then
solves. CSV records iteration count, replay wall time and average phase times;
the final solve is independently residual-checked. The replay stops before
its predicted next iteration exceeds a 0.9 s internal budget. Choose a smaller
profiling case if the warmup already consumes that budget. Standard runner
invocations still have a one-second outer limit. Profiling modes are excluded
from admission; instrumentation overhead is not an ordinary performance result.

## Performance PNGs

`plot_spd.py` uses Matplotlib with the Okabe–Ito palette and redundant markers.
It writes Pareto, scaling and CPU/CUDA PNGs plus a manifest and compact numeric
CSV. Keep PNGs outside the repository and publish them through slopbox.

```sh
python3 bench/plot_spd.py /var/tmp/fortsparse-mvp/mvp-comparison.json \
  /var/tmp/fortsparse-mvp/cuda-comparison.json \
  --outdir /var/tmp/fortsparse-mvp/plots/latest --workload warm \
  --cases grid2:7,grid3:7 --scaling-cases grid2:7,grid3:7 \
  --retain-csv bench/results/plot-inputs.csv
```

Pareto panels show elapsed solver phases and whole-workflow CPU time versus
parent peak RSS. Scaling curves use the same executable within each backend;
a newer binary replaces that backend's complete thread series. Legends state
the actual trial counts. Warmup CPU work is included in whole-workflow CPU
time, while elapsed warm factor/solve excludes warmup. UMFPACK's parent and
helper RSS stay separate, and CUDA device allocations are reported separately
when available. Failed and timed-out configurations remain in the manifest
and annotations. Preliminary snapshots from different revisions are labeled
as checkpoints and do not establish MVP admission.

The direct `bench_spd` executable accepts an optional seventh integer argument
for the public ordering option (1 = AMD, 2 = METIS ND). This is a diagnostic
common-ordering override; ordinary runner comparisons use each backend's
unchanged default policy. The CSV reports the actual native ordering, symbolic
fill, numeric factor payload and retained workspace separately from peak RSS.

The retained compact CSV is also a direct plotting input, so source history
can reproduce the measured curves without the larger raw JSON files:

```sh
python3 bench/plot_spd.py bench/results/plot-inputs-progress.csv \
  --outdir /var/tmp/fortsparse-mvp/plots/replay --workload rhs1 \
  --cases grid2:7,grid3:7 --gpu-cases grid2:1,grid2:3
```

## Frozen real-world structural matrix

The optional `--matrix` case is HB/bcsstk09: the 1083-by-1083 stiffness matrix
of a clamped square plate, supplied by John Lewis (Boeing) in 1982. Both the
[SuiteSparse collection](https://sparse.tamu.edu/HB/bcsstk09) and
[NIST Matrix Market](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk09.html)
identify it as real symmetric positive definite. Selection used only that
primary metadata and bounded order, before native performance evaluation.

The exact archive is available from
[NIST](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc1/bcsstk09.mtx.gz).
The local research archive is `~/code/sparse-dev/matrices/HB/bcsstk09/`, with
its library card and catalogue under `~/code/sparse-dev/`. Raw matrix files
stay out of source history. `fixtures/bcsstk09.json` freezes the source URI,
archive/matrix SHA256 and dimensions. The runner rejects different bytes.

```sh
python3 bench/run_spd.py --bin-dir build/bin --output /var/tmp/bcsstk09-external.json \
  --matrix "$HOME/code/sparse-dev/matrices/HB/bcsstk09/bcsstk09.mtx" \
  --backends superlu,umfpack,mumps,mkl-pardiso
```

Calibrate the external case first, retaining failures and the one-second
invocation limit, then freeze that evidence before native tuning. This case
is separate from the synthetic ladder. The narrow Fortran reader accepts
lower-triangle real-symmetric coordinate MatrixMarket, mirrors off-diagonal
entries, sorts columns and sums duplicate records. Both benchmark executables
use it. Known-solution RHS construction, residual checks and the existing
derivative oracle/tolerances are unchanged; collection SPD metadata does not
replace those checks.

## Low-memory CPU comparison

`native-low-memory` selects the native backend with
`FORTSPARSE_MEMORY_LOW`. The current supported schedule is natural-order CPU
real-SPD simplicial Cholesky; requested factor threads can therefore exceed the
reported active team. This is a memory-policy choice, not an external fallback.
The optional ninth `bench_spd` argument sets `host_memory_limit_bytes` (zero
means unlimited); argument seven remains ordering and eight remains matrix path.
Explicit AMD/ND with this policy is currently unsupported. See
[the memory contract](../docs/memory.md) for the exact cap scope.

The existing CSV columns retain their names/order, with numeric memory columns
appended. `plan_s` measures the separate public planning call;
`memory_plan_status` must be zero before using `memory_prediction_quality`
(1 exact layout, 2 bound, 3 estimate) or `predicted_host_peak_bytes`.
Unavailable predictions are `-1`. `host_current_bytes` and `host_peak_bytes`
are attributed solver payload measured after the workload, before freeing
factors. `host_memory_limit_bytes`, `memory_policy`, `strict_memory_supported`
(0/1), and required/available bytes on allocation rejection remain separate.
None predicts or limits process RSS. Planning excludes caller A/B/X, runtime,
allocator and object overhead; its time is included in total wall/CPU time,
not factor/solve time. Uncapped production and external paths can report an
unsupported planner and unknown attributed payload without failing the solve.
A failed factor prints `memory_failure` diagnostic CSV and exits unsuccessfully;
the runner retains its failure output rather than treating it as a timing sample.

After builds and unrelated measurements have stopped, use a declared narrow
case and five fresh processes per backend. Existing input/oracle/tolerance,
one-second process limit, affinity, source/binary/library identities and
parent/helper CPU/RSS accounting remain unchanged:

```sh
python3 bench/run_spd.py --bin-dir build/bin \
  --output /var/tmp/spd-memory-small.json --threads 1 --trials 5 \
  --families band,grid2 --levels 1 --workloads rhs8 \
  --backends native,native-low-memory,native-reference,superlu,umfpack,mumps,mkl-pardiso
```

The raw JSON records retain all appended fields; existing generic performance
summaries do not yet aggregate these memory columns. Read predicted/current/peak
payload alongside parent peak RSS, not as its replacement. Native reference is
an independent schedule label, and unavailable external executables remain
missing-comparator records. Report medians and observed ranges for complete
wall/CPU/RSS and factor/solve costs, with background load and reuse stated.
Shared load, failures or missing comparators prevent a performance-admission
claim. These diagnostics cannot establish a performance win from lower payload.

The ninth-argument cap can also be checked directly, keeping positional slots:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 1s \
  build/bin/bench_spd native-low-memory band 1 1 rhs8 training 0 '' 1048576
ctest --test-dir build -j 2 -R '^spd_benchmark_memory$' --output-on-failure
```

The registered check exercises the existing CSV reader, independently checked
eight-RHS solves, an exact predicted cap and rejection one byte below it. It is
correctness evidence, not ordinary performance measurement. Keep raw comparison
JSONs and generated plots outside source history; publish plots only through
the controller's existing slopbox workflow.
