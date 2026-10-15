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
