# Block-tridiagonal symmetric indefinite factorizations

A symmetric block-tridiagonal matrix with dense blocks, the shape of a radial
finite-element stability pencil, is factored by the supernodal LDL^T before
and after the setup changes of this record: natural order by default for
patterns that fill their envelope, a linear symmetry check, parallel panel
assembly without a second validation, and a panel chain factored outside the
task team. A first factorization now costs a quarter of what it did, and a
refactorization, which keeps the analysis, is on par with a LAPACK block
elimination at block width 104 and faster at width 312. Inertia is exact and
reported as reliable in every run.

## Case

`bench_block_tridiagonal BLOCKS WIDTH DECADES MODE REPEATS` builds
A = S L D L^T S with a unit lower block-bidiagonal L of random entries, a
diagonal D of prescribed signs and a diagonal row scaling S that spans DECADES
orders of magnitude, so Sylvester's law gives the inertia exactly. One process
times one MODE: its first factorization, then the fastest of REPEATS
refactorizations and solves. MODE `lapack` eliminates the same blocks with
DSYTRF, DSYTRS and DGEMM, as a block-tridiagonal solver would, and has no
analysis, so its two columns are the first and the fastest repeat.

`run_block_tridiagonal.py` starts five fresh processes per configuration,
pinned to the first THREADS allowed cores, with DECADES 12 and a one-second
limit per process. The shapes are 64 x 104 (three repeats), 256 x 104 and
64 x 312 (one repeat each).

## Result

Wall time in milliseconds, median [range] of five processes. "before" is the
same tree with `src/` of c128826 (main before the change).

| blocks x width | threads | first factor before | first factor after | refactor before | refactor after | LAPACK blocks |
|---|---:|---:|---:|---:|---:|---:|
| 64 x 104 | 1 | 116 [116-117] | 30 [30-45] | 18 [18-19] | 17 [17-17] | 16 [16-16] |
| 64 x 104 | 8 | 117 [116-118] | 26 [25-26] | 18 [18-19] | 13 [13-13] | 15 [15-17] |
| 256 x 104 | 1 | 461 [459-461] | 117 [116-117] | 74 [73-74] | 69 [69-70] | 64 [64-64] |
| 256 x 104 | 8 | 461 [460-463] | 99 [99-100] | 74 [74-74] | 52 [52-52] | 62 [62-70] |
| 64 x 312 | 1 | 1255 [1250-1259] (a) | 295 [293-295] | 185 [184-188] (a) | 175 [174-175] | 309 [307-312] |
| 64 x 312 | 8 | 1256 [1252-1264] (a) | 223 [222-233] | 185 [185-186] (a) | 106 [106-117] | 129 [128-130] |

(a) Over the one-second limit: all five processes were killed in the regular
run (`before/`), and the numbers come from a rerun without the limit
(`before-over-limit/`).

A first factorization still costs 1.8 times the block elimination at width 104
on one thread, because the symmetry check, the symbolic analysis (9 ms of the
30 ms at 64 x 104) and the panel assembly each pass over all entries. The
numeric phase alone takes 14 ms there, against 16 ms for the block elimination.

Refactorization at 1, 2, 4, 8 and 16 threads after the change, wall time and
process CPU time in milliseconds:

| blocks x width | | 1 | 2 | 4 | 8 | 16 |
|---|---|---:|---:|---:|---:|---:|
| 64 x 104 | native wall | 17 | 15 | 13 | 13 | 14 |
| | native CPU | 17 | 30 | 55 | 104 | 215 |
| | LAPACK wall | 16 | 15 | 15 | 15 | 23 |
| 256 x 104 | native wall | 69 | 61 | 54 | 52 | 55 |
| | native CPU | 69 | 113 | 192 | 358 | 756 |
| | LAPACK wall | 64 | 63 | 61 | 62 | 94 |
| 64 x 312 | native wall | 175 | 137 | 114 | 106 | 120 |
| | native CPU | 174 | 244 | 371 | 646 | 1496 |
| | LAPACK wall | 309 | 203 | 147 | 129 | 275 |

The speedup from 1 to 8 threads is 1.3 at width 104 and 1.6 at width 312, so
the factorization of a single panel chain stays limited by its serial parts,
and CPU time grows almost in proportion to the thread count. Before the
change the thread count had no effect on the wall time.

The process peak of resident memory at 64 x 104 on one thread is 65 MiB for
the native factorization, which holds the assembled matrix and the factor,
and 23 MiB for the LAPACK elimination, which holds only the blocks; at
64 x 312 it is 533 MiB against 193 MiB.

Other solvers on one thread, first factorization / refactorization in
milliseconds and the componentwise backward error of one solve:

| blocks x width | native supernodal | native scalar with BLAS-3 panels (`blocked`) | SuperLU 7.0.0 | UMFPACK 6.3.5 (SuiteSparse 7.8.3) |
|---|---:|---:|---:|---:|
| 64 x 104 | 30 / 17, 1.2e-15 | 115 / 111, 1.7e-15 | 73 / 60, 1.7e-7 | 95 / 78, 1.4e-15 |
| 256 x 104 | 117 / 69, 1.4e-15 | 464 / 454 (b), 2.0e-15 | 280 / 232, 7.3e-7 | 362 / 315, 1.8e-15 |
| 64 x 312 | 295 / 175, 2.1e-15 | 19796 / 19802 (b), 3.0e-15 | 1118 / 1040 (b), 4.1e-6 | 857 / 769 (b), 4.5e-15 |

(b) Killed at the one-second limit in `after/`; numbers from
`after-over-limit/`.

SuperLU and UMFPACK factor the matrix as a general one and report no inertia.
SuperLU with its default settings does not reach the accuracy of the others
on this row scaling, so its times are not at matched accuracy. MUMPS, the
external solver with a symmetric indefinite factorization, is missing: the
FortSparse MUMPS backend needs an initialized MPI, which this driver does not
provide, and the machine has no MUMPS built with OpenMP for this compiler. Its
row is kept as FAILED.

## Files

`before/`, `after/` and the two `-over-limit` directories each hold `raw.csv`
(every process, with CPU time of the timed sections, CPU time of the whole
process, peak resident memory, inertia and backward error), `summary.csv`
(median, minimum and maximum per configuration) and `manifest.txt` (source
revision, compiler, flags, command line, binary hash, processor and linked
libraries). The peak memory of the UMFPACK rows is that of the calling
process, not of the helper process that holds the factor.

## Reproduction

One node of CINECA Pitagora (two AMD EPYC 9745, 16 cores allocated to the
job), GNU Fortran 12.3.0, OpenBLAS 0.3.26 built with OpenMP, CMake Release
build with SuperLU and UMFPACK built from source by the FortSparse build.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DFORTSPARSE_BUILD_BENCHMARKS=ON
cmake --build build --target bench_block_tridiagonal fortsparse_umfpack_helper
python3 bench/run_block_tridiagonal.py \
    --binary build/bin/bench_block_tridiagonal --output after
# rerun configurations that hit the limit, as MODE:THREADS:BLOCKS:WIDTH
python3 bench/run_block_tridiagonal.py \
    --binary build/bin/bench_block_tridiagonal --output after-over-limit \
    --timeout 120 --select blocked:1:64:312,superlu:1:64:312
```

For "before", replace `src/` first: `git checkout c128826 -- src`.
