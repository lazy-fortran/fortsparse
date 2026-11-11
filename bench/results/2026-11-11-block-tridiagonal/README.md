# Block-tridiagonal symmetric indefinite factorizations

`bench_block_tridiagonal BLOCKS WIDTH DECADES MODE REPEATS` factors
A = S L D L^T S with a unit lower block-bidiagonal L, a diagonal D of prescribed
signs and a diagonal row scaling S that spans DECADES orders of magnitude, so
Sylvester's law gives the inertia exactly. The pattern is that of a radial
finite-element stability pencil: dense blocks on three block diagonals. The
driver reports the best of REPEATS for a first factorization, a same-pattern
refactorization and one solve, and for an elimination of the same blocks with
LAPACK DSYTRF, DSYTRS and DGEMM, which a block-tridiagonal solver would use.

`before.csv` is the state before the setup changes (AMD by default, one
symmetry search per entry, serial panel assembly, task team for the panel
chain), `after.csv` the state with them. Both use the row-local inertia
thresholds, and every run returned the exact inertia as RELIABLE with a
componentwise backward error below 3e-15.

Wall time in milliseconds, MODE `supernodal`, DECADES 12:

| blocks x width | threads | first factor before | after | refactor before | after | LAPACK blocks |
|---|---:|---:|---:|---:|---:|---:|
| 64 x 104 | 1 | 123 | 27 | 19 | 16 | 16 |
| 256 x 104 | 1 | 494 | 106 | 77 | 64 | 64 |
| 64 x 312 | 1 | 1341 | 278 | 197 | 169 | 301 |
| 64 x 104 | 8 | 121 | 22 | 20 | 12 | 16 |
| 256 x 104 | 8 | 491 | 88 | 76 | 49 | 69 |
| 64 x 312 | 8 | 1271 | 208 | 180 | 101 | 129 |

The first factorization still costs 1.7 times the LAPACK block elimination at
width 104, because the symmetry check, the symbolic analysis and the panel
assembly each pass over all entries; the numeric phase alone takes 13 ms. A
refactorization, which keeps the analysis, matches the block elimination at
width 104 and is 1.8 times faster at width 312 on one thread.

MODE `blocked` (scalar Bunch--Kaufman with BLAS-3 panels) changes with the
default order only: 145 to 107 ms at 64 x 104 and 21.7 to 19.5 s at 64 x 312
on one thread.

Measured on one node of CINECA Pitagora (two AMD EPYC 9745, 16 cores allocated),
GNU Fortran 12.3.0, OpenBLAS 0.3.26 built with OpenMP, `OMP_PLACES=cores`,
`OMP_PROC_BIND=close`.
