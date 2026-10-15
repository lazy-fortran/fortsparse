# SPD numerical source

`gen_spd_leaves.f90` retains the scalar square-root, division and Schur-update
definitions. FortSym differentiates these expressions and emits CPU and native
CUDA values/JVP/VJP from the same mathematics. The Fortran elemental function
wrappers contain no independent numerical formula. The domains are positive
inputs for the square root and nonzero denominators for division; matrix and
pivot validation belong to the calling solver.

`spd_schedule.f90.inc` retains allocation-free CPU panel factorization,
right-side triangular solve and rectangular Schur-update schedules. The solver
owns task parallelism; the panel schedules use SIMD without nested teams.
Large contractions may use BLAS rather than these reference schedules.

From the repository root, with a built FortSym at the retained revision:

```sh
bash gen/regenerate_spd.sh /path/to/fortsym
git diff --exit-code -- src/generated/
```

The generator uses FortSym's native simplifier. Ordinary library builds use
the committed output and require neither FortSym nor a CAS. Compile the
generated Fortran module into `fortsparse`; register
`test/solver/test_fortsparse_spd_leaves.f90` with the existing focused test
helper. The generated CUDA header contains device leaves only: launch geometry,
residency and stream ownership belong to the optional native CUDA adapter.

The focused test uses known values, independent central differences and
adjoint pairings, a constructed triangular-factor oracle, and dense `matmul`
residual/update oracles. Numerical probes establish behavior, not a proof or
an outward-rounded enclosure.
