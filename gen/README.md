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

The optional adapter in `src/native/fsparse_spd_cuda.cu` owns one retained
device arena, nonblocking stream, cuBLAS handle and explicit 1 MiB cuBLAS
workspace. Its thin C interface uses arena offsets, so factor panels and
update operands can remain resident between calls. Allocation occurs at
context creation. The reported allocation bytes cover the arena, pivot info
and explicit workspace; CUDA context/cuBLAS internal memory must additionally
be measured by the benchmark process/device accounting.

Small diagonal panels use generated numerical leaves in one cooperative
CUDA block. Larger triangular solves and Schur contractions lower to cuBLAS.
Upload/download calls enqueue copies; callers must synchronize before reusing
host buffers and check factor info before consuming a failed panel. Each
context belongs to one caller and one stream. The supernode engine's optional
device argument selects a resident sparse factorization pipeline: setup
retains compact integer scatter maps, each factorization uploads its panel
input once and downloads its completed factors once, and CPU solves use the
same host allocation. Panels and Schur updates stay resident between those
transfers. This establishes an experimental execution path; competitiveness
still requires transfer-inclusive measurements on the declared workload.

For optional integration, define a CUDA option defaulting off, include
`cmake/FortSparseNativeCUDA.cmake` after creating the library target, and call
`fortsparse_add_cuda_tests()` inside the testing branch. The helper respects
`CMAKE_CUDA_ARCHITECTURES`. The two device fixtures test emitted derivatives
and resident panel/triangular/Schur operations against independent host
oracles; they require an actual CUDA device.

`fortsparse_native_cuda.F90` exposes ISO_C_BINDING interfaces and the
`native_cuda_available` build capability. CUDA interfaces are guarded by
`FORTSPARSE_HAVE_NATIVE_CUDA`; compiling the module without that definition
requires no CUDA compiler, headers or runtime. The optional CMake helper adds
the bindings and definition. The Fortran device fixture verifies pointer,
size and status interoperability through an independently known 2-by-2 factor.
The supernode fixture accepts `--cuda` to exercise sparse/block-boundary
factorization, retained-context refactorization and recovery after a failed
pivot on a real device. The CPU fixture exercises the same mathematical cases
without CUDA.
