# Compressed CSR and GMRES iterative refinement

`fortsparse_precision_csr:precision_csr_t` is an additive real square CSR
interface. Existing CSC/direct-solver contexts are unchanged. It reuses
[FortNum](https://github.com/lazy-fortran/fortnum)'s precision vectors, stage
policies, refinement controller and FGMRES implementation.

`%store(n,row,column,values,storage,ok)` accepts one-based CSR, real64 values
and a storage identity from `fortnum_precision`. It validates monotonic row
pointers and column bounds before replacing the previous matrix. Duplicate
columns are summed by the ordinary SpMV accumulation. `%matvec(x,y,compute,ok)`
accepts real64 boundaries and fp32/fp64 arithmetic; fp32 arithmetic normalizes
input and matrix values before accumulation. Nonfinite output fails explicitly.
`%bytes()` counts value and index payload, excluding object and allocator
metadata. Large row loops use OpenMP. Matrix conversion allocates during setup;
SpMV allocates no workspace.

`csr_gmres_ir(a,b,x,policy,tolerance,max_steps,inner_iterations,restart,`
`info,steps,residual)` requires original **fp64** CSR storage. The original
operator is retained for every outer residual and final acceptance. A separate
rounded matrix follows `STAGE_PRECONDITION`, and is rebuilt after escalation.
Each correction uses existing restarted flexible GMRES with a selectable
low-precision operator and fp64 basis/orthogonalization, then fp64 defect
correction. This is not a low-precision direct factorization or an entirely
fp16/fp32 Krylov basis. The inner solve targets a 5% correction residual;
the caller chooses the outer tolerance. See FortNum for return statuses.

CSR fp32 computation with original fp64 storage requires entries representable
in fp32; use narrowed, scaled storage for larger magnitudes. A single matrix
scale does not preserve arbitrarily small entries in a large dynamic range.
No native BF16 dot-product, native FP16 arithmetic, CUDA backend, hard allocation
cap or performance admission is claimed by this API. Existing capped direct
solvers retain their separate contracts. Tests use an independently evaluated
nonsymmetric tridiagonal system, both arithmetic kinds, all four storage
formats and scaled coefficients far above the binary16 range.
