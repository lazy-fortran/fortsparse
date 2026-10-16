# LU triangular solve accumulation

General native LU keeps its factorization, row equilibration, pivot policy and
requested column ordering. Normal solves accumulate triangular updates with
Kahan correction state. Transpose and adjoint gathers use the same ordered
scalar stages, with Cartesian compensation for complex values. This improves
the arithmetic used to solve the retained triangular systems; it does not
perform original-matrix iterative refinement or change `refine = .false.`.

The authoritative scalar stages and indexed schedules are retained in
`gen/gen_lu.f90` and `gen/lu_compensated_schedule.f90.inc`. FortSym emits each
addition, subtraction and product stage independently with CSE disabled. The
stages must not be composed and simplified over exact reals: doing so would
erase the rounding-error state. The generated LU arithmetic has reassociation
and contraction disabled for GNU compilers. The shared SPD leaves are unchanged.

For an incoming update `increment`, the ordered stages are
`adjusted = increment - correction`, `total = value + adjusted`,
`lost = total - value`, and `correction = lost - adjusted`, followed by
`value = total`. A completed component uses `value - correction` before its
triangular division or subsequent scatter. These are floating-point schedules,
not a certified enclosure or an exact-summation claim. Equation-level JVP/VJP
contracts still differentiate the original solve equation, and the existing
independent derivative oracles exercise the changed solve path.

Each factor reserves one correction vector before numeric factorization: `8*n`
bytes for real and `16*n` for complex factors. Fixed-workspace cap checks,
retained/peak accounting, partial-allocation cleanup, same-pattern reuse and
free include this payload. Solves allocate no memory. The allocation observer
and cap/fault-injection tests check those properties independently.

The public exact-dyadic arrow regression has 64 unit contributions and a `1e16`
contribution. Its exact stored right-hand side fixes a small first component
that the former sequential triangular scatter lost. It covers real/complex
normal, transpose and adjoint blocks, repeated solves and fixed-memory limits.
This observable test supplements the existing dense, sparse, pivot, singular,
reuse, derivative and allocation-failure tests.

A private captured NEO-2 matrix also isolated triangular accumulation from
factorization error by retaining identical rounded factors. Compensation
reduced discrepancies against a higher-precision residual-corrected reference;
that floating reference is not a certificate. The original final transport
checks improve, while the whole application golden comparison still reports
intermediate and metadata differences. No tolerance has been weakened and no
full application or performance admission follows from this change.
