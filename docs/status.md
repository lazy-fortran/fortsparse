# Implementation status

FortSparse is experimental research software. Its tests exercise the matrix
classes and operations below; application-scale accuracy, memory and speed
still require checks on the application's original matrix and right-hand side.
A successful factorization alone does not certify a physical observable.

| Path | Independent checks in the test suite |
| --- | --- |
| Real and complex solves | Known small systems, manufactured Poisson systems and original-matrix residuals |
| Retained factors | Repeated right-hand sides, refactorization, context copies and explicit destruction |
| Native positive-definite factors | Real Cholesky and complex Hermitian systems, invalid input and nonpositive pivots |
| Native general and indefinite factors | Known solutions, transpose/adjoint actions, pivoting and inertia cases |
| Implicit solve derivatives | Central finite differences and real adjoint pairings for real and complex matrices |
| Memory controls | Solver-owned allocation counters, planned layouts and rejection before exceeding a configured cap |
| Optional external backends | Backend-specific solve, status and interprocess round-trip tests when built |

See [native matrix classes](native-cpu.md) for ownership and compiler
limitations, [precision](precision.md) for refinement admission and
[memory](memory.md) for the difference between an exact layout, a bound and an
estimate. Tests that need an optional backend do not verify it when that backend
is absent.

The native CUDA path covers selected real positive-definite systems and remains
experimental. Unsupported matrix classes, device combinations and backends
return a status; they are not general GPU solver support. Native distributed
factorization, multi-GPU execution and disk-backed factors are outside the
implemented scope. Whole-application parity and performance are separate from
the library regression tests.
