# Bounded stability precision measurements

Build `bench_csr_precision` and `bench_mixed_lu` with
`FORTSPARSE_BUILD_BENCHMARKS=ON` in a Release build. Both use independent
manufactured-solution checks and stop on accuracy failures.

`bench_csr_precision` measures a fixed 262144-row nonsymmetric tridiagonal
operator in fp64 and scaled fp32 storage/compute. It excludes conversion/setup,
uses two warmups and seven samples of four products. Reported seconds are per
product. `bench_mixed_lu` measures setup (including conditioning estimates and
optional storage conversion) and a retained block of sixteen solves on a
600-row, half-bandwidth-nine nonsymmetric matrix. A separate warm solve checks
the setup result. FP64 factorization is used for both storage modes.

Run five fresh processes per configuration, pin threads/affinity, and retain
CSV, hardware/compiler information, binary/source hashes and completed-process
wall/user/system/peak-RSS observations. `cpu_s` in CSV measures the named
in-process phase, not completed-process lifetime CPU. `%bytes()` is retained
array payload; it is not RSS and excludes transient setup overlap.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close \
  ./build/bin/bench_csr_precision > csr.csv
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close \
  ./build/bin/bench_mixed_lu > direct.csv
```

Sweep 1/2/4/8 threads locally and 1/2/4/8/16 where resources permit remotely.
CSR rows and independent system batches parallelize. One reference LU context
remains serial; more requested threads do not parallelize its factorization.
These bounded synthetic observations do not establish consumer or external-solver
performance admission. Compressed storage may save memory while refinement
increases solve time; retain both outcomes without claiming a universal speedup.
