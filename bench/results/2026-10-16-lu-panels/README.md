# Experimental real LU target panels: three matched comparisons

All 60 requested ordinary children passed the independent original-COO
real128 residual gates; none timed out. Width 8 remains experimental. On the
captured 5,525-row matrix, its median complete-child CPU time exceeded the
matched scalar control in all three measurement windows. The contiguous-run
schedule lost in every captured-case pair. These results do not admit a public
default, consumer pin or external-library performance claim.

| Generated width-8 schedule | Captured CPU, scalar / panel (ms) | Captured factor, scalar / panel (ms) | Small drift CPU, scalar / panel (ms) |
| --- | ---: | ---: | ---: |
| Original arbitrary mask | 523.081 / 637.687 | 364.001 / 477.657 | 11.847 / 11.008 |
| Uniform full-active SIMD | 557.854 / 602.418 | 386.453 / 425.303 | 15.169 / 10.732 |
| Contiguous active-run SIMD | 496.523 / 518.272 | 341.817 / 359.317 | 9.320 / 8.902 |

Use `trials.csv` for every individual observation and `summary.csv` for
median/minimum/maximum. Separate windows have different controls and background
conditions: comparing their panel medians directly does not measure the effect
of a schedule change. Small-case ranges overlap. The final small-case factor
median was 2.560 / 2.488 ms, with four of five pairs faster and one slower.

The captured factor fill was 3,334,519 for both widths throughout; small-case
fill was 65,598. Panel storage adds exactly `68*n+64` requested owned bytes:
375,764 bytes and 69,696 bytes respectively. Current/peak owned payload counts
factor arrays, including allocation overlap, and excludes caller/driver/check
arrays. It is distinct from the same child's process RSS.

## Reproduce the plots from this portable directory

The CSV files and `plot.py` suffice; no matrix corpus or solver build is needed.
Python and Matplotlib are optional visualization dependencies, not solver
dependencies:

```sh
python plot.py trials.csv /tmp/fortsparse-panel-plots
```

The script draws all requested points and paired lines, without filtering,
smoothing, confidence fitting or external frontier construction. It uses an
Okabe–Ito palette plus redundant markers. Complete CPU is completed-child
`wait4` user plus system time; RSS is that same child's `ru_maxrss`. CPU/RSS
panels use common scales within each case across windows; ratio panels share
the same vertical scale. Rendered PNG/SVG files are ephemeral and kept outside
source history. Both figures were rendered and visually inspected.

`prepare.py` rebuilds these tables from the three sibling raw archives; that
audit command needs the original archives. `provenance.json` retains source,
ELF, input/order, observer and resolved library hashes. The raw matrices,
binary orders, performance recordings and solver binaries stay outside Git.

## Frozen ordinary protocol and scope

Each window contains exactly two fixed cases, five fresh alternating AB/BA
pairs per case, and widths 1/8: 20 ordinary children. The scalar control ELF is
unchanged across all windows (SHA256
`0b030e20853a9f33ab0719dab872a8e354df2120c1def0710082bfc76f103fd6`).
The first window uses prototype source `2926aad`; subsequent width-8 sources
are `0598590` and `bfc607b`. Public defaults and complex execution are unchanged.
These patches are published on main as `2152242`, `5215cbc` and `8d7f634`
respectively; the measurements retain their exact original worker revisions
and frozen binaries. Full mappings are in `provenance.json`.

The captured case uses its frozen unmatched ND order and no row preference;
the public drift-real level32 training fixture uses natural order. Both use
pivot threshold 0.1 and pruning. The probe constructs a known solution,
accumulates its RHS and original-COO residual in real128, and accepts both
normwise and componentwise backward errors at most `1e-12`. Actual orders are
exported and checked bijective, with the same order hash for both widths.

Complete-child scope includes input/order reading, CSC canonicalization,
known-RHS formation, factor, solve and independent checking. Factor/solve wall
phases are reported separately. The unchanged qualified C `posix_spawn/wait4`
observer has SHA256
`e3e0692499cc35e0863bae67aab8476345ee2c126ee8980b1ca75cfafcf5ebda`;
its CLI is `MEASURE_CPUS=30 measure TARGET LOG [ARGS]`. Python supervises this
observer; its own RSS is not inherited as a target measurement floor.

Targets use Linux logical CPU30 affinity, OMP/OpenBLAS/MKL one thread, no
`LD_LIBRARY_PATH` override, a 3-GiB address-space limit and a 0.9-second GNU
timeout (TERM, two-second kill grace). Outer timeout status and child exit are
retained separately. There are no retries or shrunken cases. The Ryzen 5950X
has 16 physical cores / 32 logical CPUs; CPU30 is core14's logical ID, with
SMT sibling14. Owned jobs were held during each window, but unrelated shared
host services and SMT siblings were not reserved or stopped. See the appended
affinity clarification; at-run receipts were not rewritten.

## Separate diagnostics, not ordinary timing

`diagnostic-masks.csv` retains the first experiment's exact width-8 call/mask
counts. It contains 91,805,117 row entries and 593,122,114 scalar-equivalent
active updates, compared with 592,492,197 scalar updates. Full-eight masks
cover 56.39% of row entries. Among nonfull entries, a single contiguous active
bit run covers 60.87%; this motivated the final schedule without a new capture.
The two counter children were excluded from ordinary timing.

`diagnostic-phases.csv` retains the two later sampling profiles, one child per
width, 1,999-Hz user CPU-clock sampling with DWARF callgraphs. These are
approximate period-weighted attribution, not a controlled timing comparison.
Numerical updates account for about 67% of scalar and 69% of panel factor
samples. Scalar `spd_update_value` attribution is inline DWARF: actual scalar
scatter assembly contains no per-entry call. Masked-loop samples motivated
the contiguous-run path; actual Release assembly confirms SIMD for full/run
paths, not end-to-end speedup.

The corrected generator stack is recorded separately in
`/var/tmp/fortsparse-panel8-codegen/ACTIVE-DEPENDENCIES.md`. It depends on
FortSym `7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6`, FortGen
`42a4c6cce9272a6bc92e0539b21f3ffc0c97a068`, and the shared
`/var/tmp/neo2-curl-gen-build` include/library closure. Both retained generators
were run twice with identical checked-in outputs. Current Debug 11/11 and
Release 3/3 affected tests passed, including all 2,048 lane/mask cases,
independent implicit derivatives and exact allocator cap/fault oracles.
