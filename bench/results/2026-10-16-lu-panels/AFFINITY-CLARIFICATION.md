# Affinity clarification, appended after both measurement windows

CPU30 denotes Linux logical CPU ID30. On the verified AMD Ryzen9 5950X host
(16 physical cores/32 logical CPUs), it belongs to physical core14 and has SMT
sibling logicalCPU14. Commands pinned the target to CPU30; they did not reserve
that physical core or exclude unrelated work on the SMT sibling. Owned builds,
tests and profiling were held during each ordinary window, but shared-host
background services remained. No claim of32 physical cores or exclusive
physical-core isolation is made. Original at-run manifests/receipts are retained
unchanged; references to a single physical CPU30 in shorthand are clarified here.
