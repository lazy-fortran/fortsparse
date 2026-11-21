# Nested-dissection ordering provenance

`src/native/fortsparse_nd.f90`, `fortsparse_nd_graph.f90`,
`fortsparse_nd_refine.f90` and `fortsparse_nd_util.f90` are a Fortran
translation of the `METIS_NodeND` path of METIS 5.2.1 (`ometis.c`, `mmd.c`,
`compress.c`, `coarsen.c`, `bucketsort.c`, `graph.c`, `initpart.c`,
`balance.c`, `fm.c`, `refine.c`, `separator.c`, `srefine.c`, `sfm.c`) and of
the GKlib random numbers, sort and priority queue it uses, for the default
options, one vertex weight and seed 0:

- METIS: [KarypisLab/METIS](https://github.com/KarypisLab/METIS) tag `v5.2.1`,
  revision `f5ae915a84d3bbf1508b529af90292dd7085b9ec`. Its
  [license notice](METIS-LICENSE.txt) is retained without edits (SHA-256
  `64ab947d7b289ad76e935adff51b31e4ac160df7dfad24480dbaa452e39bbe79`).
- GKlib: [KarypisLab/GKlib](https://github.com/KarypisLab/GKlib) tag
  `METIS-v5.1.1-DistDGL-0.5`, revision
  `3eabb216ac97e11ce7e7a9b90f4c90778d9e7c18`, the version METIS 5.2.1 builds
  with. Its [license notice](GKlib-LICENSE.txt) is retained without edits
  (SHA-256 `81118d7ac02c98ad2d49a8b7d12d24f41e3bef962342f5483bf1b80a79ce1992`).

Both notices grant the Apache License 2.0, whose complete text is retained as
[Apache-2.0.txt](Apache-2.0.txt) (SHA-256
`cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`), downloaded
from the [Apache Software Foundation](https://www.apache.org/licenses/LICENSE-2.0.txt).
Each translated file names its origin and states that it is a translation.
The changes are the language, structured control flow in place of `goto`,
array kernels for speed, and the C library's `rand` (glibc's algorithm, which
METIS uses when GKlib is built without `USE_GKRAND`) written out in Fortran.

`test_fortsparse_nd` checks that the orderings equal those of METIS 5.2.1
(Arch Linux package `metis 5.2.1-1`), recorded in
`test/data/nd_metis_5.2.1.txt` for synthetic graphs.

`src/native/fortsparse_native_nd.f90` (MIT) builds the graph from one-based
CSC and exposes the ordering as `fsparse_native_order` and
`fsparse_native_order_general`.

## Optional METIS reference

With `FORTSPARSE_ENABLE_METIS=ON`, `src/native/fsparse_metis_order.c` (MIT)
calls an external METIS library through its public `METIS_NodeND` interface as
`fsparse_metis_order`, and `test_fortsparse_nd` compares the native orderings
with it; with `FORTSPARSE_ND_WRITE` set it rewrites the recording. Discovery
reuses a provided METIS CMake target or locates `metis.h` and `libmetis`.
Redistributors of such a build must retain the linked library's notices.
