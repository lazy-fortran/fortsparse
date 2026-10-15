# Optional METIS ordering dependency

The native ordering adapter calls an external METIS library through its public
`METIS_NodeND` interface. The adapter is MIT project code; it contains no copied
METIS ordering implementation or vendored third-party headers. METIS remains
optional: a build without it returns `FORTSPARSE_BACKEND_UNAVAILABLE` from the
ordering adapter. This does not introduce MPI.

Primary source checked on 2026-10-16:

- Repository: [KarypisLab/METIS](https://github.com/KarypisLab/METIS).
- Revision: `272d4a91c5f66c92327493339a476c553ebf1f5d`.
- Exact upstream [license notice](https://github.com/KarypisLab/METIS/blob/272d4a91c5f66c92327493339a476c553ebf1f5d/LICENSE),
  retained as [METIS-LICENSE.txt](METIS-LICENSE.txt).
- The notice identifies Apache License 2.0. Its complete text is retained as
  [Apache-2.0.txt](Apache-2.0.txt), downloaded from the
  [Apache Software Foundation](https://www.apache.org/licenses/LICENSE-2.0.txt).

SHA-256 of the exact retained upstream notice:
`64ab947d7b289ad76e935adff51b31e4ac160df7dfad24480dbaa452e39bbe79`.
The installed `/usr/share/licenses/metis/LICENSE` used for the first local
verification has the same hash. The installed header reports METIS 5.2.1 with
32-bit `idx_t`. The wrapper converts indices through `idx_t`, so a correctly
matched 64-bit METIS header/library installation is also supported.

SHA-256 of the complete Apache license:
`cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`.

Discovery reuses a provided METIS CMake target or locates `metis.h` and
`libmetis`. It does not use SuiteSparse's separately bundled/modified METIS
source. Redistributors must retain the linked dependency's notices and any
additional notices supplied by their particular binary package.
