# Optional AMD ordering dependency

The MIT adapter calls AMD's public `amd_order` API through a separately linked
library. It contains no AMD algorithm source or vendored third-party headers.
AMD is optional; an unavailable build returns status 3 from the adapter.

Primary source retained on 2026-10-16:

- Local research clone: `~/code/sparse-dev/repos/SuiteSparse`.
- Revision: `34db4dfbe21c53cea7b0e3c5ee2da37bfaef369a`.
- Exact [AMD license](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/34db4dfbe21c53cea7b0e3c5ee2da37bfaef369a/AMD/Doc/License.txt),
  retained without edits as [AMD-LICENSE.txt](AMD-LICENSE.txt).
- SHA-256:
  `e50fdd9a587f4860ce1d503b94576bfff47bccc421f758d1bdc8706c8cd36108`.

The retained notice specifies BSD 3-clause. The installed AMD 3.3.4 header
also identifies BSD-3-Clause. Its SuiteSparse_config dependency identifies
BSD-3-Clause in the archived `SuiteSparse_config/README.txt`; no GPL UMFPACK
code is linked into the native solver by this adapter. CMake reuses AMD's
package target, including transitive dependencies, when available. A raw
static AMD library also needs its matching SuiteSparse_config library.

The adapter accepts full canonical symmetric one-based CSC and returns
`perm[new] = old`, `inverse[old] = new`, both one-based. Input validation is
shared with the METIS adapter. It uses AMD defaults and checks the returned
permutation before exposing it. Graph ordering remains outside numerical
kernel generation.
