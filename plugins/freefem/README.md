# Optional native FreeFem adapter

The `fortsparse_freefem` CMake MODULE target produces `fortsparse_native.so`.
FreeFem `load "fortsparse_native"` registers `FORTSPARSE` as the real and complex
`sparsesolver`/`sparsesolversym` default. It calls the public native C ABI;
assembly, finite elements, boundary conditions, field storage and launching stay
with the application and FreeFem. No solver kernel is copied into the plugin.

This component is OFF by default and introduces no FreeFem or C++ dependency to
the ordinary core build. Enable it explicitly with a matching FreeFem SDK/runtime:

```sh
cmake -S . -B build-freefem -G Ninja \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_AMD=ON \
  -DFORTSPARSE_BUILD_FREEFEM_PLUGIN=ON \
  -DFREEFEM_CXX_EXECUTABLE=/prefix/bin/ff-c++ \
  -DFREEFEM_EXECUTABLE=/prefix/bin/FreeFem++
cmake --build build-freefem --target fortsparse_freefem -j2
ctest --test-dir build-freefem -R '^freefem_native_' --output-on-failure -j2
```

AMD is an optional BSD graph-ordering dependency of the core and is required by
this adapter's concrete Maxwell consumers. Its absence is a configuration error,
with no natural-order or external-factor fallback. The module requires in-process
SuperLU/MUMPS backends disabled so its binary closure contains no external sparse
factorization. The recipe also disables the unrelated UMFPACK helper. It preserves
the native C ABI's default refinement policy; vendor iteration/pivot controls are
not translated. This component adds no cap policy or performance admission claim.

`ff-c++ -n -c` supplies the installed SDK's include paths and ABI definitions.
Its vendor link flags are not used. A private build-directory `config.h` includes
the SDK configuration and excludes unused UMFPACK/CHOLMOD umbrella-header
adapters. The installed SDK is untouched. `FREEFEM_CXX_EXECUTABLE` and
`FREEFEM_EXECUTABLE` are cache overrides; changing prefixes refreshes the derived
SDK include path. SDK and runtime must match. Linux FreeFem 4.17 is the initial
verification target; other runtime/SDK versions require their own checks. The
installed 4.17.0-2 SDK/runtime passed all six focused CTests on 2026-10-16,
including the module linkage gate. Receipts are retained under
`/home/ert/code/sparse-dev/consumer-migrations/freefem-plugin/`; the stock
interpreter remains vendor-linked as distinguished below. The separately
provisioned minimal serial runtime/SDK from FreeFem v4.17 commit
`754dae8ae8f858784ad46b2f46dfd7862829cb5a` also passes all six tests. Separate
`ldd` receipts verify its interpreter, BinaryIO and rebuilt module each exclude
UMFPACK/CHOLMOD/SuperLU/MUMPS/PETSc. This establishes that named runtime closure,
not either application's complete scientific or performance qualification.

Applications retain their launchers. Depend on `fortsparse_freefem`, obtain its
path with `$<TARGET_FILE:fortsparse_freefem>` and its directory with
`$<TARGET_FILE_DIR:fortsparse_freefem>`, and prepend that directory to `FF_LOADPATH`
before invoking the matching interpreter. FreeFem uses semicolons between loader
search paths. Hard-load `fortsparse_native` in the application script; a missing
plugin, missing AMD, invalid input, singular factor or rank count greater than
one is an explicit failure. Loading an unrelated solver plugin later can change
FreeFem's global defaults; the application owns its load sequence.

The adapter expands FreeFem half storage through `HashMatrix::UnHalf`, exports
canonical zero-based CSC and refreshes pointers after SDK updates. Real transpose
flag 1 maps to T; complex flag 1 maps to H. That distinction follows FreeFem
4.17 `VirtualSolverSparseSuite.hpp` (complex trans selects `UMFPACK_At`, whereas
`UMFPACK_Aat` is the plain transpose). The
[upstream HashMatrix implementation](https://github.com/FreeFem/FreeFem-sources/blob/v4.17/src/femlib/HashMatrix.cpp)
also defines the conjugating half-storage expansion. Repeated RHS reuse factors;
numeric and structural changes invalidate the appropriate SDK state. The context
owns its factors and is freed with the SDK solver object.

Registered oracles use literal original real/complex matrices and known vectors,
N/T or N/H solves, repeated RHS, numeric/graph updates, Hermitian half expansion,
original-matrix residuals, repeated object lifetimes and singular failures. The
small complex Hrot/RT0Ortho oracle exercises the public NEO-EQ Ampere carrier and
bilinear form on positive cylindrical radius. These are numerical checks, not
rigorous certificates. The Linux module linkage test checks DT_NEEDED closure
and unresolved symbols for vendor sparse factorization references.

**Module linkage and interpreter linkage are separate claims.** A stock FreeFem
interpreter can preload UMFPACK/MUMPS/SuperLU even when this module never uses or
links them. The module test cannot establish full runtime independence. A minimal
serial interpreter, its loaded plugins and the launcher must be audited together
for that stronger claim; SDK/compiler/runtime/source identities belong to that
qualification. No ordinary timing or memory claim follows from these smoke tests.

The authored adapter source is MIT, like FortSparse. SDK headers and template
material retain FreeFem's licensing terms and are used only in this optional
module, never in the core target. The installed 4.17.0-2 package supplies an
LGPL-3.0-or-later license file; `ufunction.hpp` states LGPL-2.1-or-later. Preserve
the matching SDK notices and applicable terms when redistributing a combined
module; the adapter's MIT notice does not replace them. No FreeFem source is
vendored into the MIT core. Local provenance snapshots/hashes are retained under
`/home/ert/code/sparse-dev/provenance/freefem_plugin/` for the investigated SDK.
