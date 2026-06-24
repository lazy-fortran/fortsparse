# fortsparse design

This document covers the backend abstraction, the CSC storage format, the
shared-memory IPC protocol behind the UMFPACK helper, the MIT/GPL process
boundary, and the roadmap.

## Backend abstraction

One public API, many factorization backends. The solver type holds the backend
polymorphically and dispatches through deferred bindings, so adding a backend
touches neither client code nor the dispatch logic.

`sparse_backend_t` (in `fortsparse_backend`) is an abstract type with five
deferred bindings:

- `factor_real(self, A, refine, status)` for a real CSC matrix,
- `factor_complex(self, A, refine, status)` for a complex CSC matrix,
- `solve_real(self, b, x, status)` and `solve_complex(self, b, x, status)`
  for one right-hand side, reusing the retained factorization,
- `free(self)` to release the factors.

A concrete backend retains its factorization across solves and maps its native
error codes onto the `fortsparse_status` codes. The status carrier reports
`FORTSPARSE_OK`, `FORTSPARSE_SINGULAR`, `FORTSPARSE_INVALID_MATRIX`,
`FORTSPARSE_BACKEND_UNAVAILABLE`, `FORTSPARSE_NOT_FACTORED`, or
`FORTSPARSE_INTERNAL_ERROR`, each with a human-readable message.

`sparse_solver_t` (in `fortsparse_solver`) carries the selection and the
handle:

```fortran
type, public :: sparse_solver_t
    integer :: backend_id = FORTSPARSE_BACKEND_SUPERLU
    logical :: refine = .false.
    logical :: factored = .false.
    class(sparse_backend_t), allocatable :: backend
end type sparse_solver_t
```

`FORTSPARSE_BACKEND_SUPERLU = 1` is the default. `FORTSPARSE_BACKEND_UMFPACK_IPC
= 2` selects the out-of-process helper. KLU and PARDISO ids are reserved.

`sparse_factor` allocates the concrete backend for `backend_id` through a
`select case` factory, releasing any prior factorization first, then dispatches
to `backend%factor_*`. An unknown id sets `FORTSPARSE_BACKEND_UNAVAILABLE` and
leaves the handle deallocated. `sparse_solve` returns `FORTSPARSE_NOT_FACTORED`
before a factor, otherwise dispatches to `backend%solve_*`. `sparse_free` calls
`backend%free`, deallocates, and clears `factored`. `sparse_solve_once` chains
factor, solve, free for a single right-hand side.

### Adding a backend

1. Write a module under `src/<name>/` with a type that extends
   `sparse_backend_t` and implements the five bindings. Hold the retained
   factorization in the type; map the library's native status onto the
   `fortsparse_status` codes.
2. Add a public `FORTSPARSE_BACKEND_<NAME>` id in `fortsparse_solver`.
3. Add the `case` to the `ensure_backend` factory that allocates the new type.
4. Wire the sources into `src/CMakeLists.txt` and, if the backend is pure
   Fortran plus permissive C, into `fpm.toml`.
5. Add solver tests under `test/solver/` that select the new id.

A backend that links a permissive library (BSD, MIT) belongs in-process, like
SuperLU. A backend that would link a copyleft library belongs out-of-process,
behind the IPC boundary, like UMFPACK.

## CSC format

The storage type is compressed sparse column, the column-major layout the
direct solvers expect. `csc_t` holds real double values, `csc_z_t` holds
complex double. Both store `nrow`, `ncol`, `nnz`, and three arrays:

- `col_ptr(ncol+1)`: start index of each column; `col_ptr(ncol+1) = nnz+1`,
- `row_idx(nnz)`: row index of each stored entry,
- `val(nnz)`: the values.

Indices are 1-based Fortran indices. `row_idx` values lie in `[1, nrow]`, and
within each column they are sorted ascending. The backends convert to the
0-based `int` or `int64` indexing their C libraries want.

`csc_from_triplet` builds a CSC matrix from coordinate triplets `(rows, cols,
vals)`. It validates the shape, orders the triplets by `(col, row)` with two
stable counting-sort passes, and compresses, summing the values of duplicate
`(row, col)` pairs into one stored entry. `csc_is_valid` checks the structural
invariants: `col_ptr` monotone with the canonical endpoints, `row_idx` lengths
consistent, every row index in range. `csc_matvec` is the reference sparse
matrix-vector product used by the tests to verify residuals.

## Shared-memory IPC protocol

The UMFPACK backend runs the GPL solver in a separate process and drives it
over shared memory. The design target is highest performance across Linux,
macOS, and Windows: named shared memory carries the bulk data with zero copies
over a wire, and two named semaphores act as a request/done doorbell. No matrix
or vector is ever serialized over a pipe or socket. The C primitives live in
`src/ipc/fsparse_ipc.c` (MIT) behind three `#ifdef` branches, `__linux__`,
`__APPLE__`, `_WIN32`. The Fortran side is `fortsparse_umfpack_ipc`.

### Platform choices

- **Named semaphores, not unnamed.** macOS supports named POSIX semaphores
  (`sem_open`) but not unnamed semaphores placed in shared memory. Two named
  semaphores work identically on Linux and macOS, so the doorbell uses named
  semaphores everywhere POSIX. Windows uses `CreateSemaphoreA`.
- **No AF_UNIX sockets.** Unix-domain sockets are unreliable on Windows, so the
  doorbell never uses them. Shared memory plus named semaphores has one
  mechanism per platform: POSIX `shm_open`/`mmap` and `sem_open` on Linux and
  macOS, Win32 `CreateFileMappingA`/`MapViewOfFile` and `CreateSemaphoreA` on
  Windows.
- **Spawn.** `posix_spawn` on POSIX, `CreateProcessA` on Windows. The parent
  passes the shared-memory name, the two semaphore names, and the byte size as
  argv; on shutdown it terminates and reaps the child.

### Mapping layout

One shared mapping holds a header at offset 0 followed by a data region. The
header records the operation, a normalized status, the matrix shape, and the
byte offset of every array inside the same mapping, so the helper reads each
array straight from the offset with no copy. The layout is in
`include/fsparse_proto.h`, shared verbatim by the library and the helper:

```c
typedef struct fsparse_shm_header {
    volatile int32_t opcode;
    volatile int32_t status;
    int64_t n;
    int64_t nnz;
    int32_t refine;
    int32_t is_complex;
    int64_t off_colptr; /* int64 colptr[n+1], 0-based */
    int64_t off_rowidx; /* int64 rowidx[nnz], 0-based */
    int64_t off_ax;     /* double values (real, or complex real parts) */
    int64_t off_az;     /* double complex imaginary parts */
    int64_t off_b;      /* double rhs (real, or complex real parts) */
    int64_t off_bz;     /* double complex rhs imaginary parts */
    int64_t off_x;      /* double solution (real, or complex real parts) */
    int64_t off_xz;     /* double complex solution imaginary parts */
    int64_t data_bytes; /* size of the data region after the header */
} fsparse_shm_header;
```

Indices are stored 0-based to match UMFPACK's `int64` CSC convention; the
Fortran side converts its 1-based `col_ptr`/`row_idx` once while writing them
into the mapping. Complex data is split: real parts at the `*_ax`/`*_b`/`*_x`
offsets, imaginary parts at the matching `*_az`/`*_bz`/`*_xz` offsets, which is
the split form `umfpack_zl_*` takes.

The `opcode` and `status` fields are `volatile`: one side writes, the other
reads, and the semaphore doorbell is the only synchronization. The offsets and
shape are set once per factorization.

### Request cycle

The opcodes are `FACTOR_REAL`, `FACTOR_COMPLEX`, `SOLVE_REAL`, `SOLVE_COMPLEX`,
`FREE`, `SHUTDOWN`. The library:

1. sizes the mapping (header, `colptr(n+1)`, `rowidx(nnz)`, values, the
   right-hand side, the solution; the complex layout adds the imaginary
   arrays),
2. starts the session with `fsparse_ipc_start`, which creates the mapping and
   semaphores and spawns the helper,
3. writes the CSC arrays into the mapping through `c_f_pointer` views,
4. sets the header fields and opcode, then calls `fsparse_ipc_call`, which
   posts the request semaphore, waits on the done semaphore, and returns the
   status the helper wrote.

The helper loops: wait on the request semaphore, dispatch on `opcode`, post the
done semaphore. It keeps the UMFPACK `Symbolic` and `Numeric` factors resident
across solves, so a factor followed by many solves pays the factorization cost
once. It writes a small normalized status into the header: `0` ok, `1`
singular, `2` error. `SHUTDOWN` breaks the loop, and the library then unmaps,
unlinks, closes, and reaps.

### Helper discovery

`fortsparse_umfpack_ipc` finds the helper from the `FORTSPARSE_UMFPACK_HELPER`
environment variable, or by scanning `PATH` for an executable named
`fortsparse_umfpack_helper`. When no helper is found, it returns
`FORTSPARSE_BACKEND_UNAVAILABLE` with a clear message. That is the expected,
testable state under the fpm/fo build, which builds no helper. It is real
behavior, not a stub: a client that selects the UMFPACK backend without a
helper gets a clean status, never a crash.

## The MIT/GPL boundary

The split is structural, enforced by what each build sees, not by convention.

- The MIT library, `libfortsparse`, links only SuperLU (BSD) and the MIT C
  shims. It contains no UMFPACK source, includes no `umfpack.h`, and links no
  `-lumfpack`.
- UMFPACK is used only from `fortsparse_umfpack_helper`, a separate standalone
  executable that links UMFPACK. It is the only file that includes `umfpack.h`.
- The two sides communicate only across the process boundary, over shared
  memory and named semaphores. `include/fsparse_proto.h`, shared by both, is
  pure data layout and opcodes; it carries no algorithm and no UMFPACK or
  SuperLU symbol, so including it imposes no license obligation.

The structure makes the separation hard to break by accident:

- fpm and fo see only `src/`. The helper lives under `helper/`, outside any
  fpm source directory, so the fpm/fo build never compiles it and the library
  it produces is pure MIT plus BSD.
- CMake builds two distinct targets. The `fortsparse` library links
  `SuperLU::SuperLU` and `Threads::Threads`; the `fortsparse_umfpack_helper`
  executable links `UMFPACK::UMFPACK`. No CMake edge links UMFPACK into the
  library. Test executables link `fortsparse` only.

A separate process behind an IPC boundary is the FSF-recognized separation that
keeps a GPL dependency from imposing copyleft on the program that drives it.
The GPL of UMFPACK stays confined to the helper binary, and building that
binary is the user's explicit choice.

### Verifying the boundary

`ldd` and `nm` on `libfortsparse` and on the test binaries must show no UMFPACK
symbol; the helper binary must show UMFPACK. If `libfortsparse` links UMFPACK,
the build is wrong.

## Roadmap

- **Own clean-room LU.** Write a clean-room sparse LU under the MIT license to
  retire the UMFPACK dependency. The helper and its GPL boundary then become
  optional rather than the only route to UMFPACK-class fill-reducing direct
  factorization.
- **Enzyme AD.** The value path stays explicit so reverse-mode autodiff
  (planned Enzyme integration) differentiates the solve without rewriting the
  numerics. The `FORTSPARSE_ENABLE_ENZYME` CMake switch is the inert hook for
  this.
- **Device-memory sharing.** Extend the backend interface and, for the
  out-of-process path, the shared mapping so factor and solve can run on
  device memory. The `FORTSPARSE_ENABLE_DEVICE_OFFLOAD` switch is the inert
  hook.
- **More backends.** KLU and PARDISO behind the same API, by their reserved
  ids. KLU is permissive and would run in-process; a PARDISO with copyleft or
  proprietary terms would run behind the IPC boundary.
- **fortnum reuse.** Adopt the sibling library fortnum for shared kinds and
  numerical utilities once its API stabilizes.
