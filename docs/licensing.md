# Licensing of fortsparse and its optional components

The fortsparse library is MIT licensed (Copyright (c) 2025 lazy-fortran). It
links only permissive code: SuperLU (BSD), the optional SuiteSparse AMD (BSD) and
the MIT C shims and MIT AMD ordering in this repository. It contains no GPL source, includes no GPL header and
links no GPL library, so linking it statically or dynamically imposes no copyleft
obligation.

## UMFPACK helper (GPL)

`helper/fortsparse_umfpack_helper.c` is a separate executable that links
UMFPACK and is therefore GPL. The library communicates with it across a process
boundary over shared memory and named semaphores. A separate process behind
such a boundary keeps the GPL of UMFPACK confined to the helper binary; the
argument is given in [design](design.md). Building the helper is optional
(`FORTSPARSE_ENABLE_UMFPACK_IPC`); a build without it never produces GPL code.

## MUMPS (CeCILL-C)

The optional MUMPS backend links MUMPS in process. MUMPS is under CeCILL-C, a
weak-copyleft licence in the LGPL family. fortsparse contains no MUMPS source;
only the MIT shim includes its public header. MUMPS stays under CeCILL-C wherever
a build redistributes it. A build with `FORTSPARSE_ENABLE_MUMPS=OFF` carries no
such term.

## Benchmark comparator (GPL)

The optional `bench_native_families_umfpack` comparator links UMFPACK directly
in its own GPL-2.0-or-later executable. Its sources and licence are kept in
`bench/native_families/umfpack_direct/`. It does not link fortsparse.

## FreeFem plugin (LGPL)

The optional FreeFem module includes FreeFem SDK headers under their own LGPL
terms. It is a separate module, disabled by default, and adds no FreeFem code to
the core library; see the [plugin README](../plugins/freefem/README.md).
