/* Shared memory layout and opcodes for the fortsparse out-of-process backend.
 *
 * This header is included by both the MIT library IPC shim and the GPL UMFPACK
 * helper. It is pure data layout plus opcodes; it carries no algorithm and no
 * UMFPACK or SuperLU symbol, so including it imposes no license obligation.
 *
 * One shared mapping holds a header followed by a data region. The header
 * records the operation, a normalized status, the matrix shape, and the byte
 * offset of each array inside the same mapping. Indices are stored 0-based to
 * match UMFPACK's int64 CSC convention. Complex data is stored split: real
 * parts at the *_ax / *_b / *_x offsets and imaginary parts at the matching
 * *_az / *_bz / *_xz offsets. */

#ifndef FSPARSE_PROTO_H
#define FSPARSE_PROTO_H

#include <stdint.h>

/* Request opcodes written by the library into header->opcode before a call. */
enum {
    FSPARSE_OP_FACTOR_REAL = 1,
    FSPARSE_OP_FACTOR_COMPLEX = 2,
    FSPARSE_OP_SOLVE_REAL = 3,
    FSPARSE_OP_SOLVE_COMPLEX = 4,
    FSPARSE_OP_FREE = 5,
    FSPARSE_OP_SHUTDOWN = 6
};

/* Normalized status codes written by the helper into header->status. */
enum {
    FSPARSE_ST_OK = 0,
    FSPARSE_ST_SINGULAR = 1,
    FSPARSE_ST_ERROR = 2
};

/* Header at offset 0 of the shared mapping. The volatile request/response
 * fields are written by one side and read by the other across the semaphore
 * doorbell; the offsets and shape are set once per factorization. */
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

#endif /* FSPARSE_PROTO_H */
