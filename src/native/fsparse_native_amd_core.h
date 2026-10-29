#ifndef FSPARSE_NATIVE_AMD_CORE_H
#define FSPARSE_NATIVE_AMD_CORE_H
#include <stdint.h>
/* MIT: project-owned approximate minimum degree on a zero-based symmetric
 * graph without diagonal entries. Same caller-workspace contract as AMD's
 * amd_2: pe/len/iw hold the adjacency lists, iw[pfree..iwlen) is elbow room
 * with iwlen >= pfree+n, and every array is destroyed. On return
 * last[k] = old index of the k-th pivot and next[old] = k. Returns 0, or a
 * nonzero code for an invalid argument or a broken internal invariant. */
int fsparse_amd_core(int32_t n, int32_t *pe, int32_t *iw, int32_t *len,
                     int32_t iwlen, int32_t pfree, int32_t *nv, int32_t *next,
                     int32_t *last, int32_t *head, int32_t *elen,
                     int32_t *degree, int32_t *w);
#endif
