/* Cross-platform shared-memory + named-semaphore IPC shim (MIT).
 *
 * The library uses these primitives to drive a separate helper process at
 * arm's length: a single named shared mapping carries bulk data, two named
 * semaphores form the request/done doorbell, and the helper is spawned as a
 * child process. No data is serialized over pipes or sockets. This file and
 * its implementation contain no UMFPACK symbol; the GPL helper is a separate
 * executable reached only through this boundary.
 *
 * Platform branches live in fsparse_ipc.c (__linux__, __APPLE__, _WIN32). */

#ifndef FSPARSE_IPC_H
#define FSPARSE_IPC_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Size in bytes of the protocol header that prefixes the shared mapping. The
 * Fortran side adds the data region after this many bytes. */
int64_t fsparse_ipc_header_bytes(void);

/* Create the shared mapping and semaphores, then spawn the helper with the
 * generated names and the byte size as argv. Returns an opaque session handle,
 * or NULL with *err != 0 when the helper is missing or the spawn fails. */
void *fsparse_ipc_start(const char *helper_path, int64_t bytes, int *err);

/* Pointer to the mapped region (the protocol header is at offset 0). */
void *fsparse_ipc_data(void *sess);

/* Post the request doorbell, wait on the done doorbell, and return the status
 * the helper wrote into the protocol header. */
int fsparse_ipc_call(void *sess);

/* Tell the helper to shut down, then unmap, unlink, close, and reap it. */
void fsparse_ipc_stop(void *sess);

#ifdef __cplusplus
}
#endif

#endif /* FSPARSE_IPC_H */
