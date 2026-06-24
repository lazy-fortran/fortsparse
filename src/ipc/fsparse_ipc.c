/* Cross-platform shared-memory + named-semaphore IPC shim (MIT, our code).
 *
 * Drives a separate helper process at arm's length. One named shared mapping
 * carries bulk data, two named semaphores (req, done) form the doorbell, and
 * the helper runs as a spawned child. Three platform branches: POSIX shared
 * memory + named POSIX semaphores + posix_spawn on Linux and macOS, and the
 * Win32 file-mapping + semaphore + CreateProcess equivalents on Windows.
 *
 * No UMFPACK symbol appears here; the GPL helper is reached only over this
 * boundary, never linked into the library. */

#include "fsparse_ipc.h"
#include "fsparse_proto.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int64_t fsparse_ipc_header_bytes(void)
{
    return (int64_t) sizeof(fsparse_shm_header);
}

/* Directory of the current executable, without trailing slash. Writes at most
 * n bytes (NUL-terminated) into buf and returns the directory length, or 0 on
 * failure. Lets a program find a helper sitting next to it with no PATH entry
 * and no environment variable. */
#if defined(_WIN32)

#include <windows.h>

size_t fsparse_self_dir(char *buf, size_t n)
{
    DWORD len;
    size_t i;

    if (buf == NULL || n == 0) return 0;
    len = GetModuleFileNameA(NULL, buf, (DWORD) n);
    if (len == 0 || len >= n) return 0;
    for (i = len; i > 0; i--) {
        if (buf[i - 1] == '\\' || buf[i - 1] == '/') {
            buf[i - 1] = '\0';
            return i - 1;
        }
    }
    return 0;
}

#elif defined(__APPLE__)

#include <mach-o/dyld.h>

size_t fsparse_self_dir(char *buf, size_t n)
{
    uint32_t size = (uint32_t) n;
    size_t len, i;

    if (buf == NULL || n == 0) return 0;
    if (_NSGetExecutablePath(buf, &size) != 0) return 0;
    len = strlen(buf);
    for (i = len; i > 0; i--) {
        if (buf[i - 1] == '/') {
            buf[i - 1] = '\0';
            return i - 1;
        }
    }
    return 0;
}

#else /* Linux and other /proc systems */

#include <unistd.h>

size_t fsparse_self_dir(char *buf, size_t n)
{
    ssize_t len;
    size_t i;

    if (buf == NULL || n == 0) return 0;
    len = readlink("/proc/self/exe", buf, n - 1);
    if (len <= 0 || (size_t) len >= n - 1) return 0;
    buf[len] = '\0';
    for (i = (size_t) len; i > 0; i--) {
        if (buf[i - 1] == '/') {
            buf[i - 1] = '\0';
            return i - 1;
        }
    }
    return 0;
}

#endif

#if defined(_WIN32)

#include <windows.h>

typedef struct {
    HANDLE map;
    void *region;
    HANDLE sem_req;
    HANDLE sem_done;
    HANDLE process;
    int64_t bytes;
    char shm_name[64];
    char req_name[64];
    char done_name[64];
} ipc_session;

static void make_names(ipc_session *s)
{
    static long counter = 0;
    DWORD pid = GetCurrentProcessId();
    long c = ++counter;
    sprintf(s->shm_name, "fsparse_shm_%lu_%ld", (unsigned long) pid, c);
    sprintf(s->req_name, "fsparse_req_%lu_%ld", (unsigned long) pid, c);
    sprintf(s->done_name, "fsparse_done_%lu_%ld", (unsigned long) pid, c);
}

void *fsparse_ipc_start(const char *helper_path, int64_t bytes, int *err)
{
    ipc_session *s;
    char cmd[1024];
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;

    *err = 0;
    s = (ipc_session *) calloc(1, sizeof(ipc_session));
    if (s == NULL) { *err = 1; return NULL; }
    s->bytes = bytes;
    make_names(s);

    s->map = CreateFileMappingA(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
                                (DWORD) (bytes >> 32), (DWORD) bytes,
                                s->shm_name);
    if (s->map == NULL) { *err = 1; free(s); return NULL; }
    s->region = MapViewOfFile(s->map, FILE_MAP_ALL_ACCESS, 0, 0,
                              (SIZE_T) bytes);
    if (s->region == NULL) {
        CloseHandle(s->map);
        *err = 1;
        free(s);
        return NULL;
    }
    memset(s->region, 0, (size_t) bytes);

    s->sem_req = CreateSemaphoreA(NULL, 0, 1, s->req_name);
    s->sem_done = CreateSemaphoreA(NULL, 0, 1, s->done_name);
    if (s->sem_req == NULL || s->sem_done == NULL) {
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }

    {
        int need = _snprintf_s(cmd, sizeof(cmd), _TRUNCATE,
                               "\"%s\" %s %s %s %lld", helper_path, s->shm_name,
                               s->req_name, s->done_name, (long long) bytes);
        if (need < 0) {
            fsparse_ipc_stop(s);
            *err = 1;
            return NULL;
        }
    }
    memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    memset(&pi, 0, sizeof(pi));
    if (!CreateProcessA(NULL, cmd, NULL, NULL, FALSE, 0, NULL, NULL, &si,
                        &pi)) {
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }
    CloseHandle(pi.hThread);
    s->process = pi.hProcess;

    /* No early-unlink is needed or possible here. Win32 named mappings and
     * semaphores are reference-counted by open handle and have no persistent
     * filesystem entry (unlike POSIX /dev/shm and named semaphores). The kernel
     * closes every handle when a process terminates, including on a hard kill,
     * so the objects disappear automatically once both processes are gone. The
     * POSIX early-unlink-after-attach handshake exists only to drop the
     * filesystem names that would otherwise survive a SIGKILL; on Windows there
     * is nothing to drop, so this branch keeps its handles for the whole
     * session and leaks nothing on a hard kill. */
    return s;
}

void *fsparse_ipc_data(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    return s == NULL ? NULL : s->region;
}

int fsparse_ipc_call(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    fsparse_shm_header *h = (fsparse_shm_header *) s->region;
    HANDLE waits[2];
    DWORD w;
    ReleaseSemaphore(s->sem_req, 1, NULL);
    /* Wake on the done doorbell, or on the helper process exiting first, so a
     * crashed helper returns an error status instead of hanging forever. */
    waits[0] = s->sem_done;
    waits[1] = s->process;
    w = WaitForMultipleObjects(2, waits, FALSE, INFINITE);
    if (w != WAIT_OBJECT_0) {
        h->status = FSPARSE_ST_ERROR;
        return FSPARSE_ST_ERROR;
    }
    return (int) h->status;
}

void fsparse_ipc_stop(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    if (s == NULL) return;
    if (s->process != NULL && s->region != NULL) {
        fsparse_shm_header *h = (fsparse_shm_header *) s->region;
        h->opcode = FSPARSE_OP_SHUTDOWN;
        if (s->sem_req != NULL) ReleaseSemaphore(s->sem_req, 1, NULL);
        WaitForSingleObject(s->process, 2000);
    }
    if (s->process != NULL) {
        TerminateProcess(s->process, 0);
        CloseHandle(s->process);
    }
    if (s->sem_req != NULL) CloseHandle(s->sem_req);
    if (s->sem_done != NULL) CloseHandle(s->sem_done);
    if (s->region != NULL) UnmapViewOfFile(s->region);
    if (s->map != NULL) CloseHandle(s->map);
    free(s);
}

#else /* POSIX: Linux and macOS */

#include <errno.h>
#include <fcntl.h>
#include <semaphore.h>
#include <spawn.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

extern char **environ;

typedef struct {
    int fd;
    void *region;
    sem_t *sem_req;
    sem_t *sem_done;
    pid_t pid;
    int64_t bytes;
    char shm_name[64];
    char req_name[64];
    char done_name[64];
} ipc_session;

static void make_names(ipc_session *s)
{
    static long counter = 0;
    long pid = (long) getpid();
    long c = ++counter;
    sprintf(s->shm_name, "/fsparse_shm_%ld_%ld", pid, c);
    sprintf(s->req_name, "/fsparse_req_%ld_%ld", pid, c);
    sprintf(s->done_name, "/fsparse_done_%ld_%ld", pid, c);
}

void *fsparse_ipc_start(const char *helper_path, int64_t bytes, int *err)
{
    ipc_session *s;
    char *argv[6];
    char bytes_arg[32];
    int rc;

    *err = 0;
    s = (ipc_session *) calloc(1, sizeof(ipc_session));
    if (s == NULL) { *err = 1; return NULL; }
    s->bytes = bytes;
    s->fd = -1;
    s->sem_req = SEM_FAILED;
    s->sem_done = SEM_FAILED;
    s->pid = -1;
    make_names(s);

    s->fd = shm_open(s->shm_name, O_CREAT | O_RDWR | O_EXCL, 0600);
    if (s->fd < 0) { *err = 1; free(s); return NULL; }
    if (ftruncate(s->fd, (off_t) bytes) != 0) {
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }
    s->region = mmap(NULL, (size_t) bytes, PROT_READ | PROT_WRITE, MAP_SHARED,
                     s->fd, 0);
    if (s->region == MAP_FAILED) {
        s->region = NULL;
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }
    memset(s->region, 0, (size_t) bytes);

    s->sem_req = sem_open(s->req_name, O_CREAT | O_EXCL, 0600, 0);
    s->sem_done = sem_open(s->done_name, O_CREAT | O_EXCL, 0600, 0);
    if (s->sem_req == SEM_FAILED || s->sem_done == SEM_FAILED) {
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }

    sprintf(bytes_arg, "%lld", (long long) bytes);
    argv[0] = (char *) helper_path;
    argv[1] = s->shm_name;
    argv[2] = s->req_name;
    argv[3] = s->done_name;
    argv[4] = bytes_arg;
    argv[5] = NULL;
    rc = posix_spawn(&s->pid, helper_path, NULL, NULL, argv, environ);
    if (rc != 0) {
        s->pid = -1;
        fsparse_ipc_stop(s);
        *err = 1;
        return NULL;
    }

    /* Await the helper's READY post on the done semaphore, then unlink all
     * three names immediately. The mapping and the open semaphore handles stay
     * valid in both processes (the kernel keeps each object alive until its
     * last reference closes), but the names are gone from /dev/shm, so a hard
     * kill (SLURM/HTCondor SIGKILL on timeout or OOM) of either process leaves
     * no residue on the node. The READY post is consumed here, before any
     * request, so the per-operation req/done doorbell stays in step. */
    while (sem_wait(s->sem_done) != 0) {
        if (errno != EINTR) {
            fsparse_ipc_stop(s);
            *err = 1;
            return NULL;
        }
    }
    sem_unlink(s->req_name);
    sem_unlink(s->done_name);
    shm_unlink(s->shm_name);
    return s;
}

void *fsparse_ipc_data(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    return s == NULL ? NULL : s->region;
}

int fsparse_ipc_call(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    fsparse_shm_header *h = (fsparse_shm_header *) s->region;
    sem_post(s->sem_req);
    /* Wait on the done doorbell, but poll for helper death so a crashed or
     * killed helper returns an error status instead of hanging forever. */
    for (;;) {
        struct timespec ts;
        int wstatus;
        pid_t r;
        clock_gettime(CLOCK_REALTIME, &ts);
        ts.tv_nsec += 100000000L; /* 100 ms */
        if (ts.tv_nsec >= 1000000000L) {
            ts.tv_nsec -= 1000000000L;
            ts.tv_sec += 1;
        }
        if (sem_timedwait(s->sem_done, &ts) == 0) return (int) h->status;
        if (errno == EINTR) continue;
        if (errno != ETIMEDOUT) return FSPARSE_ST_ERROR;
        if (s->pid <= 0) continue;
        r = waitpid(s->pid, &wstatus, WNOHANG);
        if (r == s->pid || (r < 0 && errno == ECHILD)) {
            s->pid = -1; /* reaped: keep fsparse_ipc_stop from blocking */
            h->status = FSPARSE_ST_ERROR;
            return FSPARSE_ST_ERROR;
        }
    }
}

void fsparse_ipc_stop(void *sess)
{
    ipc_session *s = (ipc_session *) sess;
    int wstatus;
    if (s == NULL) return;
    if (s->pid > 0 && s->region != NULL && s->sem_req != SEM_FAILED) {
        fsparse_shm_header *h = (fsparse_shm_header *) s->region;
        h->opcode = FSPARSE_OP_SHUTDOWN;
        sem_post(s->sem_req);
        waitpid(s->pid, &wstatus, 0);
        s->pid = -1;
    }
    if (s->pid > 0) waitpid(s->pid, &wstatus, 0);
    if (s->sem_req != SEM_FAILED) sem_close(s->sem_req);
    if (s->sem_done != SEM_FAILED) sem_close(s->sem_done);
    if (s->region != NULL) munmap(s->region, (size_t) s->bytes);
    if (s->fd >= 0) close(s->fd);
    /* In the normal path fsparse_ipc_start already unlinked all three names
     * once the helper attached, so these calls just return ENOENT. They remain
     * for the early-failure path (a shm_open/mmap/sem_open error before the
     * READY handshake), where the names still exist and must be removed.
     * Unlinking an already-unlinked name is harmless, so stop stays idempotent
     * and never leaves residue. */
    sem_unlink(s->req_name);
    sem_unlink(s->done_name);
    shm_unlink(s->shm_name);
    free(s);
}

#endif

/* Pool of process-persistent helper sessions. Spawning the helper reloads
 * UMFPACK and its BLAS every time, which dominates a factor/solve/free loop when
 * it happens per factorization; keeping helpers resident and reusing their
 * mappings removes that cost. A session is held (busy) from factor to free, so
 * several live factorizations at once (NEO-2 drives a real and a complex solver)
 * each get their own session and never cross. Idle sessions are reused while
 * their mapping is large enough, and grown (respawned once) when it is not.
 *
 * The pool is per process and driven serially: fortsparse routes every solve
 * through a shared module-global solver, so concurrent factorizations come from
 * separate MPI rank processes, each with their own pool, not from threads
 * sharing this one. Overflow past the pool size falls back to an unpooled
 * session that is torn down on release. */
#define FSPARSE_POOL_SIZE 8

static struct {
    void *sess;
    int64_t bytes;
    int busy;
} g_pool[FSPARSE_POOL_SIZE];
static int g_pool_atexit = 0;

static void fsparse_ipc_atexit(void)
{
    int i;
    for (i = 0; i < FSPARSE_POOL_SIZE; i++) {
        void *s = g_pool[i].sess;
        g_pool[i].sess = NULL;
        g_pool[i].bytes = 0;
        g_pool[i].busy = 0;
        if (s != NULL) fsparse_ipc_stop(s);
    }
}

void *fsparse_ipc_acquire(const char *helper_path, int64_t bytes, int *err)
{
    int i, slot = -1;
    int64_t want;

    *err = 0;

    /* Reuse an idle session whose mapping already fits. */
    for (i = 0; i < FSPARSE_POOL_SIZE; i++) {
        if (g_pool[i].sess != NULL && !g_pool[i].busy
                && g_pool[i].bytes >= bytes) {
            g_pool[i].busy = 1;
            return g_pool[i].sess;
        }
    }
    /* Otherwise take an empty slot, or an idle one to grow (respawn larger). */
    for (i = 0; i < FSPARSE_POOL_SIZE; i++) {
        if (g_pool[i].sess == NULL) { slot = i; break; }
    }
    if (slot < 0) {
        for (i = 0; i < FSPARSE_POOL_SIZE; i++) {
            if (!g_pool[i].busy) { slot = i; break; }
        }
    }

    /* A quarter of headroom absorbs the small per-factorization variation in
     * nonzero count so a steady-state problem size respawns at most once. */
    want = bytes + bytes / 4;

    /* Every slot busy: hand back an unpooled session for this overflow; release
     * tears it down since it is not found in the pool. */
    if (slot < 0) return fsparse_ipc_start(helper_path, want, err);

    if (g_pool[slot].sess != NULL) {
        void *old = g_pool[slot].sess;
        g_pool[slot].sess = NULL;
        g_pool[slot].bytes = 0;
        fsparse_ipc_stop(old);
    }
    g_pool[slot].sess = fsparse_ipc_start(helper_path, want, err);
    if (g_pool[slot].sess == NULL) return NULL;
    g_pool[slot].bytes = want;
    g_pool[slot].busy = 1;
    if (!g_pool_atexit) {
        atexit(fsparse_ipc_atexit);
        g_pool_atexit = 1;
    }
    return g_pool[slot].sess;
}

void fsparse_ipc_release(void *sess)
{
    int i;
    if (sess == NULL) return;
    /* A pooled session goes back to idle, keeping its helper resident for the
     * next factorization. A session not in the pool is an overflow fallback and
     * is torn down here. */
    for (i = 0; i < FSPARSE_POOL_SIZE; i++) {
        if (g_pool[i].sess == sess) { g_pool[i].busy = 0; return; }
    }
    fsparse_ipc_stop(sess);
}
