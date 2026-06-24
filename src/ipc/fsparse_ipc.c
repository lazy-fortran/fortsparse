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
    sem_unlink(s->req_name);
    sem_unlink(s->done_name);
    shm_unlink(s->shm_name);
    free(s);
}

#endif
