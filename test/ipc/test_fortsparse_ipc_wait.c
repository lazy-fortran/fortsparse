/* Exercise the public IPC boundary without a numerical solver dependency. */
#include "fsparse_ipc.h"
#include "fsparse_proto.h"

#include <errno.h>
#include <fcntl.h>
#include <semaphore.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static volatile sig_atomic_t signals_received;

static void interrupt_wait(int sig)
{
    (void) sig;
    signals_received++;
    /* Bound even a broken implementation's run, then let cleanup proceed. */
    if (signals_received >= 100) alarm(0);
}

static int set_timer(int repeat)
{
    struct itimerval timer = {{0, 0}, {0, 20000}};
    if (repeat) timer.it_interval = timer.it_value;
    signals_received = 0;
    return setitimer(ITIMER_REAL, &timer, NULL);
}

static int helper(int argc, char **argv)
{
    int fd;
    size_t bytes;
    sem_t *request, *done;
    fsparse_shm_header *header;

    if (argc != 5) return 1;
    /* A killed test parent must not leave a helper blocked on its doorbell. */
    signal(SIGALRM, SIG_DFL);
    alarm(5);
    bytes = (size_t) strtoll(argv[4], NULL, 10);
    fd = shm_open(argv[1], O_RDWR, 0600);
    if (fd < 0) return 1;
    header = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if (header == MAP_FAILED) return 1;
    request = sem_open(argv[2], 0);
    done = sem_open(argv[3], 0);
    if (request == SEM_FAILED || done == SEM_FAILED) return 1;
    if (sem_post(done) != 0) return 1; /* READY */
    for (;;) {
        while (sem_wait(request) != 0) {
            if (errno != EINTR) return 1;
        }
        if (header->opcode == FSPARSE_OP_SHUTDOWN) break;
        if (header->opcode == 4) _exit(7); /* Die without a done post. */
        if (header->opcode != 1) {
            struct timespec delay = {0, 250000000L};
            while (nanosleep(&delay, &delay) != 0) {
                if (errno != EINTR) return 1;
            }
        }
        header->status = FSPARSE_ST_SINGULAR;
        if (sem_post(done) != 0) return 1;
    }
    sem_close(request);
    sem_close(done);
    munmap(header, bytes);
    return 0;
}

int main(int argc, char **argv)
{
    void *session;
    fsparse_shm_header *header;
    struct sigaction action = {0};
    int err, status, failures = 0;

    if (argc != 1) return helper(argc, argv);
    action.sa_handler = interrupt_wait;
    sigemptyset(&action.sa_mask);
    action.sa_flags = 0; /* Do not restart interrupted blocking calls. */
    if (sigaction(SIGALRM, &action, NULL) != 0) return 1;
    session = fsparse_ipc_start(argv[0], fsparse_ipc_header_bytes(), &err);
    if (session == NULL || err != 0) return 1;
    header = fsparse_ipc_data(session);

    header->opcode = 1;
    header->status = FSPARSE_ST_ERROR;
    if (fsparse_ipc_call(session) != FSPARSE_ST_SINGULAR) failures++;

    /* A 250 ms helper response must survive multiple 100 ms death polls. */
    header->opcode = 2;
    header->status = FSPARSE_ST_ERROR;
    if (fsparse_ipc_call(session) != FSPARSE_ST_SINGULAR) failures++;

    header->opcode = 3;
    header->status = FSPARSE_ST_ERROR;
    if (set_timer(0) != 0) failures++;
    status = fsparse_ipc_call(session);
    alarm(0);
    if (status != FSPARSE_ST_SINGULAR || signals_received == 0) failures++;

    header->opcode = 4;
#if defined(__APPLE__)
    /* Repeated interruptions must not postpone Darwin's absolute deadline.
     * The handler stops after 100 signals so failure still permits teardown. */
    if (set_timer(1) != 0) failures++;
#endif
    status = fsparse_ipc_call(session);
    alarm(0);
    if (status != FSPARSE_ST_ERROR) failures++;
#if defined(__APPLE__)
    if (signals_received == 0 || signals_received >= 100) failures++;
#endif
    fsparse_ipc_stop(session);
    /* Both success and failure paths above must reap the one owned child. */
    if (waitpid(-1, NULL, WNOHANG) != -1 || errno != ECHILD) failures++;
    if (failures) fprintf(stderr, "IPC wait failures: %d\n", failures);
    return failures != 0;
}
