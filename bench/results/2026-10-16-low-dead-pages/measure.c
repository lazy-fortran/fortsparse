#define _GNU_SOURCE
#include <sched.h>
#include <spawn.h>
extern char **environ;
#include <sys/resource.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <signal.h>

static volatile sig_atomic_t target_pid;
static void terminate_target(int signal_number) {
    if (target_pid > 0) (void)kill((pid_t)target_pid, signal_number);
}
static double seconds(struct timeval t) { return t.tv_sec + 1e-6*t.tv_usec; }
int main(int argc, char **argv) {
    if (argc < 3) return 2;
    int log = open(argv[2], O_CREAT|O_WRONLY|O_TRUNC, 0600);
    if (log < 0) return 3;
    cpu_set_t cpus, actual; CPU_ZERO(&cpus);
    const char *scan = getenv("MEASURE_CPUS");
    if (!scan || !*scan) scan = "0";
    while (*scan) {
        char *end;
        errno = 0;
        long cpu = strtol(scan, &end, 10);
        if (errno || end == scan || cpu < 0 || cpu >= CPU_SETSIZE) return 8;
        CPU_SET((int)cpu, &cpus);
        if (*end && *end != ',') return 8;
        scan = *end ? end + 1 : end;
        if (*end == ',' && !*scan) return 8;
    }
    if (sched_setaffinity(0, sizeof(cpus), &cpus) != 0) return 6;
    if (sched_getaffinity(0, sizeof(actual), &actual) != 0) return 6;
    if (!CPU_EQUAL(&cpus, &actual)) return 6;
    char affinity[CPU_SETSIZE * 5]; size_t offset = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) if (CPU_ISSET(cpu, &actual)) {
        int count = snprintf(affinity + offset, sizeof(affinity) - offset,
                             "%s%d", offset ? "," : "", cpu);
        if (count < 0 || (size_t)count >= sizeof(affinity) - offset) return 8;
        offset += (size_t)count;
    }
    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_adddup2(&actions, log, STDOUT_FILENO);
    posix_spawn_file_actions_adddup2(&actions, log, STDERR_FILENO);
    posix_spawn_file_actions_addclose(&actions, log);
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);
    pid_t child;
    struct sigaction termination = {0};
    termination.sa_handler = terminate_target;
    sigemptyset(&termination.sa_mask);
    if (sigaction(SIGTERM, &termination, NULL) != 0) return 9;
    char **child_argv = calloc((size_t)argc - 1, sizeof(char *));
    if (!child_argv) return 7;
    child_argv[0] = argv[1];
    for (int i = 3; i < argc; ++i) child_argv[i - 2] = argv[i];
    if (posix_spawn(&child, argv[1], &actions, NULL, child_argv, environ) != 0)
        return 4;
    target_pid = child;
    free(child_argv);
    posix_spawn_file_actions_destroy(&actions);
    close(log);
    int status; struct rusage usage;
    while (wait4(child, &status, 0, &usage) < 0) if (errno != EINTR) return 5;
    clock_gettime(CLOCK_MONOTONIC, &end);
    int exit_status = WIFEXITED(status) ? WEXITSTATUS(status) : 128+WTERMSIG(status);
    printf("{\"wall_s\":%.9f,\"user_s\":%.9f,\"system_s\":%.9f,"
           "\"cpu_s\":%.9f,\"maxrss_kib\":%ld,\"exit_status\":%d,"
           "\"cpu_affinity\":\"%s\"}\n",
           end.tv_sec-begin.tv_sec + 1e-9*(end.tv_nsec-begin.tv_nsec),
           seconds(usage.ru_utime), seconds(usage.ru_stime),
           seconds(usage.ru_utime)+seconds(usage.ru_stime), usage.ru_maxrss, exit_status, affinity);
    return exit_status;
}
