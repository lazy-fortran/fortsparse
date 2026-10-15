#include <sys/resource.h>
#include <stdio.h>
#include <string.h>

/* Linux ru_maxrss is KiB. Children are reported separately: summing their
 * peaks with the parent's would double-count shared IPC pages. */
void fsbench_process_metrics(double *rss_kib, double *child_rss_kib,
                             double *child_cpu_s) {
    struct rusage own = {0}, child = {0};
    getrusage(RUSAGE_SELF, &own);
    getrusage(RUSAGE_CHILDREN, &child);
#ifdef __APPLE__
    *rss_kib = own.ru_maxrss / 1024.0;
    *child_rss_kib = child.ru_maxrss / 1024.0;
#else
    *rss_kib = own.ru_maxrss;
    *child_rss_kib = child.ru_maxrss;
#ifdef __linux__
    /* ru_maxrss may retain the launcher's pre-exec high-water mark. VmHWM
     * resets on exec and therefore measures this benchmark executable. */
    FILE *status = fopen("/proc/self/status", "r");
    if (status) {
        char line[256];
        while (fgets(line, sizeof line, status)) {
            if (strncmp(line, "VmHWM:", 6) == 0) {
                double kib;
                if (sscanf(line + 6, "%lf", &kib) == 1) *rss_kib = kib;
                break;
            }
        }
        fclose(status);
    }
#endif
#endif
    *child_cpu_s = child.ru_utime.tv_sec + child.ru_stime.tv_sec +
        1.e-6 * (child.ru_utime.tv_usec + child.ru_stime.tv_usec);
}
