#ifndef FSBENCH_PROCESS_DEVICE_MEMORY_H
#define FSBENCH_PROCESS_DEVICE_MEMORY_H
#include <stdint.h>
/* Current own-process GPU bytes at one instant, NEVER an exact peak. */
int64_t fsbench_gpu_process_memory(int cuda_device);
#endif
