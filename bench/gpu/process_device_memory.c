#include "process_device_memory.h"
#ifdef FSBENCH_HAVE_NVML
#include <cuda_runtime_api.h>
#include <nvml.h>
#include <dlfcn.h>
#include <unistd.h>

int64_t fsbench_gpu_process_memory(int device) {
    if (device < 0) return -1;
    void *library = dlopen("libnvidia-ml.so.1", RTLD_NOW | RTLD_LOCAL);
    if (!library) return -1;
    nvmlReturn_t (*init)(void) = (nvmlReturn_t (*)(void))dlsym(library, "nvmlInit_v2");
    nvmlReturn_t (*finish)(void) = (nvmlReturn_t (*)(void))dlsym(library, "nvmlShutdown");
    nvmlReturn_t (*by_pci)(const char *, nvmlDevice_t *) =
        (nvmlReturn_t (*)(const char *, nvmlDevice_t *))dlsym(library, "nvmlDeviceGetHandleByPciBusId_v2");
    nvmlReturn_t (*processes)(nvmlDevice_t, unsigned *, nvmlProcessInfo_t *) =
        (nvmlReturn_t (*)(nvmlDevice_t, unsigned *, nvmlProcessInfo_t *))
        dlsym(library, "nvmlDeviceGetComputeRunningProcesses_v3");
    int64_t bytes = -1;
    if (init && finish && by_pci && processes && init() == NVML_SUCCESS) {
        char pci[32];
        nvmlDevice_t handle;
        if (cudaDeviceGetPCIBusId(pci, sizeof(pci), device) == cudaSuccess &&
            by_pci(pci, &handle) == NVML_SUCCESS) {
            unsigned count = 32;
            nvmlProcessInfo_t records[32];
            if (processes(handle, &count, records) == NVML_SUCCESS)
                for (unsigned i = 0; i < count; ++i)
                    if (records[i].pid == (unsigned)getpid() &&
                        records[i].usedGpuMemory != NVML_VALUE_NOT_AVAILABLE)
                        bytes = (int64_t)records[i].usedGpuMemory;
        }
        finish();
    }
    dlclose(library);
    return bytes;
}
#else
int64_t fsbench_gpu_process_memory(int device) { (void)device; return -1; }
#endif
