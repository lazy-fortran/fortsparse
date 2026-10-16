"""Optional read-only parent NVML accounting observer; creates no CUDA context."""
import ctypes as C
import time

NOT_FOUND = 6
ULL_MAX = (1 << 64) - 1


class AccountingStats(C.Structure):
    # nvmlAccountingStats_t / unversioned nvmlDeviceGetAccountingStats only.
    # nvmlAccountingStats_v2_t has a DIFFERENT layout and is never cast here.
    _fields_ = [("gpuUtilization", C.c_uint), ("memoryUtilization", C.c_uint),
                ("maxMemoryUsage", C.c_ulonglong), ("time", C.c_ulonglong),
                ("startTime", C.c_ulonglong), ("isRunning", C.c_uint),
                ("reserved", C.c_uint * 5)]


class AccountingPeak:
    def __init__(self, library=None):
        wall, cpu = time.monotonic(), time.process_time()
        self.library, self.initialized, self.handles = None, False, []
        self.metadata = dict(
            counter_reference="https://docs.nvidia.com/deploy/nvml-api/latest/api/structnvmlAccountingStats__t.html",
            memory_scope="Driver-reported process-lifetime maximum on one matched GPU, including "
                         "driver/vendor context ownership; not an exact sum of wrapper allocations "
                         "or an aggregate across GPUs.",
            stats_abi="nvmlAccountingStats_t; unversioned stats API",
                             stats_size_bytes=C.sizeof(AccountingStats), devices=[])
        try:
            if C.sizeof(AccountingStats) != 56 or C.sizeof(C.c_uint) != 4:
                raise ValueError("unsupported accounting C ABI")
            self.library = library if library is not None else C.CDLL("libnvidia-ml.so.1")
            init = self.bind("nvmlInit_v2", [], fallback="nvmlInit")
            self.shutdown = self.bind("nvmlShutdown", [])
            get_count = self.bind("nvmlDeviceGetCount_v2", [C.POINTER(C.c_uint)],
                                  fallback="nvmlDeviceGetCount")
            get_handle = self.bind("nvmlDeviceGetHandleByIndex_v2", [C.c_uint, C.POINTER(C.c_void_p)],
                                   fallback="nvmlDeviceGetHandleByIndex")
            get_uuid = self.bind("nvmlDeviceGetUUID", [C.c_void_p, C.POINTER(C.c_char), C.c_uint])
            get_mode = self.bind("nvmlDeviceGetAccountingMode", [C.c_void_p, C.POINTER(C.c_uint)])
            self.get_stats = self.bind("nvmlDeviceGetAccountingStats",
                                      [C.c_void_p, C.c_uint, C.POINTER(AccountingStats)])
            code = init()
            if code:
                raise ValueError(f"nvmlInit failed: {code}")
            self.initialized = True
            count = C.c_uint()
            code = get_count(C.byref(count))
            if code:
                raise ValueError(f"device enumeration failed: {code}")
            for index in range(count.value):
                handle, mode, identifier = C.c_void_p(), C.c_uint(), C.create_string_buffer(80)
                code = get_handle(index, C.byref(handle))
                entry = dict(index=index, status="unavailable", code=code)
                if not code:
                    code = get_mode(handle, C.byref(mode))
                    entry.update(code=code, mode=mode.value if not code else None)
                if not code and mode.value == 1:
                    code = get_uuid(handle, identifier, len(identifier))
                    entry.update(code=code)
                    if not code and identifier.value:
                        identifier = identifier.value.decode("ascii")
                        self.handles.append((handle, identifier))
                        entry.update(status="enabled", uuid=identifier)
                elif not code:
                    entry["status"] = "disabled"
                self.metadata["devices"].append(entry)
            self.metadata["status"] = "ready" if self.handles else "no_enabled_accounting_devices"
        except (OSError, AttributeError, ValueError) as error:
            self.metadata.update(status="unavailable", reason=str(error))
            self.handles = []
        self.metadata.update(initialization_wall_s=time.monotonic()-wall,
                             initialization_cpu_s=time.process_time()-cpu)

    def bind(self, name, arguments, fallback=None):
        try:
            function = getattr(self.library, name)
        except AttributeError:
            if not fallback:
                raise
            function = getattr(self.library, fallback)
        function.argtypes, function.restype = arguments, C.c_int
        return function

    def query(self, pid):
        """Call immediately after the runner's own child has terminated."""
        wall, cpu = time.monotonic(), time.process_time()
        result = dict(device_total_peak_bytes=-1, device_total_peak_source="unknown",
                      accounting_pid=pid, accounting_gpu_uuid=None,
                      accounting_is_running=-1, accounting_start_time_us=-1,
                      accounting_context_time_ms=-1)
        records, errors = [], []
        if not isinstance(pid, int) or pid <= 0 or pid > (1 << 32)-1:
            result["accounting_reason"] = "invalid_child_pid"
        elif not self.handles or not self.initialized:
            result["accounting_reason"] = ("observer_closed" if self.metadata.get("closed") else
                                            self.metadata["status"])
        else:
            for handle, identifier in self.handles:
                stats = AccountingStats()
                code = self.get_stats(handle, pid, C.byref(stats))
                if code == NOT_FOUND:
                    continue
                if code:
                    errors.append(dict(uuid=identifier, code=code))
                else:
                    records.append((identifier, stats))
            if errors:
                result.update(accounting_reason="accounting_query_error", accounting_errors=errors)
            elif len(records) > 1:
                result.update(accounting_reason="multiple_accounting_gpus",
                              accounting_observed_gpu_uuids=[identifier for identifier, _ in records])
            elif not records:
                result["accounting_reason"] = "no_accounting_record_for_child"
            else:
                identifier, stats = records[0]
                result.update(accounting_gpu_uuid=identifier, accounting_is_running=stats.isRunning,
                              accounting_start_time_us=stats.startTime,
                              accounting_context_time_ms=stats.time)
                if stats.isRunning != 0:
                    result["accounting_reason"] = "record_not_terminated"
                elif stats.maxMemoryUsage == ULL_MAX:
                    result["accounting_reason"] = "max_memory_unavailable"
                else:
                    result.update(device_total_peak_bytes=int(stats.maxMemoryUsage),
                                  device_total_peak_source="nvml_accounting", accounting_reason="ok")
        result.update(accounting_observer_wall_s=time.monotonic()-wall,
                      accounting_observer_cpu_s=time.process_time()-cpu)
        return result

    def close(self):
        wall, cpu = time.monotonic(), time.process_time()
        if self.initialized:
            self.metadata["shutdown_code"] = self.shutdown()
            self.initialized = False
        self.metadata["closed"] = True
        self.metadata.update(shutdown_wall_s=time.monotonic()-wall,
                             shutdown_cpu_s=time.process_time()-cpu)
