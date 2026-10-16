"""Parse an opt-in helper diagnostic; never infer aggregate IPC memory.

The metric is a single helper's Linux post-exec VmHWM, in KiB. It is a lower
bound on the parent-plus-helper union peak, never a replacement for that peak.
Missing or explicit unavailable receipts retain -1, including older runs.
Malformed/duplicate receipts raise ValueError and cannot supply a memory bound.
"""
import re

ENVIRONMENT = "FORTSPARSE_BENCH_HELPER_METRICS"
PREFIX = "FORTSPARSE_BENCH_HELPER_VMHWM_KIB"
SCOPE = "single_helper_postexec_VmHWM_lower_bound_on_parent_helper_union_peak"


def parse_helper_memory(stderr):
    lines = [line for line in stderr.splitlines() if line.startswith(PREFIX)]
    if len(lines) > 1:
        raise ValueError("duplicate helper post-exec VmHWM receipts")
    result = dict(helper_postexec_hwm_kib=-1, helper_memory_scope=SCOPE,
                  helper_memory_status="unavailable_no_receipt",
                  helper_memory_receipt=lines[0] if lines else None)
    if not lines:
        return result
    if lines[0] == PREFIX + "=unavailable":
        result["helper_memory_status"] = "unavailable"
        return result
    match = re.fullmatch(re.escape(PREFIX) + r"=([1-9][0-9]*)", lines[0])
    if not match or int(match[1]) > 2**63 - 1:
        raise ValueError("invalid helper post-exec VmHWM receipt")
    result.update(helper_postexec_hwm_kib=int(match[1]),
                  helper_memory_status="available")
    return result
