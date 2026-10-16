"""Standalone comparator ABI/records/link boundary and postexec RSS oracle."""
import os
from pathlib import Path
import re
import subprocess
import sys

from run import parse_csv, record_umfpack_configuration, driver_arguments

binary, mini5 = (Path(arg).resolve() for arg in sys.argv[1:3])
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
ldd = subprocess.run(["ldd", str(binary)], text=True, capture_output=True, check=True)
assert "libumfpack" in ldd.stdout, ldd.stdout
assert not re.search(r"lib(?:fortsparse|[dz]?mumps|mkl)", ldd.stdout, re.I), ldd.stdout
symbols = subprocess.run(["nm", str(binary)], text=True, capture_output=True, check=True)
assert not re.search(r"__fortsparse_(?:solver|backend)_MOD", symbols.stdout), symbols.stdout
assert "umfpack_di_numeric" in symbols.stdout and "umfpack_zi_numeric" in symbols.stdout


def run(workload, backend="umfpack-direct", **kwargs):
    arguments = driver_arguments(backend, ["drift-complex", "3", "1", workload,
        "T", "heldout", str(mini5)])
    result = subprocess.run([str(binary), *arguments], capture_output=True, text=True,
        env=env, timeout=0.8, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    row = parse_csv(result.stdout)
    assert row is not None and row["status"] == "ok", result.stdout
    assert row["backend"] == backend
    receipt = dict(status="ok", record=row, stdout=result.stdout)
    record_umfpack_configuration(receipt, backend)
    assert receipt["status"] == "ok", receipt
    controls = receipt["umfpack_controls"]
    slot = controls["irstep_index"]
    assert controls["default"][slot] == 2, controls
    assert controls["effective"][slot] == (0 if backend.endswith("-noir") else 2), controls
    assert row["backward_error"] < 1e-12 and row["forward_error"] < 1e-10
    assert row["helper_cpu_s"] == 0 and row["helper_rss_kib"] == 0, row
    assert row["parent_rss_kib"] > 0 and row["factor_mem_mb"] > 0, row
    assert row["factor_nnz"] > 0
    assert row["host_current_bytes"] == -1 and row["host_peak_bytes"] == -1
    assert row["parent_cpu_s"] >= row["parent_main_cpu_s"], row
    return row


for backend in ("umfpack-direct", "umfpack-direct-noir"):
    assert run("refactor", backend)["reused"] == "T"
    assert run("warm", backend)["prime_factor_s"] >= 0
invalid = subprocess.run([str(binary), "drift-real", "3", "1", "rhs1", "N",
    "heldout", str(mini5), "invalid-mode"], capture_output=True, text=True, env=env,
    timeout=.8)
assert invalid.returncode != 0 and parse_csv(invalid.stdout) is None
if sys.platform.startswith("linux"):
    # Force fork+exec from a live high-RSS parent. getrusage's inherited HWM
    # would exceed 96 MiB; /proc/self/status VmHWM measures the new image.
    inherited_pages = bytearray(96 * 1024 * 1024)
    row = run("rhs8", preexec_fn=lambda: None)
    assert len(inherited_pages) == 96 * 1024 * 1024
    assert row["parent_rss_kib"] < 64 * 1024, row
print("UMFPACK standalone ABI, records, GPL link isolation and postexec RSS pass")
