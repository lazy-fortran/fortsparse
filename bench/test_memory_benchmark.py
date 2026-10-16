"""Exercise the public low-memory CLI and the existing numeric CSV reader."""
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from run_spd import sample

binary = Path(sys.argv[1]).resolve()
args = SimpleNamespace(bin_dir=binary.parent, matrix=None, heldout=False)
cpus = [min(os.sched_getaffinity(0))]
record = sample(args, cpus, "native-low-memory", "band", 1, 1, "rhs8", 0)
assert record["status"] == "ok", record
required = int(record["predicted_host_peak_bytes"])
assert required > 0
assert record["memory_plan_status"] == 0
assert record["memory_prediction_quality"] == 1  # exact array layout
assert record["host_current_bytes"] == required
assert record["host_peak_bytes"] == required
assert record["strict_memory_supported"] == 1
assert record["host_memory_limit_bytes"] == 0
assert record["backward_error"] < 1e-12

# Preserve argument 7 (ordering) and 8 (matrix path), then exercise cap 9.
command = record["command"] + ["0", "", str(required)]
result = subprocess.run(command, text=True, capture_output=True, timeout=1)
assert result.returncode == 0, result.stdout + result.stderr
lines = result.stdout.splitlines()
header = next(i for i, line in enumerate(lines)
              if line.startswith("backend,family,"))
values = next(csv.DictReader(io.StringIO("\n".join(lines[header:]))))
assert int(values["host_memory_limit_bytes"]) == required
assert int(values["host_peak_bytes"]) == required
command[-1] = str(required - 1)
result = subprocess.run(command, text=True, capture_output=True, timeout=1)
assert result.returncode != 0
assert "memory_failure,plan_status," in result.stdout, result.stdout
print("Existing CSV reader, eight RHS, exact cap and one-byte-short rejection pass")
