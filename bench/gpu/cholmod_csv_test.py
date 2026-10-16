#!/usr/bin/env python3
"""Exercise the real CHOLMOD executable through the shared benchmark runner."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--workload", choices=("refactor", "warm"), default="refactor")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("shared_gpu_runner", args.runner)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    settings = SimpleNamespace(bin_dir=args.driver.resolve().parent, device=0,
        repeats=3, matrix=None, permutation=None, variant="training")
    cpus = sorted(os.sched_getaffinity(0))
    record = runner.sample(settings, cpus, "cholmod", "grid2", 1, 1,
                           args.workload, "host", "tuned", 0)
    assert record["status"] == "ok", json.dumps(record, indent=2)
    assert record["vendor_status"] == 0
    assert record["backward_error"] < 1.e-12
    assert record["forward_error"] <= 1.e-8
    assert record["primary_phase_s"] >= 0
    assert record["workflow_cpu_s"] >= 0
    assert record["device_total_peak_bytes"] == -1
    assert record["helper_rss_kib"] >= 0
    assert record["helper_cpu_s"] >= 0
    assert record["reused"] == (args.workload == "refactor")
    if args.workload == "refactor":
        assert record["refactor_s"] >= 0
    else:
        assert record["warmup_cleanup_s"] >= 0
    print(json.dumps({key: record[key] for key in (
        "status", "vendor_status", "workload", "reused", "backward_error",
        "device_total_peak_bytes", "device_process_sample_count")}, sort_keys=True))


if __name__ == "__main__":
    main()
