"""Offline, exact-configuration Nsight SQLite evidence; never a timing trial."""
import argparse
from copy import deepcopy
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3

SCHEMA = "fortsparse-gpu-offload-proof-v1"
CONFIG_KEYS = ("requested_backend", "family", "level", "workload", "variant", "io_mode",
               "threads", "ordering", "repeats", "pool_budget_bytes", "requested_device")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def uuid(value):
    return value.lower().removeprefix("gpu-")


def resolved(base, value):
    path = Path(value)
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def validate_entry(entry, base):
    """Hashes/config are capture assertions; kernel/process/device facts are queried."""
    capture, trace = entry["capture"], entry["trace"]
    configuration = capture["configuration"]
    if set(configuration) != set(CONFIG_KEYS) or any(v is None for v in configuration.values()):
        raise ValueError("complete exact configuration required")
    logical = capture["logical_device"]
    if not isinstance(logical, int) or logical < 0 or logical != configuration["requested_device"]:
        raise ValueError("logical/requested device mismatch")
    if "CUDA_VISIBLE_DEVICES" not in capture["environment"]:
        raise ValueError("capture visible-device environment required")
    if set(capture["inputs"]) != {"matrix_sha256", "permutation_sha256"}:
        raise ValueError("matrix/permutation input identities required")
    pool_variable = {"cudss": "FSBENCH_CUDSS_POOL_LIMIT_BYTES",
                     "cholmod": "CHOLMOD_GPU_POOL_BYTES"}.get(configuration["requested_backend"])
    if pool_variable and pool_variable not in capture["environment"]:
        raise ValueError("capture pool environment required")
    command = capture["command"]
    if len(command) < 4 or command[:2] != ["taskset", "--cpu-list"]:
        raise ValueError("capture must retain the full benchmark command")
    binary = resolved(base, entry["executable"]["path"])
    if sha256(binary) != entry["executable"]["sha256"]:
        raise ValueError("retained executable hash mismatch")
    database = resolved(base, trace["path"])
    if sha256(database) != trace["sha256"]:
        raise ValueError("trace hash mismatch")
    kernels = entry["kernels"]
    if not kernels or len({k["name"] for k in kernels}) != len(kernels):
        raise ValueError("distinct documented kernel names required")
    if any(not k.get("purpose") or not isinstance(k["count"], int) or k["count"] <= 0
           for k in kernels):
        raise ValueError("positive kernel counts and documented purposes required")
    try:
        with closing(sqlite3.connect(database.as_uri() + "?mode=ro&immutable=1", uri=True)) as connection:
            processes = list(connection.execute(
                "SELECT globalPid,name FROM PROCESSES WHERE pid=?", (trace["pid"],)))
            if len(processes) != 1 or Path(processes[0][1]) != Path(command[3]):
                raise ValueError("trace process executable mismatch")
            devices = list(connection.execute(
                "SELECT d.uuid,g.id,g.uuid,g.busLocation FROM TARGET_INFO_CUDA_DEVICE d "
                "JOIN TARGET_INFO_GPU g ON d.gpuId=g.id WHERE d.pid=? AND d.cudaId=?",
                (trace["pid"], logical)))
            physical = capture["physical_gpu"]
            if len(devices) != 1:
                raise ValueError("trace device mapping unavailable or ambiguous")
            device_uuid, index, gpu_uuid, bus = devices[0]
            if (index != physical["index"] or uuid(device_uuid) != uuid(physical["uuid"]) or
                    uuid(gpu_uuid) != uuid(physical["uuid"]) or
                    bus.lower() != physical["pci_bus_id"].lower()):
                raise ValueError("trace physical-device mapping mismatch")
            counts = dict(connection.execute(
                "SELECT s.value,COUNT(*) FROM CUPTI_ACTIVITY_KIND_KERNEL k "
                "JOIN StringIds s ON k.demangledName=s.id "
                "WHERE k.globalPid=? AND k.deviceId=? AND k.end>k.start GROUP BY s.value",
                (processes[0][0], logical)))
            if not counts or any(counts.get(k["name"]) != k["count"] for k in kernels):
                raise ValueError("documented GPU kernel counts not present in trace")
    except sqlite3.Error as error:
        raise ValueError("unsupported or invalid Nsight SQLite export: " + str(error)) from error
    return dict(logical_device=logical, physical_gpu=deepcopy(physical),
                documented_kernel_count=sum(k["count"] for k in kernels),
                trace_kernel_count=sum(counts.values()), trace_path=str(database),
                trace_sha256=trace["sha256"], trace_pid=trace["pid"])


def attach(document, manifest_path):
    """Return annotated copies. No timings, statuses, vendor flags or device IDs change."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != SCHEMA:
        raise ValueError("unsupported offload manifest schema")
    checked, identities = [], set()
    for entry in manifest["proofs"]:
        if not entry.get("id") or entry["id"] in identities:
            raise ValueError("unique proof ids required")
        identities.add(entry["id"])
        checked.append((entry, validate_entry(entry, manifest_path.parent)))
    result = deepcopy(document)
    manifest_hash = sha256(manifest_path)
    for row in result["trials"]:
        row["trace_verified"] = False
        row.pop("trace_actual_device", None)
        row.pop("offload_proof", None)
        command = row.get("command", [])
        identity = result.get("binary_provenance", {}).get(command[3] if len(command) > 3 else "", {})
        executable_hash = identity.get("sha256")
        if row.get("executable_sha256", executable_hash) != executable_hash:
            raise ValueError("trial and document executable hashes disagree")
        row["executable_sha256"] = executable_hash
        matches = []
        for entry, facts in checked:
            capture = entry["capture"]
            if (entry["executable"]["sha256"] == executable_hash and
                    capture["command"] == command and
                    capture["inputs"] == {
                        "matrix_sha256": (result.get("matrix_provenance") or {}).get("sha256"),
                        "permutation_sha256": (result.get("permutation_provenance") or {}).get("sha256")} and
                    all(row.get(k) == capture["configuration"][k] for k in CONFIG_KEYS) and
                    all(result.get("environment", {}).get(k) == v
                        for k, v in capture["environment"].items())):
                matches.append((entry, facts))
        if len(matches) > 1:
            raise ValueError("ambiguous proofs for one exact configuration")
        if matches:
            entry, facts = matches[0]
            row.update(trace_verified=True, trace_actual_device=facts["logical_device"],
                       offload_proof=dict(manifest_path=str(manifest_path),
                           manifest_sha256=manifest_hash, proof_id=entry["id"], **facts))
    result["offload_provenance"] = dict(path=str(manifest_path), sha256=manifest_hash,
        validated_proofs=len(checked), matched_trials=sum(r["trace_verified"] for r in result["trials"]),
        limit="SQLite verifies process/device/kernel events. Binary SHA, argv, environment and "
              "configuration are retained capture assertions, matched exactly to ordinary trials. "
              "Trace timings are diagnostic and never pooled with ordinary trials.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--offload-provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = attach(json.loads(args.input.read_text()), args.offload_provenance)
    result["offload_provenance"]["ordinary_input"] = dict(
        path=str(args.input.resolve()), sha256=sha256(args.input))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
