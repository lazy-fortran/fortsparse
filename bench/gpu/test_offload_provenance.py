"""Independent miniature Nsight exports; no CUDA/runtime dependency."""
from copy import deepcopy
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from analyze_gpu_spd import analyze, summarize
from offload_provenance import attach, sha256, CONFIG_KEYS, SCHEMA


class TraceProofTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.binary = self.base / "retained-solver"
        self.binary.write_bytes(b"independent fixture executable identity; never executed")
        self.database = self.base / "trace.sqlite"
        with closing(sqlite3.connect(self.database)) as connection:
            connection.executescript('''
                CREATE TABLE PROCESSES(globalPid INTEGER,pid INTEGER,name TEXT);
                CREATE TABLE TARGET_INFO_CUDA_DEVICE(gpuId INTEGER,cudaId INTEGER,pid INTEGER,uuid TEXT);
                CREATE TABLE TARGET_INFO_GPU(id INTEGER,uuid TEXT,busLocation TEXT);
                CREATE TABLE StringIds(id INTEGER,value TEXT);
                CREATE TABLE CUPTI_ACTIVITY_KIND_KERNEL(
                    globalPid INTEGER,deviceId INTEGER,demangledName INTEGER,start INTEGER,end INTEGER);
                INSERT INTO PROCESSES VALUES(99,42,'/capture/solver');
                INSERT INTO PROCESSES VALUES(100,43,'/other/service');
                INSERT INTO TARGET_INFO_CUDA_DEVICE VALUES(1,0,42,'abc');
                INSERT INTO TARGET_INFO_GPU VALUES(1,'abc','0000:08:00.0');
                INSERT INTO StringIds VALUES(7,'independent_numeric_kernel<double>(...)');
                INSERT INTO CUPTI_ACTIVITY_KIND_KERNEL VALUES(99,0,7,1,2);
                INSERT INTO CUPTI_ACTIVITY_KIND_KERNEL VALUES(99,0,7,3,4);
                INSERT INTO CUPTI_ACTIVITY_KIND_KERNEL VALUES(100,0,7,1,2);
            ''')
        self.command = ["taskset", "--cpu-list", "0", "/capture/solver", "grid3", "8", "1",
                        "warm", "training", "0", "host", "tuned", "-", "1", "-"]
        row = dict(requested_backend="cudss", backend="cudss", family="grid3", level=8,
                   workload="warm", variant="training", io_mode="host", threads=1,
                   ordering="tuned", repeats=1, pool_budget_bytes=128 * 2**20,
                   requested_device=0, actual_device=-1, gpu_offload_verified=0,
                   command=self.command, status="ok", primary_phase_s=.2,
                   workflow_cpu_s=.1, parent_rss_kib=100, device_total_peak_bytes=-1)
        environment = {"CUDA_VISIBLE_DEVICES": "1", "FSBENCH_CUDSS_POOL_LIMIT_BYTES": "134217728"}
        self.document = dict(trials=[row], binary_provenance={self.command[3]:
            dict(sha256=sha256(self.binary))}, environment=environment,
            matrix_provenance=None, permutation_provenance=None)
        self.entry = dict(id="one-exact-capture", executable=dict(path=self.binary.name,
            sha256=sha256(self.binary)), capture=dict(command=self.command,
                configuration={key: row[key] for key in CONFIG_KEYS}, environment=environment,
                inputs=dict(matrix_sha256=None, permutation_sha256=None), logical_device=0,
                physical_gpu=dict(index=1, uuid="GPU-ABC", pci_bus_id="0000:08:00.0")),
            trace=dict(path=self.database.name, sha256=sha256(self.database), pid=42),
            kernels=[dict(name="independent_numeric_kernel<double>(...)", count=2,
                          purpose="independent numerical-update fixture")])
        self.manifest = self.base / "proofs.json"
        self.save()

    def save(self):
        self.manifest.write_text(json.dumps(dict(schema=SCHEMA, proofs=[self.entry])))

    def trace_edit(self, statement):
        with closing(sqlite3.connect(self.database)) as connection:
            connection.execute(statement)
            connection.commit()
        self.entry["trace"]["sha256"] = sha256(self.database)
        self.save()

    def test_actual_events_attach_without_overwriting_vendor_metrics(self):
        output = attach(self.document, self.manifest)
        row = output["trials"][0]
        self.assertTrue(row["trace_verified"])
        self.assertEqual(row["trace_actual_device"], 0)
        self.assertEqual(row["gpu_offload_verified"], 0)
        self.assertEqual(row["actual_device"], -1)
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["primary_phase_s"], .2)
        self.assertEqual(row["offload_proof"]["trace_kernel_count"], 2) # other PID excluded
        self.assertNotIn("trace_verified", self.document["trials"][0])
        result = analyze(self.document, required_trials=1, offload_provenance=self.manifest)
        self.assertTrue(result["summaries"][0]["offload_verified"])
        self.assertEqual(result["summaries"][0]["trace_verified_samples"], 1)
        # Re-running analysis cannot trust a previously serialized boolean alone.
        self.assertFalse(analyze(output, required_trials=1)["summaries"][0]["offload_verified"])

    def test_binary_workload_pool_and_input_mismatch_do_not_attach(self):
        for mutation in ["binary", "workload", "pool", "matrix", "environment"]:
            with self.subTest(mutation=mutation):
                document = deepcopy(self.document)
                if mutation == "binary":
                    document["binary_provenance"][self.command[3]]["sha256"] = "0" * 64
                elif mutation == "workload":
                    document["trials"][0]["workload"] = "refactor"
                elif mutation == "pool":
                    document["trials"][0]["pool_budget_bytes"] += 1
                elif mutation == "matrix":
                    document["matrix_provenance"] = dict(sha256="1" * 64)
                else:
                    document["environment"]["CUDA_VISIBLE_DEVICES"] = "0"
                self.assertFalse(attach(document, self.manifest)["trials"][0]["trace_verified"])

    def test_changed_retained_binary_and_trace_fail_closed(self):
        self.binary.write_bytes(b"different executable")
        with self.assertRaisesRegex(ValueError, "executable hash"):
            attach(self.document, self.manifest)
        self.entry["executable"]["sha256"] = sha256(self.binary)
        self.entry["trace"]["sha256"] = "0" * 64
        self.save()
        with self.assertRaisesRegex(ValueError, "trace hash"):
            attach(self.document, self.manifest)

    def test_no_kernels_wrong_device_or_wrong_process_reject_proof(self):
        self.trace_edit("DELETE FROM CUPTI_ACTIVITY_KIND_KERNEL WHERE globalPid=99")
        with self.assertRaisesRegex(ValueError, "kernel counts"):
            attach(self.document, self.manifest)

    def test_other_logical_device_cannot_supply_events(self):
        self.trace_edit("UPDATE CUPTI_ACTIVITY_KIND_KERNEL SET deviceId=1 WHERE globalPid=99")
        with self.assertRaisesRegex(ValueError, "kernel counts"):
            attach(self.document, self.manifest)

    def test_other_executable_cannot_supply_events(self):
        self.trace_edit("UPDATE PROCESSES SET name='/other/executable' WHERE pid=42")
        with self.assertRaisesRegex(ValueError, "process executable"):
            attach(self.document, self.manifest)

    def test_physical_mapping_and_exact_kernel_counts_are_checked(self):
        for mutation in ["uuid", "count", "wildcard"]:
            with self.subTest(mutation=mutation):
                entry = deepcopy(self.entry)
                if mutation == "uuid":
                    entry["capture"]["physical_gpu"]["uuid"] = "GPU-other"
                elif mutation == "count":
                    entry["kernels"][0]["count"] = 3
                else:
                    entry["kernels"][0]["name"] = "*kernel*"
                self.manifest.write_text(json.dumps(dict(schema=SCHEMA, proofs=[entry])))
                with self.assertRaises(ValueError):
                    attach(self.document, self.manifest)

    def test_offline_command_lines_validate_without_executing_binary(self):
        ordinary = self.base / "ordinary.json"
        ordinary.write_text(json.dumps(self.document))
        for script, output in [("offload_provenance.py", "annotated.json"),
                               ("analyze_gpu_spd.py", "audit.json")]:
            command = [sys.executable, str(Path(__file__).with_name(script)), str(ordinary),
                       "--offload-provenance", str(self.manifest),
                       "--output", str(self.base / output)]
            process = subprocess.run(command, capture_output=True, text=True, timeout=1)
            self.assertEqual(process.returncode, 0, process.stderr)
        annotated = json.loads((self.base / "annotated.json").read_text())
        audit = json.loads((self.base / "audit.json").read_text())
        self.assertTrue(annotated["trials"][0]["trace_verified"])
        self.assertTrue(audit["summaries"][0]["offload_verified"])
        self.assertFalse(audit["summaries"][0]["complete"]) # one proof adds no timing trials

    def test_different_binary_series_never_pool(self):
        first = dict(self.document["trials"][0], executable_sha256="a" * 64)
        second = dict(first, executable_sha256="b" * 64)
        rows = summarize([first, second], required_trials=2)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(not row["complete"] for row in rows))


if __name__ == "__main__":
    unittest.main()
