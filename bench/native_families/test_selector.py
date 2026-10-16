#!/usr/bin/env python3
"""Behavioral size-selection checks; no solver or benchmark timings run."""
import itertools
import math
import json
import os
import sys
import unittest
from run import select_calibration_level, validate_selected_levels, invoke, EXTERNAL, binary_name
from report import summarize
from report_ab import aggregates
from report_front import assess


class ObservedFront(unittest.TestCase):
    def point(self, backend, time_lo, time_hi, rss_lo, rss_hi):
        return dict(backend=backend, threads=1, process_wall_s_min=time_lo,
                    process_wall_s_max=time_hi, parent_rss_kib_min=rss_lo,
                    parent_rss_kib_max=rss_hi)

    def test_unknown_ipc_memory_does_not_become_a_known_dominator(self):
        native = self.point('native', 2, 3, 10, 11)
        ipc = self.point('umfpack', .5, 1, 1, 2)
        self.assertEqual(assess(native, [ipc], 'process_wall_s')['status'], 'indeterminate')
        direct = dict(ipc, backend='umfpack-direct')
        self.assertEqual(assess(native, [direct], 'process_wall_s')['status'], 'dominated')

    def test_strict_time_or_safe_parent_lower_bound_excludes_possible_dominators(self):
        native = self.point('native', 1, 2, 10, 11)
        slower = self.point('superlu', 3, 4, 1, 2)
        larger_parent = self.point('umfpack', .5, 1, 12, 13)
        self.assertEqual(assess(native, [slower, larger_parent],
                                'process_wall_s')['status'], 'robust_nondominated')
        boundary = dict(larger_parent, parent_rss_kib_min=11)
        self.assertEqual(assess(native, [boundary], 'process_wall_s')['status'], 'indeterminate')


class RevisionPairs(unittest.TestCase):
    def test_revision_groups_and_trial_pairing_are_preserved(self):
        point = dict(family="drift-real", backend="native-amd", level=128,
                     threads=1, workload="rhs8", variant="training", operation="N")
        attempts = [dict(point, revision_label=label, trial=trial, status="ok",
                         returncode=0, process_wall_s=wall, process_tree_cpu_s=wall,
                         record=dict(process_wall_s=wall, process_tree_cpu_s=wall))
                    for label, trial, wall in [('scalar', 0, 4), ('indexed', 0, 2),
                                               ('indexed', 1, 1), ('scalar', 1, 4)]]
        rows, ratios = aggregates(attempts)
        self.assertEqual(len(rows), 2)
        self.assertEqual([r['attempts'] for r in rows], [2, 2])
        self.assertEqual(ratios[0]['pairs'], 2)
        self.assertEqual(ratios[0]['process_wall_s_scalar_over_indexed_median'], 3)
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            aggregates(attempts[:-1])
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            aggregates(attempts+[attempts[0]])


class ReceiptSummary(unittest.TestCase):
    def test_failed_and_unknown_metrics_do_not_become_successful_timings(self):
        point = dict(family="drift-real", backend="native", level=128, threads=1,
                     workload="rhs8", variant="training", operation="N")
        rows = [dict(point, status="ok", returncode=0,
                     record=dict(process_wall_s=wall, parent_rss_kib=-1))
                for wall in (0.1, 0.3)]
        rows += [dict(point, status="timeout", process_wall_s=0.91),
                 dict(point, status="failed", record=dict(process_wall_s=0.01))]
        result = summarize(rows)[0]
        self.assertEqual((result["attempts"], result["ok"], result["timeouts"],
                          result["failures"]), (4, 2, 1, 1))
        self.assertAlmostEqual(result["process_wall_s_median"], 0.2)
        self.assertEqual(result["process_wall_s_min"], 0.1)
        self.assertEqual(result["process_wall_s_max"], 0.3)
        self.assertEqual(result["parent_rss_kib_count"], 0)
        self.assertEqual(result["timeout_observed_wall_min_s"], 0.91)


class ProcessCPU(unittest.TestCase):
    @unittest.skipUnless(hasattr(os, "wait4"), "POSIX wait4 required")
    def test_completed_cpu_includes_waited_grandchild(self):
        # Independent fork oracle: child burns CPU, parent waits without burning.
        script = """import os, resource, time, json
pid = os.fork()
if pid == 0:
    end = time.process_time() + 0.03
    while time.process_time() < end: pass
    os._exit(0)
os.waitpid(pid, 0)
r = resource.getrusage(resource.RUSAGE_SELF)
print(json.dumps(dict(parent_cpu=r.ru_utime+r.ru_stime)))
"""
        receipt = invoke(sys.executable, ["-c", script], 1,
                         [min(os.sched_getaffinity(0))], 0.8)
        parent = json.loads(receipt["stdout"])["parent_cpu"]
        self.assertGreater(receipt["process_tree_cpu_s"], parent+0.02)
        self.assertEqual(receipt["returncode"], 0)
        self.assertEqual(receipt["process_tree_cpu_scope"],
                         "completed_process_including_reaped_descendants")

    def test_timeout_cpu_is_explicitly_incomplete(self):
        receipt = invoke(sys.executable, ["-c", "import time; time.sleep(2)"],
                         1, [min(os.sched_getaffinity(0))], 0.03)
        self.assertEqual(receipt["status"], "timeout")
        self.assertEqual(receipt["process_tree_cpu_s"], -1)
        self.assertEqual(receipt["process_tree_cpu_scope"], "incomplete_or_unavailable")


class ExternalSelection(unittest.TestCase):
    peers = ["superlu", "mumps"]
    threads = [1, 2]
    trials = 2
    candidates = [16, 32, 64]

    def sample(self, level, wall=0.1):
        return [dict(family="drift-real", level=level, backend=peer, threads=nt,
                     trial=trial, status="ok", returncode=0, process_wall_s=wall,
                     workload="rhs8", operation="N", variant="training")
                for peer, nt, trial in itertools.product(
                    self.peers, self.threads, range(self.trials))]

    def select(self, rows):
        return select_calibration_level(rows, "drift-real", self.candidates,
                                        self.peers, self.threads, self.trials)

    def manifest(self):
        return dict(fixture_source_sha256="frozen-fixture",
                    required_external_backends=self.peers,
                    thread_sweep=self.threads, fresh_process_trials=self.trials,
                    cases={"drift-real": dict(calibration_candidates=self.candidates)})

    def selection(self, peers=None):
        return dict(mode="external_only", fixture_source_sha256="frozen-fixture",
                    levels={"drift-real": 64}, backends=peers or self.peers,
                    threads=self.threads, trials=self.trials)

    def test_direct_and_ipc_umfpack_are_distinct_required_comparators(self):
        peers = ["umfpack", "umfpack-direct"]
        row = dict(family="drift-real", level=16, threads=1, trial=0,
                   status="ok", returncode=0, process_wall_s=0.1,
                   workload="rhs8", operation="N", variant="training")
        ipc = dict(row, backend="umfpack")
        direct = dict(row, backend="umfpack-direct")
        self.assertIsNone(select_calibration_level([ipc], "drift-real", [16],
                          peers, [1], 1)[0])
        self.assertEqual(select_calibration_level([ipc, direct], "drift-real",
                         [16], peers, [1], 1)[0], 16)
        self.assertIn("umfpack-direct", EXTERNAL)
        self.assertNotEqual(binary_name("umfpack"), binary_name("umfpack-direct"))

    def test_old_selection_does_not_gain_the_new_comparator_scope(self):
        rows = sum((self.sample(level) for level in self.candidates), [])
        manifest = self.manifest()
        manifest["required_external_backends"] = self.peers+["umfpack-direct"]
        with self.assertRaisesRegex(ValueError, "omits declared"):
            validate_selected_levels(self.selection(), manifest, ["drift-real"], rows)
        self.assertEqual(validate_selected_levels(self.selection(), manifest,
                         ["drift-real"], rows, diagnostic=True), "diagnostic_subset")

    def test_largest_complete_external_level(self):
        rows = sum((self.sample(level) for level in self.candidates), [])
        selected, outcomes = self.select(rows)
        self.assertEqual(selected, 64)
        self.assertTrue(all(row["complete"] for row in outcomes))
        self.assertEqual(validate_selected_levels(self.selection(), self.manifest(),
                         ["drift-real"], rows), "complete_declared_baselines")

    def test_strict_wall_gate_and_no_extrapolation(self):
        rows = self.sample(16)+self.sample(32, wall=0.5)+self.sample(64, wall=0.01)
        selected, outcomes = self.select(rows)
        self.assertEqual(selected, 16)
        self.assertEqual(outcomes[1]["status"], "wall_gate_exceeded")
        self.assertEqual(outcomes[2]["status"], "unmeasured_after_failure")

    def test_missing_duplicate_and_unavailable_are_not_size_successes(self):
        selected, outcomes = self.select(self.sample(16)[:-1])
        self.assertIsNone(selected)
        self.assertEqual(outcomes[0]["status"], "incomplete_configurations")
        rows = self.sample(16)
        rows[-1] = rows[0].copy()
        self.assertIsNone(self.select(rows)[0])
        rows = self.sample(16)
        rows[0]["status"] = "backend_unavailable"
        self.assertEqual(self.select(rows)[1][0]["status"], "comparator_unavailable")
        rows[0]["status"] = "timeout"
        self.assertEqual(self.select(rows)[1][0]["status"], "timeout")

    def test_bad_returncode_or_nonfinite_wall_rejected(self):
        for wall in (math.inf, math.nan, -0.1):
            self.assertIsNone(self.select(self.sample(16, wall))[0])
        rows = self.sample(16)
        rows[0]["returncode"] = 2
        self.assertIsNone(self.select(rows)[0])

    def test_other_workloads_do_not_prove_the_cold_rhs8_gate(self):
        rows = self.sample(16)
        for row in rows:
            row["workload"] = "rhs1"
        self.assertIsNone(self.select(rows)[0])

    def test_native_observations_cannot_select_size(self):
        rows = sum((self.sample(level) for level in self.candidates), [])
        native = dict(rows[0], backend="native", process_wall_s=100.0, status="timeout")
        self.assertEqual(self.select(rows+[native])[0], 64)
        for alias in ("native", "native-amd"):
            with self.assertRaises(ValueError):
                select_calibration_level(rows, "drift-real", self.candidates,
                                         self.peers+[alias], self.threads, self.trials)

    def test_partial_baselines_require_diagnostic_and_raw_support(self):
        rows = [r for level in self.candidates for r in self.sample(level)
                if r["backend"] == "superlu"]
        partial = self.selection(peers=["superlu"])
        with self.assertRaisesRegex(ValueError, "omits declared"):
            validate_selected_levels(partial, self.manifest(), ["drift-real"], rows)
        self.assertEqual(validate_selected_levels(partial, self.manifest(),
                         ["drift-real"], rows, diagnostic=True), "diagnostic_subset")
        with self.assertRaisesRegex(ValueError, "raw evidence"):
            validate_selected_levels(partial, self.manifest(), ["drift-real"],
                                     rows[:-1], diagnostic=True)
        with self.assertRaisesRegex(ValueError, "unknown selected family"):
            validate_selected_levels(self.selection(), self.manifest(), ["not-frozen"], rows)


if __name__ == "__main__":
    unittest.main()
