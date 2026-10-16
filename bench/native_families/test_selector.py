#!/usr/bin/env python3
"""Behavioral size-selection checks; no solver or benchmark timings run."""
import itertools
import math
import unittest
from run import select_calibration_level, validate_selected_levels


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
        with self.assertRaises(ValueError):
            select_calibration_level(rows, "drift-real", self.candidates,
                                     self.peers+["native"], self.threads, self.trials)

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
