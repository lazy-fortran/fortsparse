import csv
import tempfile
import unittest
from pathlib import Path
from benchmark_analysis import add_objectives, summarize, analyze_admission, export_csv


def observations(backend, times, memory, threads=1, workload="rhs1"):
    records = []
    for trial, seconds in enumerate(times):
        record = dict(backend=backend, family="grid2", level=1, threads=threads,
                      workload=workload, trial=trial, variant="training", status="ok",
                      factor_s=seconds, solve_s=0.0, refactor_s=0.25,
                      parent_cpu_s=seconds, helper_cpu_s=0.0, parent_rss_kib=memory,
                      helper_rss_kib=memory if backend == "umfpack" else 0.0,
                      jvp_s=0.125, vjp_s=0.25)
        add_objectives(record)
        records.append(record)
    return records


class AdmissionTests(unittest.TestCase):
    def analyze(self, *groups, backends=None):
        records = [record for group in groups for record in group]
        backends = backends or list(dict.fromkeys(r["backend"] for r in records))
        return analyze_admission(summarize(records), backends)

    def test_same_comparator_required_for_both_axes(self):
        result = self.analyze(observations("native", [1.5], 15),
                              observations("superlu", [1], 20),
                              observations("mkl-pardiso", [2], 10))
        self.assertEqual(result["overall"]["native"],
                         "not_competitive_over_observed_ranges")

    def test_best_thread_configuration_is_on_frontier(self):
        result = self.analyze(observations("native", [0.9], 5),
                              observations("superlu", [1], 10),
                              observations("superlu", [0.8], 10, threads=4))
        self.assertEqual(result["cases"][0]["known_frontier"],
                         [dict(backend="superlu", threads=4)])
        self.assertEqual(result["overall"]["native"],
                         "not_competitive_over_observed_ranges")

    def test_overlapping_noise_does_not_become_tolerance(self):
        result = self.analyze(observations("native", [0.9, 1.1], 5),
                              observations("superlu", [1.0, 1.2], 10))
        self.assertEqual(result["overall"]["native"], "indeterminate")
        self.assertEqual(result["cases"][0]["status"], "indeterminate_overlap")

    def test_slow_umfpack_does_not_block_known_frontier(self):
        result = self.analyze(observations("native", [0.5], 5),
                              observations("superlu", [1], 10),
                              observations("umfpack", [2], 3))
        self.assertEqual(result["overall"]["native"], "pass_over_observed_ranges")

    def test_fast_umfpack_unknown_memory_can_block_frontier(self):
        result = self.analyze(observations("native", [0.5], 5),
                              observations("superlu", [1], 10),
                              observations("umfpack", [0.8], 3))
        self.assertEqual(result["overall"]["native"], "indeterminate")

    def test_umfpack_parent_lower_bound_can_protect_frontier(self):
        result = self.analyze(observations("native", [0.5], 5),
                              observations("superlu", [1], 10),
                              observations("umfpack", [0.8], 20))
        self.assertEqual(result["overall"]["native"], "pass_over_observed_ranges")

    def test_missing_comparison_and_failed_native_are_explicit(self):
        good = observations("native", [0.5], 5)
        ext = observations("superlu", [1], 10)
        result = self.analyze(good, ext, backends=["native", "superlu", "mkl-pardiso"])
        self.assertEqual(result["overall"]["native"], "indeterminate")
        failed = dict(good[0], status="timeout")
        result = self.analyze(good + [failed], ext)
        self.assertEqual(result["overall"]["native"], "failed_or_unavailable")

    def test_cuda_request_does_not_prove_device_execution(self):
        result = self.analyze(observations("native-cuda", [0.5], 5),
                              observations("superlu", [1], 10))
        self.assertEqual(result["overall"]["native-cuda"],
                         "indeterminate_device_execution")

    def test_refactor_and_derivative_objectives_include_their_work(self):
        refactor = observations("native", [5], 5, workload="refactor")[0]
        derivative = observations("native", [5], 5, workload="derivatives")[0]
        self.assertEqual(refactor["admission_s"], 0.25)
        self.assertEqual(derivative["admission_s"], 5.375)
        self.assertEqual(derivative["derivative_reuse_s"], 0.375)

    def test_profiling_runs_are_excluded_from_admission(self):
        result = self.analyze(observations("native", [0.5], 5,
                                          workload="profile-factor"),
                              observations("superlu", [1], 10,
                                           workload="profile-factor"))
        self.assertEqual(result["cases"], [])
        self.assertEqual(result["overall"]["native"], "indeterminate")

    def test_csv_preserves_failed_trials_and_missing_statistics(self):
        records = observations("native", [1], 5)
        records.append(dict(records[0], status="timeout"))
        result = dict(records=records, summaries=summarize(records))
        with tempfile.TemporaryDirectory() as directory:
            paths = export_csv(result, Path(directory) / "evidence.json")
            with open(paths["trials"]) as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r["status"] for r in rows], ["ok", "timeout"])
            with open(paths["summaries"]) as stream:
                summary = next(csv.DictReader(stream))
            self.assertEqual(summary["factor_mem_mb"], "")


if __name__ == "__main__":
    unittest.main()
