import csv
import json
from pathlib import Path
import tempfile
import unittest

from report_spd import build_report, expected_scope, load_runs, write_report
from test_benchmark_analysis import observations


class ReportTests(unittest.TestCase):
    def write_run(self, directory, name, groups, binary="binary-one", **metadata):
        path = Path(directory) / (name + ".json")
        document = dict(records=[r for group in groups for r in group],
                        revision=name, platform="test-platform", cpu_model="test-cpu",
                        physical_cpu_ids=[0, 1, 2, 3],
                        fixed_environment={"OPENBLAS_NUM_THREADS": "1"},
                        binary_provenance={"bench_spd": {"sha256": binary}})
        document.update(metadata)
        path.write_text(json.dumps(document))
        return path

    def test_new_binary_replaces_entire_series_and_keeps_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            old = self.write_run(directory, "old", [observations("native", [1], 5),
                observations("native", [0.5], 5, threads=4)])
            failed = dict(observations("native", [0.25], 5)[0], status="timeout")
            new = self.write_run(directory, "new", [[failed]], binary="binary-two")
            summaries, sources = load_runs([old, new])
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0]["threads"], 1)
        self.assertEqual(summaries[0]["successes"], 0)
        self.assertEqual(summaries[0]["statuses"], {"timeout": 1})
        self.assertEqual(summaries[0]["binary_sha256"], "binary-two")
        self.assertEqual(sources[0]["input_status_counts"], {"ok": 2})
        self.assertEqual(sources[1]["input_status_counts"], {"timeout": 1})

    def test_same_binary_disjoint_threads_merge_but_trials_do_not_pool(self):
        with tempfile.TemporaryDirectory() as directory:
            old = self.write_run(directory, "old", [observations("native", [1, 2], 5)])
            extra = self.write_run(directory, "extra", [observations("native", [0.5], 5, threads=4)])
            new = self.write_run(directory, "new", [observations("native", [3], 5)])
            summaries, _ = load_runs([old, extra, new])
        self.assertEqual([s["threads"] for s in summaries], [1, 4])
        self.assertEqual(summaries[0]["samples"], 1)
        self.assertEqual(summaries[0]["admission_s"]["median"], 3)
        self.assertEqual(summaries[0]["source_revision"], "new")

    def test_changed_environment_prevents_mixed_series(self):
        with tempfile.TemporaryDirectory() as directory:
            old = self.write_run(directory, "old", [observations("native", [1], 5),
                observations("native", [0.5], 5, threads=4)])
            new = self.write_run(directory, "new", [observations("native", [0.25], 5)],
                                 fixed_environment={"OPENBLAS_NUM_THREADS": "4"})
            summaries, _ = load_runs([old, new])
        self.assertEqual(len(summaries), 1)

    def test_report_retains_same_comparator_metrics_actual_team_and_provenance(self):
        native = observations("native-nd", [0.5, 0.6], 5, threads=4)
        for r in native:
            r.update(requested_threads=4, active_threads=2, ordering_used=3,
                     parent_cpu_s=0.8, helper_cpu_s=0.1)
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_run(directory, "evidence", [native,
                observations("superlu", [1, 1.1], 10, threads=4)])
            summaries, sources = load_runs([source])
            scope = expected_scope([("grid2", 1)], ["rhs1"], ["training"])
            report = build_report(summaries, sources, scope, ["native-nd", "superlu"], [4], 2)
            output = Path(directory) / "report.json"
            write_report(report, output)
            with output.with_suffix(".cases.csv").open() as stream:
                row = next(csv.DictReader(stream))
        self.assertEqual(report["scope_admission_status"], "pass_over_observed_ranges")
        self.assertEqual(report["cases"][0]["native"]["active_threads"]["median"], 2)
        self.assertEqual(row["native_backend"], "native-nd")
        self.assertEqual(row["comparator_backend"], "superlu")
        self.assertAlmostEqual(float(row["native_workflow_cpu_s_median"]), 0.9)
        self.assertEqual(row["native_source_revision"], "evidence")
        self.assertEqual(row["native_binary_sha256"], "binary-one")
        self.assertEqual(report["mvp_claim"], "not_evaluated")

    def test_missing_case_or_trial_cannot_become_complete_scope_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_run(directory, "evidence", [observations("native", [0.5], 5),
                observations("superlu", [1], 10)])
            summaries, sources = load_runs([source])
        scope = expected_scope([("grid2", 1), ("grid3", 1)], ["rhs1"], ["training"])
        report = build_report(summaries, sources, scope, ["native", "superlu"], [1], 5)
        self.assertEqual(report["scope_admission_status"], "indeterminate_incomplete_evidence")
        self.assertEqual(report["coverage"]["missing_case_groups"], [scope[1]])
        self.assertEqual(len(report["coverage"]["trial_count_mismatches"]), 2)
        self.assertEqual(report["case_status_counts"]["indeterminate_missing_case"], 1)

    def test_archived_matrix_has_no_invented_heldout_variant(self):
        scope = expected_scope([("bcsstk09", 1)], ["warm"], ["training", "heldout"])
        self.assertEqual(scope, [dict(family="bcsstk09", level=1, workload="warm", variant="training")])

    def test_explicit_scope_retains_extra_rows_without_changing_its_admission(self):
        extra = [dict(observations("superlu", [1], 10)[0], level=7)]
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_run(directory, "evidence", [observations("native", [0.5], 5),
                observations("superlu", [1], 10), extra])
            summaries, sources = load_runs([source])
        scope = expected_scope([("grid2", 1)], ["rhs1"], ["training"])
        report = build_report(summaries, sources, scope, ["native", "superlu"], [1], 1)
        self.assertEqual(report["scope_admission_status"], "pass_over_observed_ranges")
        self.assertEqual(len(report["configurations"]), 3)
        self.assertEqual(report["out_of_scope_case_groups"][0]["level"], 7)

    def test_out_of_scope_fast_native_configuration_cannot_supply_a_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_run(directory, "evidence", [observations("native", [2], 15),
                observations("native", [0.5], 5, threads=4),
                observations("native-nd", [0.5], 5), observations("superlu", [1], 10)])
            summaries, sources = load_runs([source])
        scope = expected_scope([("grid2", 1)], ["rhs1"], ["training"])
        report = build_report(summaries, sources, scope, ["native", "superlu"], [1], 1)
        self.assertEqual(report["scope_admission_status"], "not_competitive_over_observed_ranges")
        self.assertEqual(report["out_of_scope_configurations"], 2)


if __name__ == "__main__":
    unittest.main()
