import copy
import unittest

from analyze_gpu_spd import OBJECTIVES, analyze, audit_candidate, summarize


def trials(backend, wall, rss, vram):
    return [dict(backend=backend, requested_backend=backend, family="grid3", level=7,
                 workload="warm", variant="training", io_mode="host", requested_device=0,
                 threads=1, ordering="tuned", repeats=1, status="ok", actual_device=0,
                 gpu_offload_verified=1, primary_phase_s=t, parent_rss_kib=rss,
                 device_total_peak_bytes=vram, workflow_cpu_s=t)
            for t in wall]


class FrontAudit(unittest.TestCase):
    def setUp(self):
        self.native = summarize(trials("native-cuda", [1, 2, 2, 2, 3], 20, 30))[0]
        self.external = summarize(trials("cudss", [4, 5, 5, 5, 6], 10, 20))[0]

    def audit(self, objectives="wall_host_device"):
        return audit_candidate(self.native, [self.external], ["cudss"], OBJECTIVES[objectives])

    def test_tradeoff_is_front_member_without_same_point_win(self):
        result = self.audit()
        self.assertEqual(result["status"], "observed_range_nondominated")
        self.assertEqual(result["same_comparator_match_or_beat"], [])

    def test_unknown_native_vram_blocks_device_claim_but_not_host_claim(self):
        self.native["metrics"]["device_total_peak_bytes"] = None
        self.assertEqual(self.audit()["status"], "unknown_native_objectives")
        self.assertEqual(self.audit("wall_host_rss")["status"], "observed_range_nondominated")

    def test_unknown_external_memory_can_be_excluded_by_time(self):
        self.external["metrics"]["device_total_peak_bytes"] = None
        self.assertEqual(self.audit()["status"], "observed_range_nondominated")
        self.native["metrics"]["primary_phase_s"]["maximum"] = 5
        self.assertEqual(self.audit()["status"], "possible_dominator")

    def test_overlap_does_not_become_a_win(self):
        self.native["metrics"]["primary_phase_s"]["maximum"] = 4
        self.assertEqual(self.audit()["status"], "possible_dominator")

    def test_all_external_configurations_must_be_excluded(self):
        other = copy.deepcopy(self.external)
        other["threads"] = 2
        other["metrics"]["primary_phase_s"]["minimum"] = 0.5
        result = audit_candidate(self.native, [self.external, other], ["cudss"],
                                 OBJECTIVES["wall_host_device"])
        self.assertEqual(result["status"], "possible_dominator")
        self.assertEqual(result["possible_dominators"][0]["threads"], 2)

    def test_offload_failure_and_missing_comparator_block(self):
        result = audit_candidate(self.native, [], [], OBJECTIVES["wall_host_rss"])
        self.assertEqual(result["status"], "no_external_comparators")
        self.external["offload_verified"] = False
        self.assertEqual(self.audit()["status"], "external_offload_unverified")
        self.external["offload_verified"] = True
        self.external["complete"] = False
        self.assertEqual(self.audit()["status"], "external_incomplete")
        result = audit_candidate(self.native, [], ["cholmod"], OBJECTIVES["wall_host_rss"])
        self.assertEqual(result["status"], "missing_comparators")

    def test_shared_device_and_short_trial_sets_are_unadmitted(self):
        raw = trials("native-cuda", [1, 2], 20, 30) + trials("cudss", [4, 5], 10, 20)
        result = analyze(dict(trials=raw, measurement_class="diagnostic-shared"))
        self.assertFalse(result["admission_eligible"])
        self.assertEqual(result["cases"][0]["candidates"][0]["objectives"]
                         ["wall_host_rss"]["status"], "native_incomplete")

    def test_pool_configurations_never_pool_trials(self):
        small = trials("cudss", [1, 2, 2, 2, 3], 20, 30)
        large = copy.deepcopy(small)
        for row in small:
            row["pool_limit_bytes"] = 64 * 2**20
        for row in large:
            row["pool_limit_bytes"] = 256 * 2**20
        groups = summarize(small + large)
        self.assertEqual(len(groups), 2)
        self.assertEqual([r["samples"] for r in groups], [5, 5])


if __name__ == "__main__":
    unittest.main()
