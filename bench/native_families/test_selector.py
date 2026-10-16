#!/usr/bin/env python3
"""Behavioral size-selection checks; no solver or benchmark timings run."""
import argparse
import itertools
import math
import json
import os
import sys
import unittest
from unittest.mock import patch
from run import select_calibration_level, validate_selected_levels, invoke, EXTERNAL, NATIVE, binary_name, positive_i8, memory_environment, record_memory_configuration, parse_csv, driver_arguments, record_umfpack_configuration
from report import summarize, available_backends
from report_ab import aggregates
from report_front import assess, frontier, PEERS, SIX_PEERS, validate_coverage, validate_helper_lower_bounds
from report_private_real import normalize, requested_scaling


class PrivateReport(unittest.TestCase):
    def test_private_configuration_stays_exact_and_distinct(self):
        cap = 2**53+1
        wrapper = dict(backend='native-amd-capped', threads=1, trial=0, status='pass',
            process_receipt=dict(status='ok', returncode=0,
                record=dict(memory_policy=0, host_memory_limit_bytes=cap)),
            independent_audit=dict(status='pass', rhs=[]))
        first, safe = normalize(wrapper, 'fixture')
        self.assertEqual(first['record']['host_memory_limit_bytes'], cap)
        self.assertEqual(safe['host_memory_limit_bytes'], cap)
        wrapper['process_receipt']['record']['host_memory_limit_bytes'] = cap+1
        second, _ = normalize(wrapper, 'fixture')
        with self.assertRaisesRegex(ValueError, 'mixed memory'):
            summarize([first, second])
        point = {k: first[k] for k in ('family', 'backend', 'level', 'workload',
                                       'variant', 'operation')}
        rows = [dict(point, threads=p, memory_policy=0, host_memory_limit_bytes=limit,
                     process_wall_s_count=5, process_wall_s_median=wall)
                for p, limit, wall in [(1, cap, .1), (2, cap+1, .05)]]
        requested_scaling(rows)
        self.assertEqual(rows[1]['requested_budget_wall_speedup_vs_1'], -1)

    def test_requested_budget_speedup_uses_each_backends_own_baseline(self):
        point = dict(family='fixture', level=1, workload='rhs4', variant='training', operation='N')
        rows = [dict(point, backend=b, threads=p, process_wall_s_count=count,
                     process_wall_s_median=wall)
                for b, p, count, wall in [('native-amd', 1, 5, .1),
                    ('native-amd', 4, 5, .05), ('superlu', 1, 5, .2),
                    ('superlu', 4, 5, .05), ('superlu', 8, 0, '')]]
        requested_scaling(rows)
        self.assertEqual(rows[1]['requested_budget_wall_speedup_vs_1'], 2)
        self.assertEqual(rows[1]['requested_budget_wall_efficiency_vs_1'], .5)
        self.assertEqual(rows[3]['requested_budget_wall_speedup_vs_1'], 4)
        self.assertEqual(rows[4]['requested_budget_wall_efficiency_vs_1'], -1)

    def test_wrapper_audit_failure_overrides_successful_child_without_numeric_export(self):
        wrapper = dict(backend='native-amd', threads=1, trial=0, status='failed',
            candidate_path='/private/numeric-export.mtx',
            process_receipt=dict(status='ok', returncode=0, process_wall_s=.1,
                                 record=dict(factor_s=.01, matrix_values=[99, 98])),
            independent_audit=dict(status='failed', rhs=[dict(normwise_backward=.02,
                row_scaled_backward=.03, residual_inf_exact='PRIVATE_FRACTION',
                forward_relative_inf=dict(known_x=.04))]))
        normalized, safe = normalize(wrapper, 'private-fixture')
        self.assertEqual(normalized['status'], 'failed')
        self.assertEqual(safe['child_status'], 'ok')
        self.assertEqual(safe['audit_status'], 'failed')
        self.assertEqual(safe['audit_normwise_backward_max'], .02)
        self.assertNotIn('candidate_path', safe)
        self.assertNotIn('matrix_values', normalized['record'])
        self.assertNotIn('residual_inf_exact', safe)
        wrapper.update(status='pass')
        wrapper['independent_audit']['status'] = 'pass'
        self.assertEqual(normalize(wrapper, 'private-fixture')[0]['status'], 'ok')


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


class MemorySelections(unittest.TestCase):
    def test_capped_alias_requires_exact_positive_i8(self):
        cap=2**53+1
        self.assertEqual(memory_environment('native-amd-capped',cap),
                         {'FORTSPARSE_BENCH_HOST_LIMIT_BYTES':str(cap)})
        for value in (None,0,-1,2**63,'12junk','1.5'):
            with self.assertRaises(argparse.ArgumentTypeError):
                memory_environment('native-amd-capped',value)
        self.assertEqual(memory_environment('superlu',cap),{})
        self.assertEqual(memory_environment('native-amd-low'),{})
        self.assertEqual(memory_environment('native-amd-low',cap),
                         {'FORTSPARSE_BENCH_HOST_LIMIT_BYTES':str(cap)})
        for cap in (2**53+1,2**63-1):
            record=parse_csv('backend,family,level,threads,memory_policy,host_memory_limit_bytes\n'
                             f'native-amd-capped,drift-real,3,1,0,{cap}\n')
            self.assertEqual(record['host_memory_limit_bytes'],cap)
            self.assertIsInstance(record['host_memory_limit_bytes'],int)

    def test_receipt_must_confirm_actual_policy_and_cap(self):
        good=dict(status='ok',record=dict(memory_policy=0,host_memory_limit_bytes=123))
        record_memory_configuration(good,'native-amd-capped',123)
        self.assertEqual(good['status'],'ok')
        self.assertEqual(good['requested_memory_configuration']['host_memory_limit_bytes'],123)
        low=dict(status='ok',record=dict(memory_policy=1,host_memory_limit_bytes=2**63-1))
        record_memory_configuration(low,'native-amd-low',2**63-1)
        self.assertEqual(low['status'],'ok')
        self.assertEqual(low['requested_memory_configuration']['host_memory_limit_bytes'],2**63-1)
        bad=dict(status='ok',record=dict(memory_policy=0,host_memory_limit_bytes=0))
        record_memory_configuration(bad,'native-amd-capped',123)
        self.assertEqual(bad['status'],'configuration_mismatch')
        missing=dict(status='ok',record={})
        record_memory_configuration(missing,'native-amd-low')
        self.assertEqual(missing['status'],'configuration_mismatch')

    def test_plots_include_both_aliases_and_mixed_caps_are_rejected(self):
        self.assertEqual(set(available_backends([{'backend':b} for b in NATIVE])),set(NATIVE))
        with self.assertRaisesRegex(ValueError,'plot backend'):
            available_backends([{'backend':'unknown-native'}])
        point=dict(family='drift-real',backend='native-amd-capped',level=3,threads=1,
                   workload='rhs8',variant='training',operation='N',status='ok')
        attempts=[dict(point,record=dict(memory_policy=0,host_memory_limit_bytes=cap))
                  for cap in (100,200)]
        with self.assertRaisesRegex(ValueError,'mixed memory'):
            summarize(attempts)


class HelperLowerBounds(unittest.TestCase):
    def test_helper_is_a_separate_opt_in_lower_bound_never_an_aggregate(self):
        native=dict(backend='native',threads=1,process_wall_s_min=2,process_wall_s_max=3,
                    parent_rss_kib_min=10,parent_rss_kib_max=11)
        ipc=dict(backend='umfpack',threads=1,process_wall_s_min=.5,process_wall_s_max=1,
                 parent_rss_kib_min=1,parent_rss_kib_max=2,ok=5,
                 helper_postexec_hwm_kib_count=5,helper_postexec_hwm_kib_min=12)
        self.assertEqual(assess(native,[ipc],'process_wall_s')['status'],'indeterminate')
        result=assess(native,[ipc],'process_wall_s',True)
        self.assertEqual(result['status'],'robust_nondominated')
        self.assertEqual(result['excluded_by_helper_lower_bound'],'umfpack:1')
        self.assertEqual(result['excluded_by_memory'],'')
        for changed in (dict(ipc,helper_postexec_hwm_kib_count=4),
                        dict(ipc,helper_postexec_hwm_kib_min=11),
                        dict(ipc,helper_postexec_hwm_kib_count=0)):
            self.assertEqual(assess(native,[changed],'process_wall_s',True)['status'],
                             'indeterminate')
        # Even a smaller helper is never an upper bound proving IPC domination.
        ipc.update(helper_postexec_hwm_kib_min=1)
        self.assertEqual(assess(native,[ipc],'process_wall_s',True)['status'],'indeterminate')

    def test_positive_helper_bounds_require_raw_matching_opt_in_receipts(self):
        row=dict(backend='umfpack',record=dict(helper_postexec_hwm_kib=1234),
                 stderr='FORTSPARSE_BENCH_HELPER_VMHWM_KIB=1234\n',
                 environment=dict(FORTSPARSE_BENCH_HELPER_METRICS='1'))
        validate_helper_lower_bounds([row])
        for changed in (dict(row,stderr=''),dict(row,environment={}),
                        dict(row,stderr=row['stderr']*2)):
            with self.assertRaises(ValueError):validate_helper_lower_bounds([changed])
        validate_helper_lower_bounds([dict(backend='umfpack',record={})])


class UmfpackControls(unittest.TestCase):
    def receipt(self, backend, effective=None):
        controls = dict(schema=1, version=[6, 3, 8], irstep_index=1,
                        default=[.1, 2, .01], effective=effective or [.1, 2, .01])
        return dict(status="ok", record=dict(backend=backend),
                    stdout="umfpack_controls,"+json.dumps(controls)+"\n")

    def test_distinct_commands_and_only_irstep_override(self):
        self.assertEqual(binary_name('umfpack-direct-noir'), binary_name('umfpack-direct'))
        self.assertEqual(driver_arguments('umfpack-direct', ['case']), ['case'])
        self.assertEqual(driver_arguments('umfpack-direct-noir', ['case']),
                         ['case','no-refinement'])
        for backend, effective in [('umfpack-direct',[.1,2,.01]),
                                   ('umfpack-direct-noir',[.1,0,.01])]:
            result=record_umfpack_configuration(self.receipt(backend,effective),backend)
            self.assertEqual(result['status'],'ok')
            self.assertEqual(result['umfpack_controls']['effective'],effective)
        for backend, effective in [('umfpack-direct',[.1,0,.01]),
                                   ('umfpack-direct-noir',[.2,0,.01]),
                                   ('umfpack-direct-noir',[.1,2,.01])]:
            self.assertEqual(record_umfpack_configuration(self.receipt(backend,effective),
                             backend)['status'],'configuration_mismatch')

    def test_missing_duplicate_and_malformed_controls_refused(self):
        for stdout in ('', 'umfpack_controls,broken',
                       self.receipt('umfpack-direct')['stdout']*2):
            row=self.receipt('umfpack-direct');row['stdout']=stdout
            self.assertEqual(record_umfpack_configuration(row,'umfpack-direct')['status'],
                             'configuration_mismatch')

    def test_current_coverage_requires_complete_provenance_and_controls(self):
        selection=dict(backends=list(SIX_PEERS),calibration_scope='complete_declared_baselines')
        provenance=[dict(contract=dict(schema_version=4,required_external_backends=list(SIX_PEERS)))]
        attempts=[]
        for backend, effective in [('umfpack-direct',[.1,2,.01]),
                                   ('umfpack-direct-noir',[.1,0,.01])]:
            attempts.append(dict(self.receipt(backend,effective),backend=backend))
        self.assertEqual(validate_coverage('contract-v4-six',selection,provenance,attempts),SIX_PEERS)
        for bad_selection in (dict(selection,backends=list(PEERS)),
                              dict(selection,calibration_scope='diagnostic_subset')):
            with self.assertRaises(ValueError):
                validate_coverage('contract-v4-six',bad_selection,provenance,attempts)
        attempts[-1]['stdout']=''
        with self.assertRaisesRegex(ValueError,'controls not confirmed'):
            validate_coverage('contract-v4-six',selection,provenance,attempts)

    def test_six_peer_front_requires_the_additional_configuration(self):
        point=dict(family='drift-real',level=3,workload='rhs8',variant='training',
                   operation='N',process_wall_s_min=.1,process_wall_s_max=.2,
                   process_tree_cpu_s_min=.1,process_tree_cpu_s_max=.2,
                   parent_rss_kib_min=100,parent_rss_kib_max=200)
        rows=[dict(point,backend=b,threads=n) for b in PEERS for n in (1,2,4,8,16)]
        rows.append(dict(point,backend='native',threads=1))
        self.assertEqual(len(frontier(rows)),2)
        with self.assertRaisesRegex(ValueError,'incomplete 6-comparator'):
            frontier(rows,SIX_PEERS)
        rows.extend(dict(point,backend='umfpack-direct-noir',threads=n)
                    for n in (1,2,4,8,16))
        self.assertEqual(len(frontier(rows,SIX_PEERS)),2)


class ProcessCPU(unittest.TestCase):
    def test_opt_in_helper_diagnostics_do_not_change_numeric_status(self):
        for suffix, expected in [('=1234',1234),('=0',-1),
                                 ('=1234\nFORTSPARSE_BENCH_HELPER_VMHWM_KIB=2345',-1)]:
            script="import os,sys; assert os.getenv('FORTSPARSE_BENCH_HELPER_METRICS')=='1'; " + \
                "print("+repr('FORTSPARSE_BENCH_HELPER_VMHWM_KIB'+suffix)+",file=sys.stderr); "+ \
                "print('backend,family,level,threads,parent_cpu_s,helper_cpu_s,status'); "+ \
                "print('native,drift-real,3,1,0,0,ok')"
            receipt=invoke(sys.executable,['-c',script],1,[min(os.sched_getaffinity(0))],.8)
            self.assertEqual(receipt['returncode'],0)
            self.assertEqual(receipt['status'],'ok')
            self.assertEqual(receipt['record']['helper_postexec_hwm_kib'],expected)
            self.assertEqual(receipt['helper_postexec_hwm_kib'],expected)
            self.assertEqual(receipt['environment']['FORTSPARSE_BENCH_HELPER_METRICS'],'1')

    def test_unrequested_inherited_cap_is_removed_from_child(self):
        with patch.dict(os.environ, FORTSPARSE_BENCH_HOST_LIMIT_BYTES="17"):
            receipt = invoke(sys.executable, ["-c",
                "import os; print(os.getenv('FORTSPARSE_BENCH_HOST_LIMIT_BYTES','unset'))"],
                1, [min(os.sched_getaffinity(0))], .8)
        self.assertEqual(receipt["stdout"].strip(), "unset")

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
        peers = ["umfpack", "umfpack-direct", "umfpack-direct-noir"]
        row = dict(family="drift-real", level=16, threads=1, trial=0,
                   status="ok", returncode=0, process_wall_s=0.1,
                   workload="rhs8", operation="N", variant="training")
        ipc = dict(row, backend="umfpack")
        direct = dict(row, backend="umfpack-direct")
        self.assertIsNone(select_calibration_level([ipc], "drift-real", [16],
                          peers, [1], 1)[0])
        self.assertIsNone(select_calibration_level([ipc,direct], "drift-real",
                          [16],peers,[1],1)[0])
        noir=dict(row,backend="umfpack-direct-noir")
        self.assertEqual(select_calibration_level([ipc,direct,noir], "drift-real",
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
        for alias in NATIVE:
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
