import argparse
from pathlib import Path
import random
import subprocess
import unittest
from unittest.mock import patch, MagicMock

from run_gpu_spd import parse_metrics, phase_objective, sample, sample_block, build_sources, unavailable


class RunnerTests(unittest.TestCase):
    def arguments(self):
        return argparse.Namespace(bin_dir=Path('/test/bin'), matrix=None, permutation=None,
                                  repeats=8, device=1, variant='training')

    def test_parser_preserves_order_hash_and_unavailable_memory(self):
        row = parse_metrics('noise\nbackend,family,ordering_hash,solver_version,device_total_peak_bytes\n'
                            'cudss,grid2,FFFFFFFFFFFFFFFF,0.8.0.10,-1\n')
        self.assertEqual(row['ordering_hash'], 'FFFFFFFFFFFFFFFF')
        self.assertEqual(row['solver_version'], '0.8.0.10')
        self.assertEqual(row['device_total_peak_bytes'], -1)
        self.assertEqual(parse_metrics('backend,device_total_peak_source\ncudss,nvml_accounting\n')
                         ['device_total_peak_source'], 'nvml_accounting')
        with self.assertRaises(ValueError):
            parse_metrics('backend,factor_s\ncudss,nan\n')
        for value, expected in [('T', True), ('F', False), ('.TRUE.', True), ('0', False)]:
            self.assertEqual(parse_metrics('backend,reused\ncholmod,' + value + '\n')['reused'], expected)

    def test_objectives_distinguish_reuse_and_repeated_totals(self):
        row = dict(workload='rhs1', factor_s=3., solve_s=2., refactor_s=.5, repeats=8)
        self.assertEqual(phase_objective(row), 5.)
        self.assertEqual(phase_objective(dict(row, workload='refactor')), 2.5)
        self.assertEqual(phase_objective(dict(row, workload='repeated-rhs')), 2.)
        self.assertIsNone(phase_objective(dict(row, factor_s=-1.)))

    def test_unsupported_resident_and_common_paths_stay_visible(self):
        self.assertIsNotNone(unavailable('native-cuda', 'resident', 'tuned'))
        self.assertIsNotNone(unavailable('cholmod', 'resident', 'tuned'))
        self.assertIsNotNone(unavailable('native-cuda', 'host', 'common'))
        self.assertIsNone(unavailable('cudss', 'resident', 'common'))
        with patch('run_gpu_spd.subprocess.Popen') as process:
            record = sample(self.arguments(), [0], 'native-cuda', 'grid2', 1, 1,
                            'rhs1', 'resident', 'tuned', 0)
        self.assertEqual(record['status'], 'unavailable')
        process.assert_not_called()

    def test_native_aliases_keep_cpu_and_gpu_devices_distinct(self):
        with patch.object(Path, 'is_file', return_value=False):
            for backend, device, ordering in [('native', -1, 'tuned'),
                    ('native-nd', -1, 'nd'), ('native-cuda', 1, 'tuned'),
                    ('native-cuda-nd', 1, 'nd')]:
                record = sample(self.arguments(), [0], backend, 'grid2', 1, 1,
                                'rhs1', 'host', 'tuned', 0)
                self.assertEqual(record['requested_device'], device)
                self.assertEqual(record['command'][9], str(device))
                self.assertEqual(record['ordering'], ordering)
                self.assertEqual(record['requested_backend'], backend)

    def test_cpu_build_override_keeps_gpu_executable_separate(self):
        args = self.arguments()
        args.cpu_bin_dir = Path('/cpu/bin')
        with patch.object(Path, 'is_file', return_value=False):
            cpu = sample(args, [0], 'native', 'grid2', 1, 1, 'warm', 'host', 'tuned', 0)
            gpu = sample(args, [0], 'native-cuda', 'grid2', 1, 1, 'warm', 'host', 'tuned', 0)
        self.assertEqual(cpu['command'][3], '/cpu/bin/bench_spd_native_gpu')
        self.assertEqual(gpu['command'][3], '/test/bin/bench_spd_native_gpu')

    def test_reference_processes_share_randomized_blocks_and_backend_names(self):
        args = self.arguments()
        args.reference_bin_dir = Path('/frozen/bin')
        args.reference_backends = 'native-cuda,native-cuda-nd'
        backends = ['native-cuda', 'native-cuda-nd', 'cudss']
        with patch.object(Path, 'is_file', return_value=False):
            records = sample_block(args, [0], backends, 'grid3', 8, 1, 'warm',
                                   'host', 'tuned', 2, random.Random(7))
        self.assertEqual(len(records), 5)
        self.assertEqual({r['requested_backend'] for r in records}, set(backends))
        self.assertEqual([r['series'] for r in records],
                         ['reference', 'candidate', 'reference', 'candidate', 'candidate'])
        self.assertTrue(all(r['trial'] == 2 and r['threads'] == 1 for r in records))
        for row in records:
            expected = '/frozen/bin' if row['series'] == 'reference' else '/test/bin'
            self.assertEqual(str(Path(row['command'][3]).parent), expected)
        self.assertEqual(args.bin_dir, Path('/test/bin')) # copying reference args is isolated

    def test_reference_cpu_copy_overrides_candidate_cpu_build_only_for_reference(self):
        args = self.arguments()
        args.cpu_bin_dir = Path('/candidate-cpu/bin')
        args.reference_bin_dir = Path('/reference-cpu/bin')
        args.reference_backends = 'native'
        with patch.object(Path, 'is_file', return_value=False):
            records = sample_block(args, [0], ['native'], 'grid2', 1, 1, 'rhs1',
                                   'host', 'tuned', 0, random.Random(0))
        self.assertEqual({(r['series'], r['command'][3]) for r in records}, {
            ('candidate', '/candidate-cpu/bin/bench_spd_native_gpu'),
            ('reference', '/reference-cpu/bin/bench_spd_native_gpu')})
        self.assertTrue(all(r['requested_device'] == -1 for r in records))

    def test_reference_build_metadata_keeps_separate_cache_source(self):
        args = self.arguments()
        args.cpu_bin_dir = Path('/cpu/bin')
        args.reference_bin_dir = Path('/reference/bin')
        with patch('run_gpu_spd.build_provenance', side_effect=lambda path: {'cache_path': str(path)}):
            sources = build_sources(args)
        self.assertEqual(sources['build_provenance']['cache_path'], '/test/bin')
        self.assertEqual(sources['cpu_build_provenance']['cache_path'], '/cpu/bin')
        self.assertEqual(sources['reference_build_provenance']['cache_path'], '/reference/bin')
        self.assertIn('not an attestation', sources['reference_build_provenance']['source_hash_scope'])

    def test_failed_vendor_retains_partial_numeric_output(self):
        process = MagicMock(returncode=2)
        process.communicate.return_value = ('backend,vendor_status,backward_error,factor_s\n'
                                             'cudss,6,-1,0.1\n', 'allocation failed')
        with patch.object(Path, 'is_file', return_value=True), \
                patch('run_gpu_spd.subprocess.Popen', return_value=process):
            record = sample(self.arguments(), [0], 'cudss', 'grid2', 1, 1,
                            'rhs1', 'host', 'tuned', 0)
        self.assertEqual(record['status'], 'failed')
        self.assertEqual(record['factor_s'], .1)
        self.assertEqual(record['stderr'], 'allocation failed')

    def test_timeout_kills_only_own_process_group_and_retains_receipt(self):
        process = MagicMock(pid=123)
        process.communicate.side_effect = [subprocess.TimeoutExpired('driver', 1),
                                            ('partial output', 'partial error')]
        with patch.object(Path, 'is_file', return_value=True), \
                patch('run_gpu_spd.subprocess.Popen', return_value=process), \
                patch('run_gpu_spd.os.killpg') as kill:
            record = sample(self.arguments(), [0], 'cudss', 'grid2', 1, 1,
                            'rhs1', 'host', 'tuned', 0)
        self.assertEqual(record['status'], 'timeout')
        self.assertEqual(record['stdout'], 'partial output')
        kill.assert_called_once()
        self.assertEqual(kill.call_args.args[0], 123)


if __name__ == '__main__':
    unittest.main()
