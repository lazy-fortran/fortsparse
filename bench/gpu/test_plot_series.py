"""Series identity checks without requiring matplotlib or rendering."""
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import MagicMock, patch


class PlotSeriesTests(unittest.TestCase):
    def setUp(self):
        module = importlib.util.spec_from_file_location(
            'plot_series_fixture', Path(__file__).with_name('plot_gpu_spd.py'))
        self.plot = importlib.util.module_from_spec(module)
        replacements = {name: MagicMock() for name in
                        ('matplotlib', 'matplotlib.pyplot', 'matplotlib.lines')}
        with patch.dict(sys.modules, replacements):
            module.loader.exec_module(self.plot)

    def test_scaling_never_connects_distinct_binary_or_reference_series(self):
        rows = []
        for series, digest in [('candidate', 'a'*64), ('reference', 'a'*64), ('candidate', 'b'*64)]:
            for threads in [1, 2]:
                rows.append(dict(requested_backend='native-cuda', ordering='tuned',
                    pool_budget_bytes=0, series=series, executable_sha256=digest, threads=threads))
        groups = self.plot.scaling_groups(rows)
        self.assertEqual(len(groups), 3)
        self.assertTrue(all(sorted(r['threads'] for r in group) == [1, 2]
                            for group in groups.values()))

    def test_cpu_and_gpu_labels_identify_reference_and_binary(self):
        for backend, name in [('native', 'CPU'), ('native-cuda', 'CUDA')]:
            row = dict(requested_backend=backend, series='reference', executable_sha256='f'*64)
            label = self.plot.series_label(row)
            self.assertIn(name, label)
            self.assertIn('reference', label)
            self.assertIn('ffffffff', label)


if __name__ == '__main__':
    unittest.main()
