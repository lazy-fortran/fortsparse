"""Mock NVML legacy ABI and parent/child isolation, without loading a driver."""
import ctypes as C
import struct
import argparse
import subprocess
import unittest
from unittest.mock import MagicMock, patch
from pathlib import Path

from accounting_peak import AccountingPeak, AccountingStats, ULL_MAX
from run_gpu_spd import sample


class Function:
    def __init__(self, callback):
        self.callback = callback

    def __call__(self, *args):
        return self.callback(*args)


class Library:
    def __init__(self, modes=(1,), records=None):
        self.modes, self.records = modes, records or {}
        self.calls, self.shutdowns = [], 0
        self.nvmlInit_v2 = Function(lambda: 0)
        self.nvmlShutdown = Function(self.shutdown)
        self.nvmlDeviceGetCount_v2 = Function(self.count)
        self.nvmlDeviceGetHandleByIndex_v2 = Function(self.handle)
        self.nvmlDeviceGetUUID = Function(self.uuid)
        self.nvmlDeviceGetAccountingMode = Function(self.mode)
        self.nvmlDeviceGetAccountingStats = Function(self.stats)
        self.nvmlDeviceSetAccountingMode = Function(self.forbidden)

    def forbidden(self, *args):
        raise AssertionError('observer must never change device state')

    def shutdown(self):
        self.shutdowns += 1
        return 0

    def count(self, output):
        C.cast(output, C.POINTER(C.c_uint)).contents.value = len(self.modes)
        return 0

    def handle(self, index, output):
        C.cast(output, C.POINTER(C.c_void_p)).contents.value = 100 + index
        return 0

    def mode(self, handle, output):
        C.cast(output, C.POINTER(C.c_uint)).contents.value = self.modes[handle.value-100]
        return 0

    def uuid(self, handle, output, size):
        identifier = ('GPU-' + str(handle.value-100)).encode() + b'\0'
        C.memmove(output, identifier, len(identifier))
        return 0

    def stats(self, handle, pid, output):
        self.calls.append((handle.value-100, pid))
        item = self.records.get((handle.value-100, pid))
        if item is None:
            return 6 # NVML_ERROR_NOT_FOUND; do not enumerate another process
        if isinstance(item, int):
            return item
        memory, running = item
        # Independent native C layout oracle, not constructed through AccountingStats.
        blob = struct.pack('@IIQQQI5I', 13, 17, memory, 31, 123456, running, *([0]*5))
        C.memmove(output, blob, len(blob))
        return 0


class AccountingTests(unittest.TestCase):
    def test_legacy_layout_and_64_bit_maximum(self):
        self.assertEqual(C.sizeof(AccountingStats), 56)
        self.assertEqual(AccountingStats.maxMemoryUsage.offset, 8)
        self.assertEqual(AccountingStats.isRunning.offset, 32)
        library = Library(records={(0, 123): (2**40+17, 0)})
        observer = AccountingPeak(library)
        result = observer.query(123)
        self.assertEqual(result['device_total_peak_bytes'], 2**40+17)
        self.assertEqual(result['device_total_peak_source'], 'nvml_accounting')
        self.assertEqual(result['accounting_gpu_uuid'], 'GPU-0')
        self.assertEqual(result['accounting_start_time_us'], 123456)
        self.assertEqual(result['accounting_is_running'], 0)
        self.assertEqual(library.calls, [(0, 123)])
        observer.close()
        self.assertEqual(library.shutdowns, 1)

    def test_legacy_initialization_symbols_are_supported_without_v2_struct_cast(self):
        library = Library(records={(0, 123): (4096, 0)})
        for base in ['nvmlInit', 'nvmlDeviceGetCount', 'nvmlDeviceGetHandleByIndex']:
            setattr(library, base, getattr(library, base+'_v2'))
            delattr(library, base+'_v2')
        observer = AccountingPeak(library)
        self.assertEqual(observer.query(123)['device_total_peak_bytes'], 4096)
        observer.close()
        library = Library()
        del library.nvmlDeviceGetAccountingStats
        library.nvmlDeviceGetAccountingStats_v2 = Function(library.forbidden)
        observer = AccountingPeak(library)
        self.assertEqual(observer.query(123)['device_total_peak_bytes'], -1)
        self.assertEqual(observer.metadata['status'], 'unavailable')

    def test_disabled_or_missing_nvml_remains_unknown(self):
        library = Library(modes=(0,))
        observer = AccountingPeak(library)
        result = observer.query(123)
        self.assertEqual(result['device_total_peak_bytes'], -1)
        self.assertEqual(result['accounting_reason'], 'no_enabled_accounting_devices')
        self.assertEqual(library.calls, [])
        observer.close()
        with patch('accounting_peak.C.CDLL', side_effect=OSError('unavailable')):
            observer = AccountingPeak()
        self.assertEqual(observer.query(123)['device_total_peak_source'], 'unknown')

    def test_running_unavailable_and_driver_error_are_unknown(self):
        for item, reason in [((4096, 1), 'record_not_terminated'),
                             ((ULL_MAX, 0), 'max_memory_unavailable'),
                             (3, 'accounting_query_error')]:
            observer = AccountingPeak(Library(records={(0, 123): item}))
            result = observer.query(123)
            self.assertEqual(result['device_total_peak_bytes'], -1)
            self.assertEqual(result['accounting_reason'], reason)
            observer.close()

    def test_other_pid_never_supplies_peak_and_multi_gpu_never_sums(self):
        library = Library(records={(0, 999): (2**40, 0)})
        observer = AccountingPeak(library)
        self.assertEqual(observer.query(123)['accounting_reason'], 'no_accounting_record_for_child')
        self.assertEqual(library.calls, [(0, 123)])
        observer.close()
        observer = AccountingPeak(Library(modes=(1, 1), records={
            (0, 123): (4096, 0), (1, 123): (8192, 0)}))
        result = observer.query(123)
        self.assertEqual(result['accounting_reason'], 'multiple_accounting_gpus')
        self.assertEqual(result['device_total_peak_bytes'], -1)
        observer.close()

    def test_runner_queries_own_terminated_pid_and_preserves_vendor_and_cpu_fields(self):
        args = argparse.Namespace(bin_dir=Path('/test/bin'), matrix=None, permutation=None,
                                  repeats=8, device=1, variant='training')
        args.accounting_observer = AccountingPeak(Library(records={(0, 123): (8192, 0)}))
        process = MagicMock(returncode=0, pid=123)
        process.communicate.return_value = ('backend,vendor_status,backward_error,parent_cpu_s,device_total_peak_bytes,pid\n'
                                             'cudss,0,1e-15,0.125,-1,999\n', '')
        with patch.object(Path, 'is_file', return_value=True), \
                patch('run_gpu_spd.subprocess.Popen', return_value=process):
            row = sample(args, [0], 'cudss', 'grid2', 1, 1, 'rhs1', 'host', 'tuned', 0)
        self.assertEqual(row['status'], 'ok')
        self.assertEqual(row['accounting_pid'], 123)
        self.assertEqual(row['device_total_peak_bytes'], 8192)
        self.assertEqual(row['vendor_device_total_peak_bytes'], -1)
        self.assertEqual(row['parent_cpu_s'], .125)
        self.assertEqual(row['workflow_cpu_s'], .125)
        self.assertGreaterEqual(row['accounting_observer_cpu_s'], 0)
        args.accounting_observer.close()


    def test_timeout_failure_and_invalid_output_keep_status_and_query_exited_pid(self):
        for expected in ['timeout', 'failed', 'invalid_output']:
            args = argparse.Namespace(bin_dir=Path('/test/bin'), matrix=None, permutation=None,
                                      repeats=8, device=1, variant='training')
            observer = AccountingPeak(Library(records={(0, 123): (8192, 0)}))
            args.accounting_observer = observer
            process = MagicMock(returncode=2, pid=123)
            if expected == 'timeout':
                process.communicate.side_effect = [subprocess.TimeoutExpired('driver', 1),
                    ('backend,device_total_peak_bytes\ncudss,-1\n', '')]
            elif expected == 'invalid_output':
                process.communicate.return_value = ('backend,factor_s\ncudss,nan\n', '')
            else:
                process.communicate.return_value = ('backend,vendor_status\ncudss,6\n', '')
            with patch.object(Path, 'is_file', return_value=True), \
                    patch('run_gpu_spd.subprocess.Popen', return_value=process), \
                    patch('run_gpu_spd.os.killpg'):
                row = sample(args, [0], 'cudss', 'grid2', 1, 1, 'rhs1', 'host', 'tuned', 0)
            self.assertEqual(row['status'], expected)
            self.assertEqual(row['device_total_peak_bytes'], 8192)
            self.assertEqual(row['accounting_pid'], 123)
            self.assertEqual(process.communicate.call_count, 2 if expected == 'timeout' else 1)
            observer.close()

    def test_unavailable_binary_has_no_child_or_accounting_query(self):
        library = Library(records={(0, 123): (8192, 0)})
        args = argparse.Namespace(bin_dir=Path('/test/bin'), matrix=None, permutation=None,
                                  repeats=8, device=1, variant='training',
                                  accounting_observer=AccountingPeak(library))
        with patch.object(Path, 'is_file', return_value=False), \
                patch('run_gpu_spd.subprocess.Popen') as launch:
            row = sample(args, [0], 'cudss', 'grid2', 1, 1, 'rhs1', 'host', 'tuned', 0)
        self.assertEqual(row['status'], 'unavailable')
        self.assertEqual(library.calls, [])
        launch.assert_not_called()
        args.accounting_observer.close()


if __name__ == '__main__':
    unittest.main()
