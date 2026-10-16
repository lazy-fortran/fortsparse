"""Receipt parser and actual helper exec/HWM behavior; no performance claim."""
import ctypes
import mmap
import os
from pathlib import Path
import struct
import sys
import tempfile
import time

from helper_memory import ENVIRONMENT, PREFIX, parse_helper_memory
from run import MeasuredProcess


def test_parser():
    assert parse_helper_memory("")["helper_postexec_hwm_kib"] == -1
    assert parse_helper_memory(PREFIX + "=unavailable")["helper_postexec_hwm_kib"] == -1
    assert parse_helper_memory("diagnostic\n" + PREFIX + "=1234\n")[
        "helper_postexec_hwm_kib"] == 1234
    for value in ("0", "-1", "nan", "1.5", " 123", "123 trailing", str(2**63)):
        try:
            parse_helper_memory(PREFIX + "=" + value)
        except ValueError:
            pass
        else:
            raise AssertionError(value)
    for text in (PREFIX + "=123\n" + PREFIX + "=123",
                 PREFIX + "=unavailable\n" + PREFIX + "=123", PREFIX + ":123"):
        try:
            parse_helper_memory(text)
        except ValueError:
            pass
        else:
            raise AssertionError(text)


def helper_shutdown(binary, enabled):
    """Exercise the actual helper READY/SHUTDOWN path without factor storage."""
    libc = ctypes.CDLL(None, use_errno=True)
    libc.sem_open.restype = ctypes.c_void_p
    libc.sem_open.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_uint,
                             ctypes.c_uint]
    for name in ("sem_trywait", "sem_post", "sem_close"):
        getattr(libc, name).argtypes = [ctypes.c_void_p]
    libc.sem_unlink.argtypes = [ctypes.c_char_p]
    fd, path = tempfile.mkstemp(prefix="fs-helper-hwm-", dir="/dev/shm")
    request_name = ("/" + Path(path).name + "-req").encode()
    done_name = ("/" + Path(path).name + "-done").encode()
    semaphores, region, process = [], None, None
    try:
        os.ftruncate(fd, 4096)
        region = mmap.mmap(fd, 4096)
        for name in (request_name, done_name):
            handle = libc.sem_open(name, os.O_CREAT | os.O_EXCL, 0o600, 0)
            assert handle not in (None, ctypes.c_void_p(-1).value)
            semaphores.append(handle)
        env = os.environ.copy()
        env.pop(ENVIRONMENT, None)
        env["OPENBLAS_NUM_THREADS"] = "1"
        if enabled:
            env[ENVIRONMENT] = "1"
        # Force fork so the helper inherits the parent's live resident pages
        # before exec. wait4 retains that highwater; /proc VmHWM resets on exec.
        process = MeasuredProcess([str(binary), "/" + Path(path).name,
            request_name.decode(), done_name.decode(), "4096"],
            stdout=-1, stderr=-1, text=True, env=env, preexec_fn=lambda: None,
            start_new_session=True)
        deadline = time.monotonic() + 0.5
        while libc.sem_trywait(semaphores[1]) != 0:
            if time.monotonic() >= deadline or process.poll() is not None:
                raise AssertionError("helper READY handshake failed")
            time.sleep(0.001)
        struct.pack_into("=i", region, 0, 6)  # existing SHUTDOWN opcode
        assert libc.sem_post(semaphores[0]) == 0
        stdout, stderr = process.communicate(timeout=0.5)
        assert process.returncode == 0, stderr
        assert not stdout
        return stderr, process.completed_usage
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.communicate()
        for handle in semaphores:
            libc.sem_close(handle)
        for name in (request_name, done_name):
            libc.sem_unlink(name)
        if region is not None:
            region.close()
        os.close(fd)
        os.unlink(path)


def test_actual_helper(binary):
    if not sys.platform.startswith("linux"):
        return  # The opt-in metric is explicitly unavailable elsewhere.
    memory = bytearray(128 * 1024 * 1024)
    for offset in range(0, len(memory), 4096):
        memory[offset] = 1
    stderr, usage = helper_shutdown(binary, True)
    receipt = parse_helper_memory(stderr)
    assert receipt["helper_memory_status"] == "available", stderr
    assert usage is not None
    assert usage.ru_maxrss >= len(memory) // 1024
    assert 0 < receipt["helper_postexec_hwm_kib"] < len(memory) // 2048
    print("Witness: parent_live_bytes=", len(memory),
          "inherited_wait4_maxrss_kib=", usage.ru_maxrss,
          "helper_postexec_VmHWM_kib=", receipt["helper_postexec_hwm_kib"])
    stderr, _ = helper_shutdown(binary, False)
    assert PREFIX not in stderr  # Product/default behavior emits no metric.
    del memory


if __name__ == "__main__":
    test_parser()
    if len(sys.argv) > 1:
        test_actual_helper(Path(sys.argv[1]).resolve())
    print("Helper receipt validation and requested post-exec HWM checks pass")
