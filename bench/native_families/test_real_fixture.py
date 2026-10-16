"""Public tiny data only: independent audit and optional external driver check."""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from real_fixture_audit import audit, real_array, real_rows, sha256
from run import EXTERNAL, binary_name
from run_private_real import verify_qualification
import run_private_real

MATRIX = """%%MatrixMarket matrix coordinate real general
3 3 6
1 1 4
3 1 1
1 2 1
2 2 2
2 3 -1
3 3 3
"""
X = [[1., 0., 2.], [0., 2., 1.], [1., 0., 2.], [-.5, .25, .75]]
B = [[4., -2., 7.], [2., 3., 3.], [4., -2., 7.], [-1.75, -.25, 1.75]]
POLICY = dict(gates=dict(normwise_backward_max=1e-12,
    row_scaled_backward_max=1e-12, physical_relative_inf_max=1e-10,
    known_x_relative_inf_max=1e-10))


def write_array(path, columns):
    path.write_text("%%MatrixMarket matrix array real general\n" +
        f"{len(columns[0])} {len(columns)}\n" +
        "".join(f"{v:.17g}\n" for column in columns for v in column))


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)
        (self.path/"matrix.mtx").write_text(MATRIX)
        self.rows = real_rows(self.path/"matrix.mtx")
        self.refs = dict(umfpack=X[:3], superlu=X[:3])

    def tearDown(self):
        self.temp.cleanup()

    def test_exact_acceptance_and_wrong_answer(self):
        result = audit(self.rows, B, X[3], self.refs, X, POLICY)
        self.assertEqual(result["status"], "pass")
        self.assertTrue(all(r["residual_inf"] == 0 for r in result["rhs"]))
        wrong = [list(column) for column in X]
        wrong[0][0] += .001
        result = audit(self.rows, B, X[3], self.refs, wrong, POLICY)
        self.assertEqual(result["status"], "accuracy_failure")
        self.assertIn("rhs1:forward:umfpack", result["failures"])
        self.assertIn("rhs1:forward:superlu", result["failures"])

    def test_both_references_and_known_vector_are_required(self):
        wrong_refs = dict(self.refs, superlu=[[v+.01 for v in col] for col in X[:3]])
        result = audit(self.rows, B, X[3], wrong_refs, X, POLICY)
        self.assertIn("rhs1:forward:superlu", result["failures"])
        result = audit(self.rows, B, [0., 0., 0.], self.refs, X, POLICY)
        self.assertIn("rhs4:forward:known_x", result["failures"])

    def test_row_scaled_zero_constraint_gate(self):
        rows = [[(0, F(1))], [(1, F(1_000_000))]]
        candidates = [[1e-8, 1.]]*4
        rhs = [[0., 1_000_000.]]*4
        refs = dict(umfpack=[[0., 1.]]*3, superlu=[[0., 1.]]*3)
        result = audit(rows, rhs, [0., 1.], refs, candidates, POLICY)
        self.assertLess(result["rhs"][0]["normwise_backward"], 1e-12)
        self.assertIn("rhs1:row_scaled_backward", result["failures"])

    def test_rounding_residual_is_not_called_zero(self):
        rows = [[(0, F(.1))]]
        candidate = [[.3]]*4
        rhs = [[float(F(.1)*F(.3))]]*4
        refs = dict(umfpack=candidate[:3], superlu=candidate[:3])
        result = audit(rows, rhs, [.3], refs, candidate, POLICY)
        self.assertEqual(result["status"], "pass")
        self.assertGreater(result["rhs"][0]["residual_inf"], 0)

    def test_malformed_and_nonfinite_input(self):
        path = self.path/"bad.mtx"
        for text in [MATRIX.rsplit("3 3 3\n", 1)[0],
                     MATRIX+"1 1 4\n", MATRIX.replace("1 1 4", "0 1 4"),
                     MATRIX.replace("3 1 1", "1 1 1"),
                     MATRIX.replace("1 1 4", "1 1 nan")]:
            path.write_text(text)
            with self.assertRaises(ValueError):
                real_rows(path)
        path.write_text("%%MatrixMarket matrix array real general\n1 1\nnan\n")
        with self.assertRaises(ValueError):
            real_array(path)

    def test_qualification_cannot_be_asserted_or_partial(self):
        class Fixture:
            manifest_sha256, policy_sha256 = "matrix", "policy"
        path = self.path/"qualification.json"
        path.write_text(json.dumps(dict(mode="qualify", status="pass",
            manifest_sha256="matrix", policy_sha256="policy", attempts=[])))
        with self.assertRaisesRegex(ValueError, "all declared"):
            verify_qualification(path, Fixture(), self.path)

    def test_split_qualification_and_changed_binary(self):
        class Fixture:
            manifest_sha256, policy_sha256 = "matrix", "policy"
            base = self.path

            @staticmethod
            def check(candidate):
                return dict(status="pass")
        build = self.path/"build"
        (build/"bin").mkdir(parents=True)
        for name in {binary_name(peer) for peer in EXTERNAL}:
            (build/"bin"/name).write_text("public fake identity; never executed")
        candidate = self.path/"candidate.mtx"
        write_array(candidate, X)
        attempts = [dict(backend=peer, status="pass", candidate_path=str(candidate),
            candidate_sha256=sha256(candidate),
            binary_sha256=sha256(build/"bin"/binary_name(peer))) for peer in EXTERNAL]
        paths = [self.path/"first.json", self.path/"second.json"]
        for path, group in zip(paths, [attempts[:4], attempts[4:]]):
            path.write_text(json.dumps(dict(mode="qualify", status="pass",
                manifest_sha256="matrix", policy_sha256="policy", attempts=group)))
        self.assertEqual(len(verify_qualification(paths, Fixture(), build)), 2)
        (build/"bin"/binary_name("superlu")).write_text("changed")
        with self.assertRaisesRegex(ValueError, "binary changed"):
            verify_qualification(paths, Fixture(), build)


class CampaignTests(unittest.TestCase):
    def test_failed_trial_is_retained_and_later_trials_still_run(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            class Fixture:
                manifest = {"id": "public-test"}
                manifest_sha256, policy_sha256 = "matrix", "policy"
            launched = []
            def child(backend, build, stage, folder, fixture, threads, cores, timeout, host_limit_bytes=None):
                self.assertEqual(host_limit_bytes,2**63-1)
                launched.append((backend, threads, tuple(cores), folder))
                return dict(backend=backend,
                            status="timeout" if len(launched) == 1 else "pass")
            argv = ["run_private_real.py", "--fixture", str(path/"fixture.json"),
                    "--policy-sha256", "policy", "--build", str(path/"build"),
                    "--output", str(path/"output"), "--mode", "measure",
                    "--backend", "native-amd-capped", "--native-host-limit-bytes", str(2**63-1),
                    "--qualification", str(path/"qual.json"),
                    "--thread-sweep", "1,2", "--trials", "2", "--window-note", "unit test"]
            def stage(fixture, output):
                output.mkdir()
                return output
            with mock.patch.object(sys, "argv", argv), \
                 mock.patch.object(run_private_real, "PrivateFixture", return_value=Fixture()), \
                 mock.patch.object(run_private_real, "physical_cores", return_value=[2, 4]), \
                 mock.patch.object(run_private_real, "verify_qualification", return_value={"qual": "hash"}), \
                 mock.patch.object(run_private_real, "prepare", side_effect=stage), \
                 mock.patch.object(run_private_real, "attempt", side_effect=child):
                self.assertEqual(run_private_real.main(), 2)
            receipt = json.loads((path/"output/receipt.json").read_text())
            raw = [json.loads(s) for s in (path/"output/attempts.jsonl").read_text().splitlines()]
            self.assertEqual(receipt["status"], "failed")
            self.assertEqual(receipt["requested_native_host_limit_bytes"],2**63-1)
            self.assertEqual(len(raw), 4*(len(EXTERNAL)+1))
            self.assertEqual(sum(a["status"] == "timeout" for a in raw), 1)
            self.assertEqual(receipt["attempts"], raw)
            self.assertEqual(len({entry[3] for entry in launched}), 4*(len(EXTERNAL)+1))
            for nt in (1, 2):
                for trial in (0, 1):
                    self.assertEqual({a["backend"] for a in raw
                                      if a["threads"] == nt and a["trial"] == trial},
                                     set(EXTERNAL)|{"native-amd-capped"})
            self.assertTrue(all(c == (2,) if nt == 1 else c == (2, 4)
                                for _, nt, c, _ in launched))


def driver_check(binary, backend):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        (path/"matrix.mtx").write_text(MATRIX)
        for name, values in [("physical-rhs.mtx", B[:3]), ("known-rhs.mtx", B[3:]),
                             ("known-x.mtx", X[3:]), ("physical-reference.mtx", X[:3])]:
            write_array(path/name, values)
        target = path/"candidate.mtx"
        args = ["real-mm", "1", "1", "rhs4", "N", "training", str(path)]
        from run import driver_arguments
        args = driver_arguments(backend, args)
        env = dict(os.environ, FORTSPARSE_BENCH_SOLUTION=str(target))
        result = subprocess.run([binary, *args], env=env, capture_output=True,
                                text=True, timeout=.8)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        audit_result = audit(real_rows(path/"matrix.mtx"), B, X[3],
            dict(umfpack=X[:3], superlu=X[:3]), real_array(target), POLICY)
        if audit_result["status"] != "pass":
            raise AssertionError(audit_result)
        print("PASS external real-mm driver, original RHS and independent exact audit")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        driver_check(sys.argv[1], sys.argv[2])
    else:
        unittest.main()
