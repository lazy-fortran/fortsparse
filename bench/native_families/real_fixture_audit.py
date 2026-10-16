"""Independent exact discrete residual audit; no NumPy or solver dependency."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import time


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records(path, header):
    with Path(path).open() as stream:
        if stream.readline().strip() != header:
            raise ValueError("unsupported MatrixMarket header")
        for line in stream:
            text = line.strip()
            if text and not text.startswith("%"):
                yield text.split()


def real_array(path):
    lines = records(path, "%%MatrixMarket matrix array real general")
    shape = next(lines)
    if len(shape) != 2:
        raise ValueError("invalid array shape")
    n, columns = map(int, shape)
    if not 0 < n <= 10_000_000 or not 0 < columns <= 16:
        raise ValueError("array shape outside benchmark reader bounds")
    values = []
    for fields in lines:
        if len(fields) != 1:
            raise ValueError("one MatrixMarket array scalar per record required")
        value = float(fields[0])
        if not math.isfinite(value):
            raise ValueError("nonfinite array value")
        values.append(value)
        if len(values) > n * columns:
            raise ValueError("surplus array values")
    if len(values) != n * columns:
        raise ValueError("truncated array")
    return [values[j*n:(j+1)*n] for j in range(columns)]


def real_rows(path):
    lines = records(path, "%%MatrixMarket matrix coordinate real general")
    shape = next(lines)
    if len(shape) != 3:
        raise ValueError("invalid coordinate shape")
    nr, nc, nz = map(int, shape)
    if not 0 < nr == nc <= 10_000_000 or not 0 <= nz <= 100_000_000:
        raise ValueError("matrix outside sparse benchmark reader bounds")
    rows = [[] for _ in range(nr)]
    previous = (0, 0)
    count = 0
    for fields in lines:
        if len(fields) != 3:
            raise ValueError("invalid coordinate record")
        i, j = map(int, fields[:2])
        value = float(fields[2])
        if not 1 <= i <= nr or not 1 <= j <= nc or not math.isfinite(value):
            raise ValueError("invalid index or nonfinite matrix value")
        if (j, i) <= previous:
            raise ValueError("canonical column/row order without duplicates required")
        rows[i-1].append((j-1, Fraction(value)))
        previous = (j, i)
        count += 1
        if count > nz:
            raise ValueError("surplus matrix entries")
    if count != nz:
        raise ValueError("truncated matrix")
    return rows


def relative_inf(x, y, symmetric=False):
    if len(x) != len(y) or not x:
        raise ValueError("forward comparison shape mismatch")
    x, y = list(map(Fraction, x)), list(map(Fraction, y))
    numerator = max(abs(a-b) for a, b in zip(x, y))
    denominator = max(map(abs, y))
    if symmetric:
        denominator = max(denominator, max(map(abs, x)))
    if not denominator:
        return Fraction(0) if not numerator else None
    return numerator / denominator


def audit(rows, rhs, known_x, reference_sets, candidates, policy):
    """All acceptance comparisons use exact rationals, including thresholds."""
    start = time.perf_counter()
    n = len(rows)
    if len(candidates) != 4 or len(rhs) != 4 or len(known_x) != n:
        raise ValueError("original three physical plus one manufactured RHS required")
    if any(len(column) != n for column in candidates + rhs):
        raise ValueError("candidate/RHS dimension mismatch")
    if set(reference_sets) != {"umfpack", "superlu"}:
        raise ValueError("both independent physical reference centers required")
    for values in reference_sets.values():
        if len(values) != 3 or any(len(column) != n for column in values):
            raise ValueError("reference dimension mismatch")
    if not all(math.isfinite(v) for col in candidates for v in col):
        raise ValueError("nonfinite candidate")
    gates = policy["gates"]
    rownorms = [sum((abs(v) for _, v in row), Fraction(0)) for row in rows]
    anorm = max(rownorms)
    metrics, failures = [], []
    for column, (candidate, b) in enumerate(zip(candidates, rhs)):
        x, b = list(map(Fraction, candidate)), list(map(Fraction, b))
        xn, bn = max(map(abs, x)), max(map(abs, b))
        residual = [abs(bi - sum((v*x[j] for j, v in row), Fraction(0)))
                    for row, bi in zip(rows, b)]
        def ratio(numerator, denominator):
            if denominator:
                return numerator / denominator
            return Fraction(0) if not numerator else None
        normwise = ratio(max(residual), anorm*xn + bn)
        scaled = [ratio(r, rowsum*xn + abs(bi))
                  for r, rowsum, bi in zip(residual, rownorms, b)]
        rowwise = max(scaled) if all(v is not None for v in scaled) else None
        for name, value in [("normwise_backward", normwise),
                            ("row_scaled_backward", rowwise)]:
            if value is None or value > Fraction(str(gates[name+"_max"])):
                failures.append(f"rhs{column+1}:{name}")
        forward = ({name: relative_inf(candidate, values[column], True)
                    for name, values in reference_sets.items()} if column < 3
                   else {"known_x": relative_inf(candidate, known_x)})
        key = "physical_relative_inf_max" if column < 3 else "known_x_relative_inf_max"
        for name, value in forward.items():
            if value is None or value > Fraction(str(gates[key])):
                failures.append(f"rhs{column+1}:forward:{name}")
        metrics.append(dict(rhs=column+1, residual_inf=float(max(residual)),
            residual_inf_exact=str(max(residual)),
            normwise_backward=float(normwise) if normwise is not None else None,
            normwise_backward_exact=str(normwise),
            row_scaled_backward=float(rowwise) if rowwise is not None else None,
            row_scaled_backward_exact=str(rowwise),
            forward_relative_inf={key: float(v) if v is not None else None
                                  for key, v in forward.items()}))
    return dict(status="pass" if not failures else "accuracy_failure",
                failures=failures, rhs=metrics,
                audit_wall_s=time.perf_counter()-start,
                arithmetic="exact rational on stored binary64 data and candidates",
                claim="discrete residual audit and empirical forward regression; "
                      "inverse-norm/forward certificate unavailable")


class PrivateFixture:
    def __init__(self, manifest, expected_policy_sha256):
        self.path = Path(manifest).resolve()
        self.base = self.path.parent
        self.manifest = json.loads(self.path.read_text())
        self.manifest_sha256 = sha256(self.path)
        self.policy_path = self.checked_file(self.manifest["accuracy_policy"])
        if sha256(self.policy_path) != expected_policy_sha256:
            raise ValueError("accuracy policy differs from explicitly frozen hash")
        self.policy = json.loads(self.policy_path.read_text())
        self.policy_sha256 = expected_policy_sha256
        self.reader = self.manifest["reader_format"]
        self.matrix = self.checked_file(self.reader["matrix_market"])
        self.physical_rhs = self.checked_file(self.reader["physical_rhs"])
        if sha256(self.matrix) != self.policy["matrix_sha256"] or sha256(
                self.physical_rhs) != self.policy["rhs_sha256"]:
            raise ValueError("matrix or physical RHS differs from frozen policy")
        self.rows = real_rows(self.matrix)
        self.rhs = real_array(self.physical_rhs)
        if len(self.rhs) != 3:
            raise ValueError("three original physical columns required")
        manufactured = real_array(self.checked_file(self.reader["manufactured_rhs"]))
        known = real_array(self.checked_file(self.reader["manufactured_known_x"]))
        if len(manufactured) != 1 or len(known) != 1:
            raise ValueError("one manufactured RHS/known solution required")
        self.rhs += manufactured
        self.known = known[0]
        self.references = {solver: real_array(self.checked_file(
            self.manifest["references"][solver]["physical"]))
            for solver in ("umfpack", "superlu")}

    def checked_file(self, entry):
        path = (self.base / entry["path"]).resolve()
        if not path.is_relative_to(self.base):
            raise ValueError("fixture file escapes private archive")
        if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            raise ValueError("fixture file differs from manifest hash/size")
        return path

    def check(self, candidate):
        return audit(self.rows, self.rhs, self.known, self.references,
                     real_array(candidate), self.policy)
