"""Independent complex discrete audit; Euclidean modulus, enclosing rationals."""
from fractions import Fraction as Q
from math import isfinite, isqrt
import json
from pathlib import Path
import time
from real_fixture_audit import records, sha256, checked_file


def check_deadline(deadline):
    if deadline is not None and time.perf_counter() >= deadline:
        raise TimeoutError("bounded independent complex audit exceeded its budget")


def scalar(fields):
    if len(fields) != 2:
        raise ValueError("complex scalar requires separate real and imaginary values")
    z = complex(*map(float, fields))
    if not isfinite(z.real) or not isfinite(z.imag):
        raise ValueError("nonfinite complex value")
    return z


def exact(z):
    return Q(z.real), Q(z.imag)


def complex_array(path, deadline=None):
    lines = records(path, "%%MatrixMarket matrix array complex general")
    shape = next(lines)
    if len(shape) != 2:
        raise ValueError("invalid complex array shape")
    n, columns = map(int, shape)
    if not 0 < n <= 10_000_000 or not 0 < columns <= 16:
        raise ValueError("complex array outside benchmark bounds")
    values = []
    for fields in lines:
        if len(values) % 128 == 0:
            check_deadline(deadline)
        values.append(scalar(fields))
        if len(values) > n*columns:
            raise ValueError("surplus complex array values")
    if len(values) != n*columns:
        raise ValueError("truncated complex array")
    return [values[j*n:(j+1)*n] for j in range(columns)]


def complex_rows(path):
    lines = records(path, "%%MatrixMarket matrix coordinate complex general")
    shape = next(lines)
    if len(shape) != 3:
        raise ValueError("invalid complex coordinate shape")
    nr, nc, nz = map(int, shape)
    if not 0 < nr == nc <= 10_000_000 or not 0 <= nz <= 100_000_000:
        raise ValueError("complex matrix outside benchmark bounds")
    rows, previous, count = [[] for _ in range(nr)], (0, 0), 0
    for fields in lines:
        if len(fields) != 4:
            raise ValueError("invalid complex coordinate record")
        i, j = map(int, fields[:2])
        if not 1 <= i <= nr or not 1 <= j <= nc or (j, i) <= previous:
            raise ValueError("canonical unique complex CSC indices required")
        rows[i-1].append((j-1, exact(scalar(fields[2:]))))
        previous, count = (j, i), count+1
        if count > nz:
            raise ValueError("surplus complex matrix entries")
    if count != nz:
        raise ValueError("truncated complex matrix")
    return rows


def modulus(z, bits=128):
    """Exact enclosing bounds for sqrt(re²+im²); relative width < 2**(1-bits)."""
    squared = z[0]*z[0]+z[1]*z[1]
    if not squared:
        return Q(0), Q(0)
    exponent = (squared.numerator.bit_length()-squared.denominator.bit_length())//2
    scale = Q(2)**(bits-exponent)
    q = squared*scale*scale
    lower = isqrt(q.numerator//q.denominator)
    upper = lower if lower*lower*q.denominator == q.numerator else lower+1
    return Q(lower)/scale, Q(upper)/scale


def ratio(numerator, denominator):
    if not numerator[1]:
        return Q(0), Q(0)  # Includes the declared zero-over-zero convention.
    if not denominator[1]:
        return None, None
    return (numerator[0]/denominator[1],
            numerator[1]/denominator[0] if denominator[0] else None)


def upper_max(bounds):
    return max(v[0] for v in bounds), max(v[1] for v in bounds)


def audit(rows, rhs, centers, candidates, policy, deadline=None):
    start = time.perf_counter()
    if deadline is None:
        deadline = start+policy['audit_budget_s']
    bits = policy['exact_audit']['sqrt_relative_bits']
    if bits != 128 or len(rhs) != 2 or len(centers) != 2 or len(candidates) != 2:
        raise ValueError("frozen complex physical1/known1 and 128-bit enclosure required")
    n = len(rows)
    if any(len(c) != n for c in rhs+centers+candidates):
        raise ValueError("complex candidate/RHS/center shape mismatch")
    failures, ambiguous, metrics = [], [], []
    try:
        rownorms = []
        terms = 0
        def tick():
            nonlocal terms
            if terms % 64 == 0:
                check_deadline(deadline)
            terms += 1
        def bounded_modulus(value):
            tick()
            return modulus(value,bits)
        def bounded_exact(value):
            tick()
            return exact(value)
        for row in rows:
            lo = hi = Q(0)
            for _, value in row:
                lower, upper = bounded_modulus(value)
                lo, hi = lo+lower, hi+upper
            rownorms.append((lo, hi))
        anorm = upper_max(rownorms)
        for column, (values, bvalues, target) in enumerate(zip(candidates,rhs,centers)):
            x, b, center = [list(map(bounded_exact, v)) for v in (values,bvalues,target)]
            xn = upper_max([bounded_modulus(v) for v in x])
            bn = upper_max([bounded_modulus(v) for v in b])
            residual, scaled = [], []
            for i, row in enumerate(rows):
                rr, ri = b[i]
                for j, (ar, ai) in row:
                    tick()
                    rr, ri = rr-ar*x[j][0]+ai*x[j][1], ri-ar*x[j][1]-ai*x[j][0]
                r = bounded_modulus((rr,ri))
                residual.append(r)
                ab = bounded_modulus(b[i])
                scaled.append(ratio(r, (rownorms[i][0]*xn[0]+ab[0],
                                        rownorms[i][1]*xn[1]+ab[1])))
            normwise = ratio(upper_max(residual),
                             (anorm[0]*xn[0]+bn[0],anorm[1]*xn[1]+bn[1]))
            rowwise = (max((v[0] for v in scaled),default=None),
                       max(v[1] for v in scaled) if all(v[1] is not None for v in scaled) else None)
            distance = upper_max([bounded_modulus((v[0]-c[0],v[1]-c[1]))
                                  for v,c in zip(x,center)])
            forward = ratio(distance, upper_max([bounded_modulus(v) for v in center]))
            report = dict(rhs=column+1, role='physical' if column == 0 else 'known_x')
            for name, bound, key in [('normwise_backward',normwise,'normwise_backward_max'),
                ('row_scaled_backward',rowwise,'row_scaled_backward_max'),
                ('forward_relative_inf',forward,'physical_relative_inf_max' if column == 0
                 else 'known_x_relative_inf_max')]:
                gate = Q(str(policy['gates'][key]))
                if bound[0] is None or bound[0] > gate:
                    failures.append(f'rhs{column+1}:{name}')
                elif bound[1] is None or bound[1] > gate:
                    ambiguous.append(f'rhs{column+1}:{name}')
                report[name+'_lower_exact'] = str(bound[0])
                report[name+'_upper_exact'] = str(bound[1])
            metrics.append(report)
        check_deadline(deadline)
    except TimeoutError as error:
        return dict(status='audit_timeout',reason=str(error),rhs=metrics,
                    audit_wall_s=time.perf_counter()-start)
    return dict(status='accuracy_failure' if failures else 'audit_indeterminate' if ambiguous
                else 'pass',failures=failures,indeterminate=ambiguous,rhs=metrics,
                audit_wall_s=time.perf_counter()-start,
                arithmetic='exact Gaussian rationals; enclosing Euclidean moduli with 128 relative bits',
                claim='discrete residual and empirical center regression; forward enclosure unavailable')


class ComplexFixture:
    driver_case, workload = 'complex-mm', 'rhs2'
    audit_module = 'complex_fixture_audit.py'

    def __init__(self, manifest, expected_policy_sha256, fixture_id=None):
        self.path = Path(manifest).resolve()
        self.base = self.path.parent
        collection = json.loads(self.path.read_text())
        if collection.get('schema') != 'mephit-cylinder-complex-mm-v1' or not fixture_id:
            raise ValueError('explicit complex collection schema and fixture ID required')
        matches = [f for f in collection['fixtures'] if f['id'] == fixture_id]
        if len(matches) != 1:
            raise ValueError('complex fixture ID missing or ambiguous')
        self.selection = matches[0]
        self.manifest = dict(collection,id=fixture_id)
        self.manifest_sha256 = sha256(self.path)
        self.policy_path = checked_file(self.base,collection['accuracy_policy'],require_size=False)
        if sha256(self.policy_path) != expected_policy_sha256:
            raise ValueError('complex policy differs from explicitly frozen hash')
        self.policy_sha256 = expected_policy_sha256
        self.policy = json.loads(self.policy_path.read_text())
        if (self.policy.get('schema') != 'fortsparse-complex-fixture-accuracy-v1' or
                self.policy.get('operation') != 'N' or
                (self.policy.get('nrhs'),self.policy.get('physical_rhs_count'),
                 self.policy.get('known_rhs_count')) != (2,1,1) or
                self.policy['exact_audit']['sqrt_relative_bits'] != 128):
            raise ValueError('unsupported frozen complex audit policy')
        self.ordinary_workflow_budget_s = self.policy['ordinary_child_plus_audit_budget_s']
        f = self.selection
        if (f['operation'] != 'N' or (f['nrhs'],f['physical_rhs_count'],f['known_rhs_count'])
                != (2,1,1)):
            raise ValueError('complex consumer requires declared N physical1/known1')
        self.inputs = {}
        for key, name in [('matrix_file','matrix.mtx'),('rhs_file','rhs.mtx'),
                          ('centers_file','centers.mtx')]:
            filename = f[key]
            entry = dict(f['files'][filename],path=str(Path(f['directory'])/filename))
            self.inputs[name] = self.checked_file(entry)
        self.rows = complex_rows(self.inputs['matrix.mtx'])
        self.rhs = complex_array(self.inputs['rhs.mtx'])
        self.centers = complex_array(self.inputs['centers.mtx'])
        if (len(self.rhs) != 2 or len(self.centers) != 2 or
                any(len(c) != f['n'] for c in self.rhs+self.centers)):
            raise ValueError('complex RHS/center dimensions differ from frozen roles')
        if len(self.rows) != f['n'] or sum(map(len,self.rows)) != f['nnz']:
            raise ValueError('complex matrix dimensions differ from frozen manifest')

    def checked_file(self, entry):
        return checked_file(self.base,entry)

    def input_paths(self):
        return self.inputs

    def check(self, candidate, budget_s=None):
        start = time.perf_counter()
        budget = self.policy['audit_budget_s']
        if budget_s is not None:
            budget = min(budget,max(0,budget_s))
        deadline = start+budget
        try:
            values = complex_array(candidate,deadline)
            result = audit(self.rows,self.rhs,self.centers,values,self.policy,deadline)
        except TimeoutError as error:
            result = dict(status='audit_timeout',reason=str(error))
        result['audit_wall_s'] = time.perf_counter()-start
        return result
