"""Public tiny complex data; exact audit and independent dense N/T/H oracles."""
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from complex_fixture_audit import (audit, complex_array, complex_rows, exact,
                                   modulus, ComplexFixture)
from real_fixture_audit import sha256
from run import driver_arguments, parse_csv, record_umfpack_configuration
from run_private_real import verify_qualification

POLICY = dict(schema='fortsparse-complex-fixture-accuracy-v1',operation='N',
    nrhs=2,physical_rhs_count=1,known_rhs_count=1,gates=dict(normwise_backward_max=1e-12,row_scaled_backward_max=1e-12,
    physical_relative_inf_max=1e-10,known_x_relative_inf_max=1e-10),
    audit_budget_s=.5, ordinary_child_plus_audit_budget_s=.9,
    exact_audit=dict(sqrt_relative_bits=128))
A = [[4+1j,1-2j,.5j],[-.5+.25j,5-1j,1+.25j],[.25-1j,-1+2j,6+.5j]]
X = [[1+.25j,-.5+.125j,.25-.5j],[-.25+.5j,.75-.25j,-1+.125j]]


def write_array(path, columns):
    path.write_text('%%MatrixMarket matrix array complex general\n'+
        f'{len(columns[0])} {len(columns)}\n'+''.join(
            f'{v.real:.17e} {v.imag:.17e}\n' for col in columns for v in col))


def fixture(path, operation='N'):
    (path/'matrix.mtx').write_text('%%MatrixMarket matrix coordinate complex general\n3 3 9\n'+
        ''.join(f'{i+1} {j+1} {A[i][j].real:.17e} {A[i][j].imag:.17e}\n'
                for j in range(3) for i in range(3)))
    matrix = A if operation=='N' else [[A[j][i] if operation=='T' else A[j][i].conjugate()
                                      for j in range(3)] for i in range(3)]
    rhs = [[sum(v*x for v,x in zip(row,col)) for row in matrix] for col in X]
    write_array(path/'rhs.mtx',rhs)
    write_array(path/'centers.mtx',X)
    return rhs


class ComplexAuditTests(unittest.TestCase):
    def test_euclidean_enclosures_and_exact_complex_products(self):
        self.assertEqual(modulus(exact(3+4j)),(Q(5),Q(5)))
        lo, hi = modulus(exact(1+1j))
        self.assertLessEqual(lo*lo,Q(2));self.assertGreaterEqual(hi*hi,Q(2))
        self.assertLess(hi-lo,lo/Q(2)**126)
        self.assertLess(hi,Q(2))  # Would fail an undeclared real/imaginary L1 norm.
        x = [[1+2j],[-.5+.25j]]
        rhs = [[(3+4j)*c[0]] for c in x]
        result=audit([[(0,exact(3+4j))]],rhs,x,x,POLICY)
        self.assertEqual(result['status'],'pass')
        self.assertEqual(result['rhs'][0]['normwise_backward_upper_exact'],'0')
        self.assertEqual(audit([[(0,exact(3+4j))]],rhs,x,[[1+2.01j],x[1]],POLICY)
                         ['status'],'accuracy_failure')

    def test_overlapping_enclosure_is_not_accepted(self):
        lower, _ = modulus(exact(1+1j))
        gates=dict(POLICY['gates'],normwise_backward_max=10,row_scaled_backward_max=10,
                   physical_relative_inf_max=str(lower))
        policy=dict(POLICY,gates=gates)
        result=audit([[(0,exact(1+0j))]],[[1+0j]]*2,[[1+0j]]*2,
                     [[2+1j],[1+0j]],policy)
        self.assertEqual(result['status'],'audit_indeterminate')

    def test_row_scaling_zero_denominator_and_audit_deadline(self):
        rows=[[(0,exact(1000+0j))],[(1,exact(.001+0j))]]
        x=[[1+0j,1+0j]]*2;rhs=[[1000+0j,.001+0j]]*2
        candidate=[[1+0j,1.000001+0j]]*2
        result=audit(rows,rhs,x,candidate,POLICY)
        self.assertIn('rhs1:row_scaled_backward',result['failures'])
        self.assertLessEqual(Q(result['rhs'][0]['normwise_backward_upper_exact']),Q('1e-12'))
        zero=[[0j]]*2
        result=audit([[]],zero,zero,[[1e-30+0j],[0j]],POLICY)
        self.assertEqual(result['status'],'accuracy_failure')
        self.assertEqual(audit(rows,rhs,x,x,POLICY,deadline=time.perf_counter()-1)
                         ['status'],'audit_timeout')

    def test_public_manifest_identity_hashes_and_reader_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory); fixture(path)
            policy=path/'accuracy-policy.json';policy.write_text(json.dumps(POLICY))
            names=('matrix.mtx','rhs.mtx','centers.mtx')
            files={name:dict(sha256=sha256(path/name),bytes=(path/name).stat().st_size) for name in names}
            member=dict(id='public-complex-v1',directory='.',n=3,nnz=9,nrhs=2,
                physical_rhs_count=1,known_rhs_count=1,operation='N',matrix_file=names[0],
                rhs_file=names[1],centers_file=names[2],files=files)
            manifest=path/'manifest.json';manifest.write_text(json.dumps(dict(
                schema='mephit-cylinder-complex-mm-v1',id='public-only',fixtures=[member],
                accuracy_policy=dict(path=policy.name,sha256=sha256(policy)))))
            instance=ComplexFixture(manifest,sha256(policy),'public-complex-v1')
            self.assertEqual(instance.check(path/'centers.mtx')['status'],'pass')
            self.assertEqual(instance.check(path/'centers.mtx',budget_s=0)['status'],'audit_timeout')
            with self.assertRaises(ValueError):ComplexFixture(manifest,sha256(policy),'other-id')
            with self.assertRaises(ValueError):ComplexFixture(manifest,'wrong','public-complex-v1')
            qualification=path/'qual.json';qualification.write_text(json.dumps(dict(
                mode='qualify',status='pass',manifest_sha256=instance.manifest_sha256,
                policy_sha256=instance.policy_sha256,fixture_id='another-system',attempts=[])))
            with self.assertRaises(ValueError):verify_qualification(qualification,instance,path)
            original=(path/'matrix.mtx').read_text()
            for text in [original.replace('1 1 ', '0 1 ',1),original+'1 1 2 3\n',
                         original.replace('coordinate complex general','coordinate complex hermitian')]:
                (path/'matrix.mtx').write_text(text)
                with self.assertRaises(ValueError):complex_rows(path/'matrix.mtx')
            (path/'matrix.mtx').write_text(original+'%changed hash\n')
            with self.assertRaises(ValueError):ComplexFixture(manifest,sha256(policy),'public-complex-v1')
            (path/'rhs.mtx').write_text('%%MatrixMarket matrix array complex general\n3 2\n1 NaN\n')
            with self.assertRaises(ValueError):complex_array(path/'rhs.mtx')


def driver_check(binary, backend):
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
                 MKL_THREADING_LAYER='GNU')
        for operation in ('N','T','H'):
            fixture(path,operation)
            target=path/f'candidate-{operation}.mtx'
            env['FORTSPARSE_BENCH_SOLUTION']=str(target)
            args=driver_arguments(backend,['complex-mm','1','1','rhs2',operation,'training',str(path)])
            result=subprocess.run([binary,*args],env=env,capture_output=True,text=True,timeout=.8)
            row=parse_csv(result.stdout)
            if row and row['status']=='backend_unavailable':raise SystemExit(77)
            assert result.returncode==0 and row and row['status']=='ok',result.stdout+result.stderr
            receipt=dict(status='ok',record=row,stdout=result.stdout)
            record_umfpack_configuration(receipt,backend)
            assert receipt['status']=='ok',receipt
            candidate=complex_array(target)
            assert max(abs(a-b) for col,truth in zip(candidate,X) for a,b in zip(col,truth))<2e-12
            assert row['backward_error']<1e-12 and row['forward_error']<1e-10,row
        # Existing candidate files are never overwritten.
        previous=target.read_bytes()
        result=subprocess.run([binary,*args],env=env,capture_output=True,text=True,timeout=.8)
        assert result.returncode!=0 and target.read_bytes()==previous
        for damaged in ('matrix','rhs'):
            fixture(path)
            if damaged=='matrix':
                (path/'matrix.mtx').write_text('%%MatrixMarket matrix coordinate complex general\n3 3 1\n1 1 NaN 0\n')
            else:
                write_array(path/'rhs.mtx',[[1+0j,2+0j]])
            target=path/f'invalid-{damaged}-export.mtx'
            env['FORTSPARSE_BENCH_SOLUTION']=str(target)
            result=subprocess.run([binary,*args],env=env,capture_output=True,text=True,timeout=.8)
            assert result.returncode!=0 and not target.exists()
    print('PASS complex-MM original full matrix, independent N/T/H and candidate export')


if __name__=='__main__':
    if len(sys.argv)==3:driver_check(sys.argv[1],sys.argv[2])
    else:unittest.main()
