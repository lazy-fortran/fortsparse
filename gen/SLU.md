# Static-structure supernodal LU source

`gen_slu_leaves.f90` retains the scalar mathematics of the supernodal LU
(`docs/lu-supernodal.md`): the update `a-b*c`, the division `a/b` by a nonzero
pivot and the scaling product `a*b*c`. FortSym emits them as pure elemental
value leaves. No derivative leaves are emitted: derivatives differentiate the
solve equation, never the factorization.

`slu_schedule.f90.inc` retains the allocation-free execution schedules that
call only these leaves:

- `slu_block_factor`: restricted-pivot LU of a w-by-w diagonal block, packed as
  LAPACK getrf. Pivoting picks the first maximum among block rows k..w and
  exchanges whole block rows. A pivot with |p| < tau_abs becomes
  sign(p)*tau_abs with sign(0) = +. Perturbations and exchanges are counted
  exactly; min |u_kk| and max |u_ij| are accumulated.
- `slu_panel_factor`: one column panel k0..k0+kb-1 of the blocked form of
  `slu_block_factor` for wide diagonal blocks, with the same pivot,
  exchange and perturbation rule; the caller applies the BLAS-3 trailing
  update (`dtrsm`, `dgemm`).
- `slu_trsm_upper` (B := B U^-1) and `slu_trsm_unit_lower_t` (B := B L^-T): the
  panel triangular solves below the BLAS thresholds and on the strict no-BLAS
  path.
- `slu_schur_update`: the general rectangular update C -= left right^T for the
  no-BLAS path.
- `slu_apply_swaps` and the diagonal-block solves `slu_block_lower`,
  `slu_block_upper`, `slu_block_upper_t`, `slu_block_lower_t`.

Pivot selection, exchanges and the perturbation decision are comparisons and
data movement, not arithmetic leaves. Dummies are assumed-size, so callers may
pass a panel element by sequence association. The mapped Schur scatter and
the off-diagonal solve scatter/gather reuse the generated SPD leaves, which
are the same subtract-into-map arithmetic.

From the repository root, with a built FortSym at the retained revision:

```sh
bash gen/regenerate_slu.sh /path/to/fortsym [GEN_BUILD [FORTSYM_BUILD]]
git diff --exit-code -- src/generated/
```

Ordinary builds use the committed `src/generated/fortsparse_slu_leaves.f90` and
need neither FortSym nor generation. The complex twin (stage 8) is to come
from the same definitions. `test/solver/test_fortsparse_slu_leaves.f90` checks
each compiled kernel against an independent dense reference: textbook partial
pivoting restricted to the block, P A = L U reconstruction, matmul residuals,
and exact perturbation counts and signs.
