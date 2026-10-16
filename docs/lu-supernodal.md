# Static-structure supernodal LU for FortSparse: design

Status: accepted design (2026-10-16), implemented in stages; unchecked until each stage lands.
Base: `main` 16e80c3. Evidence: [static-structure diagnostic](../bench/results/2026-10-16-lu-static-structure/README.md).
Paths are repository-relative.

## 0. Why this design (evidence in one paragraph)

On both frozen KIN6D Jacobians, the pattern is 95.7% structurally symmetric. With a
maximum-product matching (MPM) plus its dual scaling and METIS ND on the matched A+Aᵀ,
the static symmetric L+U holds 0.97 M (n=5,525) and 4.94 M (n=19,989) entries. The
existing pivoted scalar LU with matching stores 4.23–4.45 M and 35.5–57.5 M entries.
Even with the weighted matching, it still stores 2.0 M and 20.6 M.

Dense emulation of exactly this numerical scheme was run on the scaled, matched,
ND-ordered matrices. Pivoting was restricted to supernode diagonal blocks, with
perturbation at 1e-13·‖Ã‖∞:

- **Perturbations:** zero.
- **Growth:** 131 (n=5,525) and 94 (n=19,989).
- **Backward error after one refinement step against the original A:**
  - n=5,525: normwise 1.1e-16, componentwise 6.5e-16.
  - n=19,989: normwise 1.9e-16, componentwise 5.2e-16.

The scheme only works with the weighted matching and scaling. The existing
Hopcroft–Karp matching leads to growth of up to 9e16 and 136–528 perturbations, and refinement never brings the componentwise backward error below 3.6e-4.

## 1. Inventory of reusable code (file:line)

Symbolic analysis
- `native_symbolic_counts(n,col_ptr,row_idx,parent,count,status)`, `src/native/fortsparse_native_symbolic.f90:36`.
  - Input is a full, sorted, symmetric one-based CSC. It builds the etree from the upper part (:72-86).
  - Output: `parent` and `count` (count excludes the diagonal).
  - Workspace: `native_symbolic_workspace_bytes` (:24).
  - The diagnostic checked it against an independent ancestor-walk count on 24 matrix/order pairs. Etrees and counts matched exactly.
- SPD analysis, `analyze` in `src/native/fortsparse_native.f90:564`.
  - Counts: :623.
  - Full `lp/li` construction: :654-685 (`row_reach` :833).
  - Blocked-path decision: :686-689.
  - There is no explicit postorder anywhere. AMD and METIS orders are already almost postordered: only 0–35 more exact supernodes form after an etree postorder (diagnostic).
- Supernode partition and relaxed amalgamation: `build_storage`, `src/native/fortsparse_native_supernodes.F90:232`.
  - Merge loop: :280-322.
  - Exact width ≤64, or ≤128 when nnz(L) > 128n.
  - Relaxed width ≤32 with padding ≤1.25.
  - Along an etree parent chain the rule depends only on `parent`/`count`, because m = width + count(last). The diagnostic replica relies on this.

Panel storage, schedule, OpenMP and memory
- Types in `src/native/fortsparse_native_supernodes.F90`:
  - `panel_t` (:34): `first,width,m,offset(c_size_t),rows(:),values(:,:)` pointer into one owning `host_values` (:76).
  - `supernodes_t` (:47-89).
- Left-looking, destination-driven numeric factorization: `factor` :488, `factor_destination` :621, `factor_panel` :639, `accumulate_update` :682.
  - Updates use relative row maps `rowmap(n,thread)`.
  - `dsyrk`/`dgemm` write directly into a contiguous destination, otherwise into the per-thread `update_work` followed by a generated mapped scatter (:698-761).
- Update lists: `setup_incoming` :982. Subtree task groups: `setup_groups` :1040, with `group_target_flops` = max(1e6, flops/(12T)).
- OpenMP task tokens with `depend(iterator)`: :548-590.
- Thread choice: :431-446. `use_blas` threshold: :428-430 (disabled under a strict budget).
- Ledger:
  - Charging: `charge` :122.
  - `fit_thread_scratch` :134 halves the thread count until scratch fits.
  - `build(...,budget)` :172 checks the ledger against `memory_bytes()` (:208-212).
  - Ledger type: `memory_budget_t` `reserve`/`release` in `src/fortsparse_memory.f90:329/:363`.
  - Plan types: `sparse_memory_plan_t` (:18-28), with `MEMORY_PLAN_EXACT_LAYOUT/BOUND/ESTIMATE` (:15-16).
  - Ordered SPD plan precedent: `native_plan_ordered_spd`, `fortsparse_native.f90:133`; public `sparse_analyze`, `src/fortsparse_solver.F90:262`.

Matching, ordering, LU and plumbing
- Structural matching: unweighted Hopcroft–Karp `structural_row_matching` (`src/native/fortsparse_native_matching.f90:16`) and `matching_graph_rows` (:154).
  - Called from `general_order` (`src/native/fortsparse_native_general.f90:218`; match at :252, graph at :266).
  - There is no weighted matching or scaling anywhere. The only scaling is the scalar LU's row max-abs `row_norm`.
- Orderings (perm[new]=old, one-based, 0 = OK):
  - `fsparse_native_order_general` (METIS ND on A+Aᵀ, `src/native/fsparse_native_order.c:141`).
  - `fsparse_native_order_amd_general` (`fsparse_native_order_amd.c:76`).
  - Bounded AMD with a caller arena: `fsparse_native_amd_bounded_size/_bounded` (`fsparse_native_amd_bounded.c:60/:77`; Fortran bindings in `fortsparse_native_amd.f90`).
  - ND is refused under a cap or LOW (`fortsparse_native_general.f90:381-385`).
- Scalar LU `lu_factor_t` (`src/native/fortsparse_native_lu.f90:31`) is **generated** by `gen/gen_lu.f90`.
  - Gilbert–Peierls with threshold row pivoting. Arguments include `row_preference` and `pruning`.
  - `solve(rhs,status,operation)` with `LU_NORMAL/TRANSPOSE/ADJOINT` (:21, :871).
  - This is the safety-net fallback.
- General adapter `native_general_backend_t` (`fortsparse_native_general.f90:48`):
  - `prepare` :~300-413; `general_factor_real` :416 (`mode`=0 real LU); `general_refactor_real` :491.
  - Retention: `retain_real` :1150 and `reserve_original` :1101.
  - `refine_real` :1235 runs ≤3 corrections to ≤32·eps backward error, otherwise SINGULAR.
  - `update_stats` :982; family strings :1039-1049.
- Residual oracle: `native_original_norms` and `native_residual_error` (`src/native/fortsparse_native_residual.f90:18/:49`), plus `backward_error` (:157).
- Public types (`src/fortsparse_backend.f90`):
  - `sparse_options_t` (:52-67): `ordering`, `structural_matching`, `pivot_tolerance`, `memory_policy`, `host_memory_limit_bytes`, `nthreads`.
  - `sparse_factor_stats_t` (:71-102): `native_factor_family`, `refinement_steps`, `backward_error`, `budget_*`. It has no pivot or perturbation counters.
- Selection: `ensure_backend` (`src/fortsparse_solver.F90:872`, :922-931) sends non-SPD native requests to `native_general_backend_t`. Matching is gated to native CPU GENERAL (:893-903).
- Derivatives are generic: `sparse_solve_jvp_real` :496 and `sparse_solve_vjp_real` :539 (transpose solve on the primal factors).
- Generated leaves:
  - `gen/gen_spd_leaves.f90` → `src/generated/fortsparse_spd_leaves.f90`: `spd_panel_factor` :397 (generated scalar Cholesky, not dpotrf), `spd_panel_solve` :426 (dtrsm above n≥96, m≥32), `spd_schur_update` :465, `spd_scatter_mapped_update` :510.
  - `gen/gen_lu.f90` → `src/generated/fortsparse_lu_leaves.f90` plus the scalar LU module. Regeneration: `gen/regenerate_*.sh` (FortSym pin in the scripts).
- Tests to keep green or extend:
  - `test_fortsparse_native_supernodes`, `_supernode_memory`, `test_spd_panel_no_blas`, `test_spd_amd_bounded`, `test_native_symbolic*`.
  - `test_fortsparse_native_lu*`, `test_lu_row_preference`, `test_fortsparse_structural_matching`, `test_matching_allocations`.
  - `test_fortsparse_native_general_api`, `_native_refinement`, `test_fortsparse_solve_derivatives`.

## 2. Pipeline

Ã = Q̃ᵀ · Dr · Pm · A · Dc · Q̃, with factors Ã = L̂ U.

1. **MPM + scaling.** Pm places the matched entries on the diagonal. Dr = exp(u) and Dc = exp(v)/colmax come from the optimal duals. Then |Ã_ij| ≤ 1 and diag(Ã) = 1.
2. **Ordering.**
   - Q from ND (METIS) or AMD on the pattern of PmA + (PmA)ᵀ. Reuse `matching_graph_rows` and the `*_general` wrappers unchanged.
   - Q̃ = Q followed by an etree postorder. Put the child with the largest count last; this is cheap and makes the chain rule robust.
3. **Static symbolic analysis.** Counts, then the partition (SPD rule), then per-supernode `rows(:)`, `incoming` lists, groups and the exact plan. No numerical values are touched.
4. **Numeric factorization.** Left-looking supernodal (§4) with pivoting restricted to each supernode's w×w diagonal block, plus perturbation.
5. **Solve and refinement.** Solve with the permutations and scalings (§5). Refine against the ORIGINAL A, and gate on independent backward error. On failure, fall back automatically to the existing scalar threshold LU and report it (§6).

## 3. Modules and types (real first)

- `src/native/fortsparse_native_mpm.f90`: `module fortsparse_native_mpm`, reentrant, with caller-owned workspace (same style as `fortsparse_native_matching`).
  - `mpm_workspace_bytes(n, nnz, cache_costs) -> integer(i8)`
  - `max_product_matching(n, ptr, rows, values, cmatch, u, v, iwork, rwork, status, stats)`, overloaded for real and complex (complex uses |a|).
  - `mpm_scaling(n, ptr, rows, values, cmatch, u, v, dr, dc)` (pure).
  - `stats`: `greedy_matched`, `augmentations`, `heap_pops`, `relaxations`, `log_range_dr/dc`.
- `src/native/fortsparse_native_supernode_schedule.F90` is extracted from `supernodes_t` without numeric changes:
  - `supernode_partition(n,parent,count,strict,starts,ends,m)`: the merge rule as a pure function of parent/count.
  - `supernode_rows(...)`.
  - `schedule_t` holding `colnode`, `incoming_*`, `group_*`, `tokens`, `rowmap`, `update_work`, thread fit and ledger charging.
  - `supernodes_t` (SPD) is re-based on `schedule_t`. Both families then share the same partition, maps and OpenMP groups.
  - Landed (stage 2) as `supernode_schedule_t`, composed privately in `supernodes_t`. The structural record `supernode_t(first,width,m,rows)` is the parent of the SPD `panel_t`. `supernode_gather_rows` forms sorted row unions from a column structure. The mapped scatter is re-exported as `supernode_scatter_mapped_update` (a use-rename of the generated SPD leaf, so nothing is regenerated). Test: `test_fortsparse_supernode_partition`.
- `src/native/fortsparse_native_slu.F90`: `module fortsparse_native_slu`, type `slu_factor_t`.
  - Public fields: `n, node_count, factor_nnz (padded L+U), symbolic_nnz, flops, perturbed_pivots, block_swaps, active_threads, peak_memory_bytes, failed_required_bytes, failed_available_bytes, max_growth, min_pivot`.
  - Private fields: `schedule_t`, `slu_panel_t(:)` (`first,width,m,offset_l,offset_u,rows(:)`, `ipiv(:)`), `host_values(:)` (one owning allocation), `cmatch(n)`, `dr(n)`, `dc(n)`, `perm(n)`, `inverse(n)`, `anorm_scaled`, `tau`.
  - Methods: `analyze(a, options, budget)`, `factor(a)`, `refactor(a)`, `solve(x, operation)`, `solve_block`, `free(budget)`, `memory_bytes()`, `plan(...) -> sparse_memory_plan_t`.
  - Landed (stage 3): `analyze(n, ptr, rows, values, status, ordering, low_memory, pivot_perturbation, budget)` on raw one-based CSC. It performs MPM, scaling, AMD (bounded caller arena when available, otherwise opaque and refused under a cap/LOW) or ND (opaque, refused under a cap/LOW), the symmetrized pattern, `native_symbolic_counts`, the largest-child-last postorder, `supernode_partition`, row unions, `supernode_schedule_t` lists/groups, and allocates every panel and scratch array. Phases MPM, graph, AMD arena, counts, rows, incoming, panels and scratch are charged before allocation; `failed_phase` names a refusal. `plan_bytes` (exact layout, computed before the panels are allocated) equals `memory_bytes()`.
  - Deviations: the retained 28n vectors are `drp, dcp` (scales in new labels) and `row_old, row_new, col_old`, which serve assembly and both solve directions without `cmatch`. `solve_work` holds n + max(m) reals (permuted RHS plus off-diagonal scratch), not 8n. Descriptors are charged as `storage_size(slu_panel_t)/8 + 8` per supernode (tokens and failure flags).
- A complex twin `slu_factor_z_t` comes later from the same generator definitions (§9).
- In the adapter, a new `mode = 5`, family `'LU_SUPERNODAL'`. `general_factor_real` dispatches on the options (§10).

## 4. Data layout and numeric scheme

Panels for supernode s, with m rows, w columns and `rows(1:w)` = first..first+w−1:
- **L panel** `L_s(m,w)`, column-major, lda=m. Rows 1:w hold the getrf-packed diagonal block (unit L below, U on and above the diagonal). Rows w+1:m hold L21.
- **Ut panel** `Ut_s(m−w,w)`, column-major, holds U12ᵀ: `Ut_s(r,c) = U(first+c−1, rows(w+r))`.
  - Because the structure is symmetric, U's row structure equals L's column structure. The same `rows(:)`, the same relative `rowmap`, the same `incoming` lists and the same GEMM shapes therefore serve both halves.
- **`ipiv(w)`** (int32) holds the local row interchanges of the diagonal block.
- **Exact payload per supernode:** `8·(m·w + (m−w)·w) + 4·m + 4·w`. The SPD panel is `8·m·w + 4·m`.

Update (left-looking, destination d, same `incoming` list as SPD):
- Rows of source s that fall in d's columns form one contiguous slice R_f..R_{f+nd−1} of `rows_s(w_s+1:)`. Rows beyond d's columns are R_{f+nd}..R_end.
- **GEMM1** updates d's L panel, including the FULL diagonal block:
  `L_d[map(R_f..R_end), map(R_f..R_{f+nd−1})] −= L_s[R_f..R_end, :] · Ut_s[R_f..R_{f+nd−1}, :]ᵀ`
- **GEMM2** updates d's Ut panel:
  `Ut_d[map(R_{f+nd}..R_end), map(R_f..R_{f+nd−1})] −= Ut_s[R_{f+nd}..R_end, :] · L_s[R_f..R_{f+nd−1}, :]ᵀ`
- When destination rows are contiguous, write directly with `dgemm` (as SPD :735-749). Otherwise compute into `update_work` and apply the mapped scatter. One per-thread buffer is reused sequentially for both GEMMs, so `maxwork` is the SPD value.

`factor_panel(d)`:
1. **Diagonal block.** Restricted-pivot LU of the w×w diagonal block, with the PARDISO rule.
   - In column k, take the row of maximum |·| among block rows k..w (partial pivoting inside the block; an optional diagonal preference uses a threshold τ_d).
   - If the chosen |p| < τ·‖Ã‖∞ (default τ = 1e-13), set p = sign(p)·τ‖Ã‖∞ (sign(0) = +) and increment `perturbed_pivots`.
   - Count swaps in `block_swaps`, and track max|U|/max|Ã| and min|u_kk|.
2. **Swaps.** Apply `ipiv` as column swaps of `Ut_d`, since d's rows are Ut's columns.
3. **Triangular solves.**
   - `L21 ← L21·U_dd⁻¹` with `dtrsm('R','U','N','N')`.
   - `Ut_d ← Ut_d·L_dd⁻ᵀ` with `dtrsm('R','L','T','U')`.

Why the swaps never need to touch descendant L entries in d's rows:
- The swaps are local to the block. Descendant L entries in d's rows stay keyed by d's pre-swap row positions.
- In the forward solve, contributions to block d accumulate in those labels before P_d is applied: L̂_dd = P_dᵀ L_dd. So no structural change or extra fill can occur.
- Row structures of different rows in one supernode are NOT equal in descendant columns, so LAPACK-style left swaps would be structurally wrong.
- The dense emulation applies full-row swaps. That is arithmetically equivalent: the same products are formed on the same values.

Flop shares from the diagnostic (relaxed, postordered, weighted+ND):
- Schur GEMMs: 77–83% of flops; 74–81% with k=w≥8.
- TRSM: 16–20%.
- Diagonal blocks: 1.2–2.2%.
- Supernodes with w≥8 carry 97–98% of LU flops, though up to half of the supernodes have width 1.

Landed (stage 4, serial): `factor(ptr, rows, values, status)` assembles M through binary searches in the static row lists (an entry outside the structure returns INVALID_MATRIX, so a changed pattern is re-analyzed) and runs the left-looking scheme above. Contiguous destinations use `dgemm` directly; otherwise `update_work` plus the mapped scatter. Rank-one sources scatter directly. Without BLAS (small flops or a strict cap), `slu_schur_update` and the generated triangular solves are used; BLAS `dtrsm` applies for w >= 96 and m-w >= 32, as for SPD. `perturbed_pivots`, `block_swaps`, `max_growth` (max|U|/max|M|) and `min_pivot` are exact for each factorization; `nonfinite` flags a non-finite panel entry (status SINGULAR). Nothing is allocated after `analyze`.

## 5. Solves, permutations, scalings

- **Normal solve** A x = b: t = Q̃ᵀ(Dr·(Pm b)); forward L̂ (per supernode: accumulate, apply P_d, unit-lower diagonal block); backward U; x = Dc·(Q̃ t).
- **Transpose** Aᵀy = c: t = Q̃ᵀ(Dc c); forward with Uᵀ; backward with L̂ᵀ, applying P_dᵀ after the diagonal-block Lᵀ solve. Then y = Pmᵀ(Dr·(Q̃ t)). ADJOINT (complex) conjugates the factors and scalings.
- The kernels are the SPD `solve`/`solve_block` loops (:765/:843) with two panels and `ipiv`. One `solve_work(n)` per call. Blocks of RHS reuse `solve_block`'s tiling.
- Retained vectors: `dr, dc` (real even for complex), `cmatch` (or its inverse), `perm`, `inverse` (int32). Total exactly 16n + 12n bytes.

## 6. Refinement, independent gate, fallback and reporting

- **Factor-time gates (cheap, deterministic).** Any of these sets `fallback_reason` and triggers a refactorization with the existing `lu_factor_t` (threshold 0.1, structural-matching preference plus AMD/ND, as today):
  - (a) non-finite panel entry;
  - (b) `perturbed_pivots > perturb_limit` (default: more than max(0, n/1000), measured on the corpus);
  - (c) `max_growth > 1/(64·√eps)`.
- **Mandatory validation.** Perturbed pivots > 0, or option `verify_factor` (default ON while opt-in), run one validation solve at factor time: b = A·x̂ for a fixed deterministic x̂. Refine with `refine_real` (:1235) against the retained ORIGINAL A, and compute `backward_error` (:157) independently.
  - If ≤3 steps do not reach the existing ≤32·eps gate, fall back (reason `REFINEMENT`).
  - Cost: two triangular solves per step plus one SpMV, i.e. O(nnz(L+U)).
- **Solve-time behaviour.** Refinement as today, plus the componentwise criterion of [native refinement](native-refinement.md): corrections continue while the componentwise backward error exceeds eps and halves, within the three-correction cap; acceptance stays the normwise gate. The validation solve above keeps the normwise-only loop. If the user RHS misses the gate and the original A is retained, one refactor-with-fallback is allowed, then the solve is repeated. Otherwise SINGULAR, as today.
- **Fallback under a cap.** The fallback must pass its own preflight. If it cannot, return RESOURCE_LIMIT with the scalar LU's `failed_required/available` and keep `fallback_reason`. Never silently exceed the cap.
- **Reporting** (new `sparse_factor_stats_t` fields; `src/fortsparse_backend.f90` is currently being edited by another worker, so add them in a separate commit):
  - `perturbed_pivots`, `block_pivot_swaps`, `pivot_growth`, `min_pivot`, `lu_method_used`, `lu_fallback` (logical), `lu_fallback_reason`;
  - `matching_used` (0 none, 1 structural, 2 weighted);
  - `scaling_log10_range`.
  - `native_factor_family` is `'LU_SUPERNODAL'` or `'LU_REAL'` after fallback.
  - Status stays OK after a successful fallback; the stats flag carries it.

Landed (stages 5 and 6, real): `slu_factor_t%gate_reason` (NONFINITE, PERTURBATION, GROWTH), adapter mode 6 (mode 5 was already complex LDLT) with the validation solve, the scalar fallback (reasons above plus REFINEMENT, SINGULAR for a structurally singular matching, SOLVE for the solve-time refactor) and the stats fields. Deviations: a perturbed factor retains the original matrix for refinement even when `refine` is false (otherwise its solves would be unrefined); if that retention does not fit the cap the factor falls back with reason PERTURBATION. The solve-time refactor-with-fallback runs only without a host cap, since the retained original moves through a temporary. Options: `lu_method` (AUTO=0 behaves as SCALAR, SCALAR=1 default, SUPERNODAL=2), `pivot_perturbation`, `verify_factor`. Tests: `test_fortsparse_slu_fallback`, `test_fortsparse_slu_api`, `test_slu_allocations` (Linux allocator oracle).

## 7. Exact memory plan (bytes; i8 checked arithmetic via `memory_checked_mul/add`)

With I = 4 and R = 8; ns = #supernodes; T = active threads.

- **MPM phase:** 7n·I (cmatch, rmatch, heap, hpos, touched, pred, state) + 4n·R (u, v, dist, colmax) + [cache_costs] nnz·R. Not retained beyond the dr/dc/cmatch output.
- **Graph phase:** matched graph nnz·I (`matching_graph_rows`), plus the ordering:
  - bounded AMD arena from `fsparse_native_amd_bounded_size`, or
  - METIS opaque (refused under a cap or LOW, as today).
- **Symbolic phase:** symmetrized pattern ≤ (2·nnz + n)·I + (n+1)·I, plus `native_symbolic_workspace_bytes(n)`, plus parent/count/post 3n·I.
- **Retained factor:**
  `P = Σ_s 8·(2·m_s·w_s − w_s²) + Σ_s 4·(m_s + w_s) + ns·(storage_size(slu_panel_t)/8 + 8) + 4n (colnode) + 12·#updates (incoming) + group arrays + T·(4n + 8·maxwork) + 8n (solve_work) + 28n (dr, dc, perm, inverse, cmatch)`
  - Every term is known after step 3 of §2, before any numeric allocation. Report it as `MEMORY_PLAN_EXACT_LAYOUT` (owned arrays only).
- **Refinement:** the existing `refinement_payload_bytes` plus the retained original (`reserve_original`).
- **Fallback:** the scalar LU preflight is separate and bounded by its own growth reserve. The plan reports it as `BOUND`, not exact.

Under a cap, charge each phase through `memory_budget_t%reserve`, before allocation:
- `fit_thread_scratch` reduces T;
- `use_blas` stays on only when not strict (BLAS internal buffers are unaccounted, same as SPD).

LOW_MEMORY:
- T=1;
- exact supernodes only (no relaxed padding; the diagnostic shows padding only adds 0.5–3.7%);
- MPM costs recomputed on the fly (no nnz·R cache);
- bounded AMD (ND refused);
- validation solve reuses the refinement buffers.

Accuracy does not change.

## 8. Refactorization and derivatives

- **Same pattern, new values.** Keep Pm, dr, dc, Q̃, the partition, the maps and the panel allocation, as PARDISO phase 22 does. Assemble, factor, and run the §6 gates.
  - On gate failure: redo the full analysis with a new MPM. If that fails too, fall back to scalar.
  - Report `last_reused_pattern` and `reused` as the adapter does now (`general_refactor_real` :491).
- **Derivatives.** Unchanged generic JVP/VJP (`fortsparse_solver.F90:496/:539`) use `solve(op=NORMAL/TRANSPOSE)` on the retained factors. The factor iterations are never differentiated. Sparse cotangents stay −b̄ᵢxⱼ on A's pattern. Refinement applies per operation.

## 9. OpenMP and generated kernels

- **OpenMP.** Reuse `setup_groups`/tokens/`depend(iterator)` unchanged through `schedule_t`. LU doubles the update flops, so `group_target_flops` uses LU flops. Record requested and actual active threads; per-thread `rowmap` and `update_work` are unchanged.
- **Determinism.** Task order changes the summation order of updates into a destination. Accept rounding-level differences, tested by original-matrix residual, not bitwise. Alternatively, group-ordered accumulation if bitwise reproducibility is later required.
- **FortSym/FortGen emission** (new `gen/gen_slu_leaves.f90` → `src/generated/fortsparse_slu_leaves.f90`; real and complex from one definition):
  1. `slu_block_factor` (w×w restricted-pivot LU with perturbation, counts, growth). Recommended as a **generated** kernel, not LAPACK `dgetrf`:
     - the perturbation and its statistics must sit inside the elimination loop, which `dgetrf` cannot do;
     - diagonal blocks are ≤128 wide and carry 1–3% of the flops;
     - it matches `spd_panel_factor`, which is also generated scalar code;
     - the complex twin comes from the same definition;
     - an exact, deterministic perturbation count is required for reporting.
     A blocked left-looking variant for w≥64 is a later, measured option only.
  2. `slu_panel_trsm_l/u` small-size fallbacks. BLAS `dtrsm` applies above the SPD-style thresholds.
  3. `slu_schur_update`, a general (not SYRK-lower) small update for the no-BLAS/strict path.
  4. Mapped scatter: reuse `spd_scatter_mapped_update` (identical subtract-into-map). Expose it under a neutral name when extracting the schedule.
  5. Supernodal solve leaves: diagonal-block forward/back with `ipiv`, transposed variants.

  Value-only kernels are needed; derivatives go through the solve equation. Independent compiled-kernel tests: dense LAPACK comparison per block; perturbation count and sign.
  - Landed (stage 4) as `gen/gen_slu_leaves.f90` + `gen/slu_schedule.f90.inc` → `src/generated/fortsparse_slu_leaves.f90` (`gen/SLU.md`, `gen/regenerate_slu.sh`): value leaves `slu_update`, `slu_divide`, `slu_scale`; schedules `slu_block_factor`, `slu_trsm_upper`, `slu_trsm_unit_lower_t`, `slu_schur_update`, `slu_apply_swaps` and the four diagonal-block solves. The off-diagonal solve scatter/gather and the mapped Schur scatter reuse the SPD leaves (identical arithmetic). Tests: `test_fortsparse_slu_leaves` (per kernel) and `test_fortsparse_slu_factor` (dgesv, residuals, emulation, transpose identity, block RHS, zero pivots, disconnected blocks, zero diagonals, BLAS and generated paths).

## 10. Options and AUTO dispatch

- **New option** `sparse_options_t%lu_method` = `FORTSPARSE_LU_AUTO` (0), `FORTSPARSE_LU_SCALAR` (1), `FORTSPARSE_LU_SUPERNODAL` (2).
  - Plus `pivot_perturbation` (real, default 1e-13).
  - The default stays SCALAR until consumer evidence exists.
  - `structural_matching` keeps its meaning for the scalar path. The supernodal path always uses MPM.
- **AUTO.** After the static symbolic analysis (ms; 1 s ≈ whole 19,989 probe on this Mac including I/O):
  - choose scalar Gilbert–Peierls when LU flops < F0 or the mean supernode width < W0;
  - otherwise choose supernodal.
  - F0 and W0 are calibrated on the frozen sub-second corpus with five fresh-process trials (measurement contract), then frozen.
  - A structurally singular MPM sends the request to scalar, which reports SINGULAR as now.

## 11. Weighted matching algorithm (MPM)

- **Problem.** Maximize ∏|a_{σ(j)j}|, which is equivalent to a min-cost perfect matching on c_ij = log colmax_j − log|a_ij| ≥ 0 over numerically nonzero entries.
- **Method.** Duff–Koster (2001) MC64 job-5 equivalent, independently implemented:
  1. Initial duals u_i = min_j c_ij and v_j = min_i(c_ij − u_i).
  2. Greedy matching on tight edges.
  3. For each free column, Dijkstra on reduced costs c_ij − u_i − v_j ≥ 0 over the alternating graph, using a binary heap with a position array and a touched list.
  4. Stop at the first free row (distance D). Update only the touched vertices: u_i −= D − d_i and v_{σ(i)} += D − d_i for finalized rows, and v_{j0} += D. Augment along `pred`.
- **Complexity.** Worst case O(n·(nnz + n) log n); typically far less.
  - KIN6D: 1,474 / 5,700 free after greedy.
  - 17,476 / 76,254 heap pops in total.
  - 7.6 ms / 31 ms unqualified on this Mac.
- **Memory.** Exactly 7n·I + 4n·R (+ nnz·R if costs are cached). No dynamic growth, so it is fully cap-accountable before the call.
- **Output certificate (tested).**
  - min reduced cost ≥ −tol;
  - matched reduced cost = 0;
  - Σc(matched) = Σu + Σv (LP duality ⇒ optimal).
  - The prototype reached 2.4e-15 / 4.9e-15 and equal objectives. On n=5,525 its matching is identical to scipy's.
- **Scaling** follows from the duals: |Dr A Dc| ≤ 1 with matched entries = 1. Verified max 1+2.7e-15. The log10 ranges of dr/dc are 4.6/7.6 (n=5,525) and 4.8/9.1 (n=19,989), so the scaling is essential.
- **Note.** scipy's `min_weight_full_bipartite_matching` took over 3.5 minutes on n=19,989 without finishing, so it is not a viable oracle at that size. Use it only for small cross-checks.

## 12. Staged implementation (small, separately testable commits)

1. **MPM module** (`fortsparse_native_mpm`, real + complex magnitudes).
   - Tests: brute-force optimality over all n! permutations for n≤8 (random sparse, ties, wide magnitude ranges); dual feasibility and complementary slackness; |DrADc|≤1 with diag=1; structurally singular → SINGULAR; numerically-zero entries excluded; non-finite → INVALID; workspace-size identity and no allocation inside (allocation test like `test_matching_allocations`).
2. **Schedule extraction** (`supernode_schedule`) from `supernodes_t`, with no numeric change.
   - Tests: existing SPD supernode, memory, no-BLAS and bounded tests unchanged. A new test checks that `supernode_partition(parent,count)` equals the `li`-scan partition on random SPD patterns, with and without a strict budget.
3. **LU static analysis + exact plan** (`slu_factor_t%analyze`, plan only).
   - Tests: nnz(L) against dense boolean elimination of PmA+(PmA)ᵀ for n≤40; postorder validity; plan bytes == `memory_bytes()` after allocation; cap exhaustion injected at each phase (MPM, graph, AMD arena, counts, rows, incoming, panels, scratch) → RESOURCE_LIMIT with required/available and zero residual ledger.
4. **Serial numeric factor + NORMAL/TRANSPOSE solve** (BLAS and generated paths; generated `slu_block_factor` leaf with independent compiled-kernel tests).
   - Tests: dense LAPACK `dgesv` comparison (n≤200, random unsymmetric with symmetric-ish patterns); known-solution original-matrix residual; disconnected blocks; matrices needing matching (zero diagonal); forced in-block swaps (compare with dense restricted emulation); exact-zero pivot → one perturbation with the expected sign and count; transpose identity yᵀ(Ax) = (Aᵀy)ᵀx; block RHS.
5. **Gates, refinement, fallback, stats** (separate stats commit coordinated with the `fortsparse_backend.f90` owner).
   - Tests: a constructed matrix whose static pivots need perturbation, with the gate passing after refinement; a constructed failure (growth or perturbation over the limit) → fallback to `lu_factor_t` with reason reported and the residual oracle passing; fallback refused under a cap → RESOURCE_LIMIT; non-finite input.
6. **Adapter/API integration** (`mode=5`, `lu_method`, refactor reuse, free/lifetime).
   - Tests: general API, refinement, lifetime and solve-derivative suites extended to SUPERNODAL (JVP/VJP against finite differences and adjoint identities); same-pattern refactor reuse plus re-analysis path.
7. **OpenMP subtree tasks.**
   - Tests: 1/2/4 threads agree with serial to a residual tolerance; actual active threads reported; budget-reduced threads.
8. **Complex generated twin** (`slu_factor_z_t`, ADJOINT).
   - Tests: complex dense references, Hermitian-adjoint identity, complex derivative tests.
9. **LOW_MEMORY** (exact supernodes, T=1, on-the-fly costs).
   - Tests: same accuracy as AUTO; lower plan; cap sweeps.
10. **AUTO calibration and frontier evidence** on the frozen sub-second corpus, then a KIN6D consumer draft PR with matched baseline/native CPU time and memory per AGENTS.md.

Later, not in the first plan:
- complete (row and column) supernode pivoting;
- the PARDISO option of keeping matching cycles inside one diagonal block;
- a 2×2 Bunch–Kaufman analogue for symmetric indefinite matrices;
- GPU panels.
