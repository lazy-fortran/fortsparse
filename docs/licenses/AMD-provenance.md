# AMD ordering provenance

`src/native/fortsparse_amd_kernel.inc` is a Fortran port of SuiteSparse AMD
3.3.3 (`amd_aat.c`, `amd_1.c`, `amd_2.c`, `amd_postorder.c`,
`amd_post_tree.c`) from the
[SuiteSparse v7.8.3 release](https://github.com/DrTimothyAldenDavis/SuiteSparse/releases/tag/v7.8.3),
with the default controls and without the statistics and debug code. It
carries the AMD copyright and BSD 3-clause notice; the
[AMD license](AMD-LICENSE.txt) is retained without edits.
`test_native_order_amd_reference` checks that its permutations equal
`amd_order`'s.

`src/native/fortsparse_amd.f90` (MIT) wraps the port for one-based CSC input:
it returns `perm[new] = old` and `inverse[old] = new`, both one-based, and
checks the permutation before exposing it. Input validation matches the METIS
adapter's. No AMD or SuiteSparse_config library is linked.
