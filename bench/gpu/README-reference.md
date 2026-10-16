# Paired frozen-binary reference runs

To compare a candidate against an earlier native executable under the same
host/GPU load window, add the frozen executable directory to the ordinary
runner. Reference copies enter the same shuffled block as candidate backends
for every case, host budget, workflow and fresh-process trial. They retain the
same requested backend name and gain `series="reference"`; ordinary entries
have `series="candidate"`. No new solver aliases or numerical semantics apply.

```sh
python3 bench/gpu/run_gpu_spd.py --bin-dir /candidate/gpu/bin \
  --reference-bin-dir /frozen/panel-stage/gpu/bin \
  --backends native-cuda,native-cuda-nd,cudss,cholmod \
  --cases grid2:9,grid3:8 --workloads warm,refactor --threads 1 --trials 5 \
  --device 0 --output /var/tmp/fortsparse-gpu/paired.json
```

`--reference-backends` defaults to `native-cuda,native-cuda-nd`; only names also
present in `--backends` get reference copies. It can include CPU native names;
the reference directory then supplies their CPU binary, while candidate CPU
runs continue using `--cpu-bin-dir`. Missing reference executables remain
unavailable records. External-only calibration still excludes native entries.

All executable paths have separate binary provenance and per-trial SHA256.
Candidate, optional CPU and reference build caches have separate metadata
keys. An archived cache may point at a source checkout that has since changed:
the observed source digest is explicitly labelled as such, not the historical
reference source identity. Retain the original build/source receipt alongside
immutable binaries; executable SHA256 identifies the measured code.

Analysis keeps candidate/reference series separate even if both hashes match,
so five trials per series never become ten pooled trials. Plots label each
series and short executable hash, use hollow reference markers, and separate
scaling lines by executable identity and series. Instrumented Nsight traces
remain separate diagnostic captures. Load changes across historical windows
cannot be resolved by comparing their independent medians.
