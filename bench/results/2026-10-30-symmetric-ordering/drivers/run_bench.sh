#!/bin/bash
# Usage: run_bench.sh BASE_BINARY NEW_BINARY OUT.csv
# Fresh process per trial; core 5; one OpenMP/BLAS thread; 16 GiB no-swap cap
# (local-cap.sh = systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0);
# 600 s timeout per run.
set -u
base=$1 new=$2 out=$3
cap=${LOCAL_CAP:-/home/ert/code/prompts/scripts/local-cap.sh}
echo "build,family,order,mode,m,n,shift,nnz,status,ordering_used,factor_nnz,analysis_s,stats_factor_s,factor_wall_s,backward_error,inertia_neg,inertia_zero,inertia_pos,exact_neg,vmhwm_kib,process_cpu_s,host_peak_bytes,plan_peak_bytes,active_threads,trial,exit,load1" > "$out"
run() {  # run TAG BINARY TRIALS FAMILY ORDER M SHIFT [MODE]
  local tag=$1 bin=$2 trials=$3; shift 3
  for t in $(seq 1 "$trials"); do
    local load; load=$(cut -d' ' -f1 /proc/loadavg)
    local line
    line=$(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$cap" taskset -c 5 timeout 600 "$bin" "$@" 2>/dev/null)
    local code=$?
    [ -z "$line" ] && line="$1,$2,${5:-auto},$3,,$4,,,,,,,,,,,,,,,,,"
    echo "$tag,$line,$t,$code,$load" >> "$out"
  done
}
for s in 0 0.3 1 2.01 3.5; do
  run base "$base" 5 ldlt natural 70 $s
  run base "$base" 5 ldlt amd 70 $s
  run new "$new" 5 ldlt default 70 $s
  run new "$new" 1 ldlt natural 70 $s
done
run base "$base" 5 kkt natural 40 0
run base "$base" 5 kkt amd 40 0
run new "$new" 5 kkt default 40 0
run base "$base" 5 kkt amd 70 0
run new "$new" 5 kkt default 70 0
run base "$base" 1 kkt natural 70 0
run base "$base" 5 hpd default 150 0
run new "$new" 5 hpd default 150 0
for s in 0.3 3.5; do
  run base "$base" 5 ldlh default 70 $s
  run new "$new" 5 ldlh default 70 $s
done
run base "$base" 5 zsym default 70 3.5
run new "$new" 5 zsym default 70 3.5
for mode in low refine cap; do
  run base "$base" 5 spd default 200 0 $mode
  run new "$new" 5 spd default 200 0 $mode
done
echo done >> "$out"
