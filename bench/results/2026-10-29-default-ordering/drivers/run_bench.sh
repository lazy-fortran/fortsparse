#!/bin/bash
# Fresh process per trial; fixed core; 1 thread; 16 GiB no-swap cap; 600 s timeout.
# Writes results.csv (renamed trials.csv in the results folder).
set -u
B=$(find build -name bench_default_ordering -type f -perm -u+x | head -1)
out=results.csv
echo "case,ordering,n,nnz,ordering_used,factor_nnz,analysis_s,stats_factor_s,factor_wall_s,solve_s,rel_residual,vmhwm_kib,process_cpu_s,active_threads,trial,exit,outer_wall_s,load1" > $out
run() {
  local kind=$1 ord=$2 m=$3 trial=$4
  local load=$(cut -d' ' -f1 /proc/loadavg)
  local t0=$(date +%s.%N)
  line=$(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -- taskset -c 2 timeout 600 $B $kind $ord $m 2>err_${kind}_${ord}_${trial}.txt)
  local code=$?
  local t1=$(date +%s.%N)
  [ -z "$line" ] && line="$kind,$ord,,,,,,,,,,,,"
  echo "$line,$trial,$code,$(echo "$t1 - $t0" | bc),$load" >> $out
}
for t in 1 2 3 4 5; do run laplace default 317 $t; run elastic default 317 $t; done
run laplace natural 317 1
run elastic natural 317 1
echo done >> $out
