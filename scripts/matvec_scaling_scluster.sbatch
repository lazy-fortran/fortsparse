#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -o /home/ert/fortsparse-omp-%j.log
# Build fortsparse, run the full ctest, then csc_matvec strong scaling 1-16 threads.
# Usage: cd $HOME/fortsparse-omp-src; sbatch scripts/matvec_scaling_scluster.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR}"
E=$HOME/fs1/env
export PATH="$E/bin:$PATH" FC="$E/bin/x86_64-conda-linux-gnu-gfortran" CC="$E/bin/x86_64-conda-linux-gnu-gcc"
export LD_LIBRARY_PATH="$E/lib:${LD_LIBRARY_PATH:-}" OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=false
git rev-parse HEAD > revision.txt 2>/dev/null
cmake -S . -B build-sc -G Ninja -DCMAKE_BUILD_TYPE=Release -DFORTSPARSE_USE_SYSTEM_DEPS=ON \
    -DCMAKE_PREFIX_PATH="$E" -DCMAKE_Fortran_FLAGS="-O3 -march=native -ffast-math" > cfg.log 2>&1 || exit 3
cmake --build build-sc -j8 > build.log 2>&1 || exit 4
( cd build-sc && OMP_NUM_THREADS=4 ctest -j2 --output-on-failure --timeout 600 ) > ctest.log 2>&1
echo "ctest exit $?" >> ctest.log
: > scaling.log
for m in 100 200; do
  for t in 1 2 4 8 16; do
    OMP_NUM_THREADS=$t build-sc/bin/bench_fortsparse_matvec $m >> scaling.log 2>&1
  done
done
